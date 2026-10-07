#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include <cuda_runtime.h>

// Opções de compilação, informadas pelo Makefile
#ifndef OPCOES
#define OPCOES "desconhecidas"
#endif

typedef std::vector<double> Vetor;

// Kernels de uma iteração do método, executados na GPU. Os vizinhos de um
// volume são sempre da outra cor, então os volumes de uma mesma cor não
// dependem uns dos outros: cada thread da GPU calcula um volume da cor atual
// (ou, no resíduo, uma coluna) e grava no próprio phi_new. Cada vetor guarda a
// malha linha por linha: o volume (i, j) fica na posição i * nVolX + j

// Faces da cor atual (0 = vermelho, i + j par; 1 = preto, i + j ímpar). A
// thread k calcula as faces Oeste e Leste da linha k e as faces Sul e Norte da
// coluna k
__global__ void bordas(int nVolY, int nVolX, const double* Ap, const double* Aw, const double* Ae, const double* As,
                       const double* An, const double* Bp, double* phi_new, int cor) {
    auto p = [nVolX](int i, int j) { return (std::size_t) i * nVolX + j; };

    int k = blockIdx.x * blockDim.x + threadIdx.x;
    int i, j;

    if (k >= 1 && k < nVolY - 1) {
        // Atualiza as bordas Oeste
        i = k;
        j = 0;
        if ((i + j) % 2 == cor) {
            phi_new[p(i, j)] = (
                - Ae[p(i, j)] * phi_new[p(i, j + 1)]
                - As[p(i, j)] * phi_new[p(i - 1, j)]
                - An[p(i, j)] * phi_new[p(i + 1, j)]
                + Bp[p(i, j)]
            ) / Ap[p(i, j)];
        }

        // Atualiza as bordas Leste
        j = nVolX - 1;
        if ((i + j) % 2 == cor) {
            phi_new[p(i, j)] = (
                - Aw[p(i, j)] * phi_new[p(i, j - 1)]
                - As[p(i, j)] * phi_new[p(i - 1, j)]
                - An[p(i, j)] * phi_new[p(i + 1, j)]
                + Bp[p(i, j)]
            ) / Ap[p(i, j)];
        }
    }

    if (k >= 1 && k < nVolX - 1) {
        // Atualiza as bordas Sul
        i = 0;
        j = k;
        if ((i + j) % 2 == cor) {
            phi_new[p(i, j)] = (
                - Aw[p(i, j)] * phi_new[p(i, j - 1)]
                - Ae[p(i, j)] * phi_new[p(i, j + 1)]
                - An[p(i, j)] * phi_new[p(i + 1, j)]
                + Bp[p(i, j)]
            ) / Ap[p(i, j)];
        }

        // Atualiza as bordas Norte
        i = nVolY - 1;
        if ((i + j) % 2 == cor) {
            phi_new[p(i, j)] = (
                - Aw[p(i, j)] * phi_new[p(i, j - 1)]
                - Ae[p(i, j)] * phi_new[p(i, j + 1)]
                - As[p(i, j)] * phi_new[p(i - 1, j)]
                + Bp[p(i, j)]
            ) / Ap[p(i, j)];
        }
    }
}

// Volumes internos da cor atual: a thread (k, i) calcula o k-ésimo volume
// dessa cor na linha i, então nenhuma thread fica sem volume para calcular
__global__ void internos(int nVolY, int nVolX, const double* Ap, const double* Aw, const double* Ae, const double* As,
                         const double* An, const double* Bp, double* phi_new, int cor) {
    auto p = [nVolX](int i, int j) { return (std::size_t) i * nVolX + j; };

    int k = blockIdx.x * blockDim.x + threadIdx.x;
    int i = blockIdx.y * blockDim.y + threadIdx.y;
    int j = 2 - (cor + i) % 2 + 2 * k;

    if (i >= 1 && i < nVolY - 1 && j < nVolX - 1) {
        phi_new[p(i, j)] = (- Aw[p(i, j)] * phi_new[p(i, j - 1)]
            - Ae[p(i, j)] * phi_new[p(i, j + 1)]
            - As[p(i, j)] * phi_new[p(i - 1, j)]
            - An[p(i, j)] * phi_new[p(i + 1, j)]
            + Bp[p(i, j)]) / Ap[p(i, j)];
    }
}

// Resíduo real do sistema linear, b - A*phi: a thread j soma os quadrados do
// resíduo dos volumes da coluna j (os cantos fantasmas não fazem parte do
// sistema)
__global__ void residuo(int nVolY, int nVolX, const double* Ap, const double* Aw, const double* Ae, const double* As,
                        const double* An, const double* Bp, const double* phi_new, double* soma_coluna) {
    auto p = [nVolX](int i, int j) { return (std::size_t) i * nVolX + j; };

    int j = blockIdx.x * blockDim.x + threadIdx.x;

    if (j < nVolX) {
        double soma = 0.0;

        for (int i = 0; i < nVolY; i++) {
            if ((i == 0 || i == nVolY - 1) && (j == 0 || j == nVolX - 1)) {
                continue;
            }

            double res = Bp[p(i, j)] - Ap[p(i, j)] * phi_new[p(i, j)];

            if (j > 0) {
                res -= Aw[p(i, j)] * phi_new[p(i, j - 1)];
            }

            if (j < nVolX - 1) {
                res -= Ae[p(i, j)] * phi_new[p(i, j + 1)];
            }

            if (i > 0) {
                res -= As[p(i, j)] * phi_new[p(i - 1, j)];
            }

            if (i < nVolY - 1) {
                res -= An[p(i, j)] * phi_new[p(i + 1, j)];
            }

            soma += res * res;
        }

        soma_coluna[j] = soma;
    }
}

// Encerra a execução se uma chamada da CUDA falhar (por exemplo, por falta de
// memória na GPU)
void confere(cudaError_t erro, const char* onde) {
    if (erro != cudaSuccess) {
        std::fprintf(stderr, "erro da CUDA em %s: %s\n", onde, cudaGetErrorString(erro));

        std::exit(1);
    }
}

// Aloca um vetor na GPU e copia para ele um vetor da memória principal
double* para_gpu(const Vetor& vetor) {
    double* vetor_gpu;

    confere(cudaMalloc(&vetor_gpu, vetor.size() * sizeof(double)), "alocação");
    confere(cudaMemcpy(vetor_gpu, vetor.data(), vetor.size() * sizeof(double), cudaMemcpyHostToDevice), "cópia para a GPU");

    return vetor_gpu;
}

double segundos(std::chrono::steady_clock::time_point inicio, std::chrono::steady_clock::time_point fim) {
    return std::chrono::duration<double>(fim - inicio).count();
}

// Uso: gauss_seidel_red_black_GPU [nVolY] [arquivo onde gravar o campo final]
//      gauss_seidel_red_black_GPU --compilacao
int main(int argc, char* argv[]) {
    if (argc > 1 && std::strcmp(argv[1], "--compilacao") == 0) {
        std::printf("nvcc %d.%d.%d, g++ %s, %s\n", __CUDACC_VER_MAJOR__, __CUDACC_VER_MINOR__, __CUDACC_VER_BUILD__,
                    __VERSION__, OPCOES);

        return 0;
    }

    // Início da medição do tempo de alocação
    auto start_allocation_time = std::chrono::steady_clock::now();

    // Configuração inicial
    double L = 2.0;  // Comprimento do domínio [m]
    double D = 0.01;  // Diâmetro do canal [m]
    double rho = 1e3;  // Densidade da água [kg/m³]
    double mu = 8.9e-4;  // Viscosidade da água [kg/m·s]
    double Gamma = 0.61 / 4200;  // Coeficiente de condução [W/(m·K)]
    double Pin = 1.0;  // Pressão na entrada [Pa]
    double Pout = 0.0;  // Pressão na saída [Pa]
    double Tin = 25.0;  // Temperatura de entrada [°C]
    double Twall = 100.0;  // Temperatura nas paredes [°C]

    // Número de volumes na direção y (primeiro argumento da linha de comando)
    int nVolY = argc > 1 ? std::atoi(argv[1]) : 400;

    // Configuração da malha
    double R = D / 2;
    double dx = D / nVolY;

    int nVolX = (int) std::nearbyint(L / dx);
    L = nVolX * dx;

    nVolX += 2;
    nVolY += 2;

    // Centros dos volumes na direção y, igualmente espaçados de -R - dx/2 a
    // R + dx/2. A solução não depende da posição em x
    Vetor y(nVolY);
    double passo = ((R + dx / 2) - (-R - dx / 2)) / (nVolY - 1);

    for (int i = 0; i < nVolY; i++) {
        y[i] = i * passo + (-R - dx / 2);
    }

    y[nVolY - 1] = R + dx / 2;

    // Campo de velocidade: só depende de y, e é o mesmo nas faces Oeste e
    // Leste e no centro do volume. Não há velocidade na direção y
    double dPdx = (Pout - Pin) / L;

    Vetor u(nVolY);

    for (int i = 0; i < nVolY; i++) {
        u[i] = 1 / (4 * mu) * (-dPdx) * (R * R - y[i] * y[i]);
    }

    double vs = 0.0;
    double vn = 0.0;

    // Coeficientes da equação de temperatura
    std::size_t tamanho = (std::size_t) nVolY * nVolX;

    Vetor Ap(tamanho, 0.0);
    Vetor Aw(tamanho, 0.0);
    Vetor Ae(tamanho, 0.0);
    Vetor As(tamanho, 0.0);
    Vetor An(tamanho, 0.0);
    Vetor Bp(tamanho, 0.0);

    auto p = [nVolX](int i, int j) { return (std::size_t) i * nVolX + j; };

    // Faces Oeste
    for (int i = 0; i < nVolY; i++) {
        Ap[p(i, 0)] = 1;
        Aw[p(i, 0)] = 0;
        Ae[p(i, 0)] = 1;
        As[p(i, 0)] = 0;
        An[p(i, 0)] = 0;
        Bp[p(i, 0)] = 2 * Tin;
    }

    // Faces Leste
    for (int i = 0; i < nVolY; i++) {
        Ap[p(i, nVolX - 1)] = 1;
        Aw[p(i, nVolX - 1)] = -1;
        Ae[p(i, nVolX - 1)] = 0;
        As[p(i, nVolX - 1)] = 0;
        An[p(i, nVolX - 1)] = 0;
        Bp[p(i, nVolX - 1)] = 0;
    }

    // Faces Sul
    for (int j = 0; j < nVolX; j++) {
        Ap[p(0, j)] = 1;
        Aw[p(0, j)] = 0;
        Ae[p(0, j)] = 0;
        As[p(0, j)] = 0;
        An[p(0, j)] = 1;
        Bp[p(0, j)] = 2 * Twall;
    }

    // Faces Norte
    for (int j = 0; j < nVolX; j++) {
        Ap[p(nVolY - 1, j)] = 1;
        Aw[p(nVolY - 1, j)] = 0;
        Ae[p(nVolY - 1, j)] = 0;
        As[p(nVolY - 1, j)] = 1;
        An[p(nVolY - 1, j)] = 0;
        Bp[p(nVolY - 1, j)] = 2 * Twall;
    }

    // Volumes internos
    for (int i = 1; i < nVolY - 1; i++) {
        for (int j = 1; j < nVolX - 1; j++) {
            Ap[p(i, j)] = dx * rho * (- std::min(0.0, u[i]) + std::max(0.0, u[i]) - std::min(0.0, vs) + std::max(0.0, vn)) + 4 * Gamma;

            Aw[p(i, j)] = -dx * rho * std::max(0.0, u[i]) - Gamma;
            Ae[p(i, j)] = dx * rho * std::min(0.0, u[i]) - Gamma;
            As[p(i, j)] = -dx * rho * std::max(0.0, vs) - Gamma;
            An[p(i, j)] = dx * rho * std::min(0.0, vn) - Gamma;
            Bp[p(i, j)] = 0;
        }
    }

    // Solução inicial de phi
    Vetor phi_new(tamanho, 0.0);

    // Resolução do sistema linear
    double residuo_iteracao = 1;
    int numero_iteracao = 0;
    int numero_maximo_iteracao = 1000000;
    double residuo_final = 1e-10;

    // Norma do termo fonte (os cantos fantasmas não fazem parte do sistema)
    double norma_b = 0.0;

    for (int i = 0; i < nVolY; i++) {
        for (int j = 0; j < nVolX; j++) {
            if ((i == 0 || i == nVolY - 1) && (j == 0 || j == nVolX - 1)) {
                continue;
            }

            norma_b += Bp[p(i, j)] * Bp[p(i, j)];
        }
    }

    norma_b = std::sqrt(norma_b);

    // Os coeficientes são montados na memória principal e copiados para a GPU
    double* Ap_gpu = para_gpu(Ap);
    double* Aw_gpu = para_gpu(Aw);
    double* Ae_gpu = para_gpu(Ae);
    double* As_gpu = para_gpu(As);
    double* An_gpu = para_gpu(An);
    double* Bp_gpu = para_gpu(Bp);

    // As duas cores são atualizadas no mesmo vetor
    double* phi_new_gpu = para_gpu(phi_new);

    // Soma dos quadrados do resíduo de cada coluna, na GPU e na memória
    // principal
    Vetor soma_coluna(nVolX, 0.0);
    double* soma_coluna_gpu = para_gpu(soma_coluna);

    // Divisão das threads da GPU em blocos: uma dimensão para as bordas e o
    // resíduo, duas para os volumes internos (em cada linha, só metade dos
    // volumes é de cada cor)
    int threads_1d = 256;
    int blocos_1d = (std::max(nVolX, nVolY) + threads_1d - 1) / threads_1d;

    dim3 threads_2d(16, 16);
    dim3 blocos_2d((nVolX / 2 + threads_2d.x - 1) / threads_2d.x, (nVolY + threads_2d.y - 1) / threads_2d.y);

    std::printf("=> Início das iterações\n");

    // A primeira iteração carrega os kernels na GPU; ela é descartada
    bool aquecimento = true;

    confere(cudaDeviceSynchronize(), "sincronização");
    auto end_allocation_time = std::chrono::steady_clock::now();
    auto start_iteration_time = end_allocation_time;
    auto end_warm_up_time = end_allocation_time;

    while (residuo_iteracao > residuo_final && numero_iteracao < numero_maximo_iteracao) {
        // Volumes vermelhos (i + j par) e depois pretos (i + j ímpar)
        for (int cor = 0; cor < 2; cor++) {
            bordas<<<blocos_1d, threads_1d>>>(nVolY, nVolX, Ap_gpu, Aw_gpu, Ae_gpu, As_gpu, An_gpu, Bp_gpu, phi_new_gpu, cor);
            internos<<<blocos_2d, threads_2d>>>(nVolY, nVolX, Ap_gpu, Aw_gpu, Ae_gpu, As_gpu, An_gpu, Bp_gpu, phi_new_gpu, cor);
        }

        residuo<<<blocos_1d, threads_1d>>>(nVolY, nVolX, Ap_gpu, Aw_gpu, Ae_gpu, As_gpu, An_gpu, Bp_gpu, phi_new_gpu, soma_coluna_gpu);

        // ||b - A*phi|| / ||b||: as somas das colunas são copiadas da GPU e
        // somadas na CPU
        confere(cudaMemcpy(soma_coluna.data(), soma_coluna_gpu, nVolX * sizeof(double), cudaMemcpyDeviceToHost), "iteração");

        double soma = 0.0;

        for (int j = 0; j < nVolX; j++) {
            soma += soma_coluna[j];
        }

        residuo_iteracao = std::sqrt(soma) / norma_b;

        numero_iteracao += 1;

        if (aquecimento) {
            aquecimento = false;
            confere(cudaMemcpy(phi_new_gpu, phi_new.data(), tamanho * sizeof(double), cudaMemcpyHostToDevice), "cópia para a GPU");
            residuo_iteracao = 1;
            numero_iteracao = 0;
            confere(cudaDeviceSynchronize(), "sincronização");
            end_warm_up_time = std::chrono::steady_clock::now();
            start_iteration_time = end_warm_up_time;
        }
    }

    confere(cudaDeviceSynchronize(), "sincronização");
    auto end_time = std::chrono::steady_clock::now();

    confere(cudaMemcpy(phi_new.data(), phi_new_gpu, tamanho * sizeof(double), cudaMemcpyDeviceToHost), "cópia do campo final");

    const char* estado = residuo_iteracao <= residuo_final ? "convergiu" : "nao_convergiu";

    std::printf("=> Resultado: nVolY=%d alocacao=%.6f iteracao=%.6f numero_iteracao=%d residuo=%.6e estado=%s aquecimento=%.6f\n",
                nVolY - 2,
                segundos(start_allocation_time, end_allocation_time),
                segundos(start_iteration_time, end_time),
                numero_iteracao,
                residuo_iteracao,
                estado,
                segundos(end_allocation_time, end_warm_up_time));

    // Grava o campo final (nVolY linhas de nVolX números de 8 bytes)
    if (argc > 2) {
        std::FILE* arquivo = std::fopen(argv[2], "wb");

        if (arquivo == NULL || std::fwrite(phi_new.data(), sizeof(double), tamanho, arquivo) != tamanho) {
            std::fprintf(stderr, "não foi possível gravar o campo em %s\n", argv[2]);

            return 1;
        }

        std::fclose(arquivo);
    }

    return 0;
}
