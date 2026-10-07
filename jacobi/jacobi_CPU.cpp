#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <utility>
#include <vector>

// Opções de compilação, informadas pelo Makefile
#ifndef OPCOES
#define OPCOES "desconhecidas"
#endif

typedef std::vector<double> Vetor;

// Uma iteração do método. Cada vetor guarda a malha linha por linha: o volume
// (i, j) fica na posição i * nVolX + j
double iteracao(int nVolY, int nVolX, const Vetor& Ap, const Vetor& Aw, const Vetor& Ae, const Vetor& As,
                const Vetor& An, const Vetor& Bp, Vetor& phi_old, Vetor& phi_new, double norma_b) {
    auto p = [nVolX](int i, int j) { return (std::size_t) i * nVolX + j; };

    int i, j;

    // Atualiza as bordas Oeste
    j = 0;
    for (i = 1; i < nVolY - 1; i++) {
        phi_old[p(i, j)] = (
            - Ae[p(i, j)] * phi_old[p(i, j + 1)]
            - As[p(i, j)] * phi_old[p(i - 1, j)]
            - An[p(i, j)] * phi_old[p(i + 1, j)]
            + Bp[p(i, j)]
        ) / Ap[p(i, j)];
    }

    // Atualiza as bordas Leste
    j = nVolX - 1;
    for (i = 1; i < nVolY - 1; i++) {
        phi_old[p(i, j)] = (
            - Aw[p(i, j)] * phi_old[p(i, j - 1)]
            - As[p(i, j)] * phi_old[p(i - 1, j)]
            - An[p(i, j)] * phi_old[p(i + 1, j)]
            + Bp[p(i, j)]
        ) / Ap[p(i, j)];
    }

    // Atualiza as bordas Sul
    i = 0;
    for (j = 1; j < nVolX - 1; j++) {
        phi_old[p(i, j)] = (
            - Aw[p(i, j)] * phi_old[p(i, j - 1)]
            - Ae[p(i, j)] * phi_old[p(i, j + 1)]
            - An[p(i, j)] * phi_old[p(i + 1, j)]
            + Bp[p(i, j)]
        ) / Ap[p(i, j)];
    }

    // Atualiza as bordas Norte
    i = nVolY - 1;
    for (j = 1; j < nVolX - 1; j++) {
        phi_old[p(i, j)] = (
            - Aw[p(i, j)] * phi_old[p(i, j - 1)]
            - Ae[p(i, j)] * phi_old[p(i, j + 1)]
            - As[p(i, j)] * phi_old[p(i - 1, j)]
            + Bp[p(i, j)]
        ) / Ap[p(i, j)];
    }

    // As bordas (e os cantos) valem para os dois vetores
    for (i = 0; i < nVolY; i++) {
        phi_new[p(i, 0)] = phi_old[p(i, 0)];
        phi_new[p(i, nVolX - 1)] = phi_old[p(i, nVolX - 1)];
    }

    for (j = 0; j < nVolX; j++) {
        phi_new[p(0, j)] = phi_old[p(0, j)];
        phi_new[p(nVolY - 1, j)] = phi_old[p(nVolY - 1, j)];
    }

    // Volumes internos: lê só a iteração anterior (phi_old), o que mantém o
    // método como Jacobi. Com um vetor só, o laço seria Gauss-Seidel
    for (i = 1; i < nVolY - 1; i++) {
        for (j = 1; j < nVolX - 1; j++) {
            phi_new[p(i, j)] = (- Aw[p(i, j)] * phi_old[p(i, j - 1)]
                - Ae[p(i, j)] * phi_old[p(i, j + 1)]
                - As[p(i, j)] * phi_old[p(i - 1, j)]
                - An[p(i, j)] * phi_old[p(i + 1, j)]
                + Bp[p(i, j)]) / Ap[p(i, j)];
        }
    }

    // Resíduo real do sistema linear: ||b - A*phi|| / ||b||
    // (os cantos fantasmas não fazem parte do sistema)
    // O laço das colunas não testa em que face o volume está: as linhas das
    // faces Sul e Norte são tratadas à parte e, nas demais, as faces Oeste e
    // Leste ficam fora do laço. A ordem da soma não muda
    double soma = 0.0;

    for (i = 0; i < nVolY; i++) {
        if (i == 0 || i == nVolY - 1) {
            // Linhas das faces Sul e Norte, sem os cantos
            for (j = 1; j < nVolX - 1; j++) {
                double res = Bp[p(i, j)] - Ap[p(i, j)] * phi_new[p(i, j)];

                res -= Aw[p(i, j)] * phi_new[p(i, j - 1)];
                res -= Ae[p(i, j)] * phi_new[p(i, j + 1)];

                if (i > 0) {
                    res -= As[p(i, j)] * phi_new[p(i - 1, j)];
                }

                if (i < nVolY - 1) {
                    res -= An[p(i, j)] * phi_new[p(i + 1, j)];
                }

                soma += res * res;
            }

            continue;
        }

        // Face Oeste
        j = 0;
        double res = Bp[p(i, j)] - Ap[p(i, j)] * phi_new[p(i, j)];

        res -= Ae[p(i, j)] * phi_new[p(i, j + 1)];
        res -= As[p(i, j)] * phi_new[p(i - 1, j)];
        res -= An[p(i, j)] * phi_new[p(i + 1, j)];

        soma += res * res;

        // Volumes internos
        for (j = 1; j < nVolX - 1; j++) {
            res = Bp[p(i, j)] - Ap[p(i, j)] * phi_new[p(i, j)];

            res -= Aw[p(i, j)] * phi_new[p(i, j - 1)];
            res -= Ae[p(i, j)] * phi_new[p(i, j + 1)];
            res -= As[p(i, j)] * phi_new[p(i - 1, j)];
            res -= An[p(i, j)] * phi_new[p(i + 1, j)];

            soma += res * res;
        }

        // Face Leste
        j = nVolX - 1;
        res = Bp[p(i, j)] - Ap[p(i, j)] * phi_new[p(i, j)];

        res -= Aw[p(i, j)] * phi_new[p(i, j - 1)];
        res -= As[p(i, j)] * phi_new[p(i - 1, j)];
        res -= An[p(i, j)] * phi_new[p(i + 1, j)];

        soma += res * res;
    }

    return std::sqrt(soma) / norma_b;
}

double segundos(std::chrono::steady_clock::time_point inicio, std::chrono::steady_clock::time_point fim) {
    return std::chrono::duration<double>(fim - inicio).count();
}

// Uso: jacobi_CPU [nVolY] [arquivo onde gravar o campo final]
//      jacobi_CPU --compilacao
int main(int argc, char* argv[]) {
    if (argc > 1 && std::strcmp(argv[1], "--compilacao") == 0) {
        std::printf("g++ %s, %s\n", __VERSION__, OPCOES);

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

    // Solução inicial de phi. O laço atualiza um vetor a partir do outro,
    // então ela é necessária nos dois
    Vetor phi_new(tamanho, 0.0);
    Vetor phi_old(tamanho, 0.0);

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

    std::printf("=> Início das iterações\n");

    auto end_allocation_time = std::chrono::steady_clock::now();
    auto start_iteration_time = end_allocation_time;

    while (residuo_iteracao > residuo_final && numero_iteracao < numero_maximo_iteracao) {
        residuo_iteracao = iteracao(nVolY, nVolX, Ap, Aw, Ae, As, An, Bp, phi_old, phi_new, norma_b);
        std::swap(phi_old, phi_new);
        numero_iteracao += 1;
    }

    auto end_time = std::chrono::steady_clock::now();

    // Depois da troca, a solução da última iteração está em phi_old
    std::swap(phi_old, phi_new);

    const char* estado = residuo_iteracao <= residuo_final ? "convergiu" : "nao_convergiu";

    std::printf("=> Resultado: nVolY=%d alocacao=%.6f iteracao=%.6f numero_iteracao=%d residuo=%.6e estado=%s\n",
                nVolY - 2,
                segundos(start_allocation_time, end_allocation_time),
                segundos(start_iteration_time, end_time),
                numero_iteracao,
                residuo_iteracao,
                estado);

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
