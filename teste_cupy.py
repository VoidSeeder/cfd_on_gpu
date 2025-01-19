import cupy as cp

try:
    # Criação de um array simples no CuPy
    x = cp.array([1, 2, 3, 4, 5])
    print(f"Array criado com CuPy: {x}")

    # Operação de soma usando a GPU
    y = x ** 2
    print(f"Resultado da operação na GPU: {y}")

    # Verificando a GPU disponível
    print(f"Dispositivo ativo: {cp.cuda.Device()}")
    

except cp.cuda.runtime.CUDARuntimeError as e:
    print("Erro ao acessar a GPU com CuPy:", e)
