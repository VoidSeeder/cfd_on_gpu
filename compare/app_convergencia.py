from gpu import gpu_script
from cpu import cpu_script
import gc
import cupy as cp
import time

# Estudo até a convergência: cada execução roda até o resíduo cair abaixo da
# tolerância; o limite de iterações é só uma barreira contra execução sem fim

# Malhas em que todas as execuções convergem em minutos (só pares: o Jacobi
# não converge com nVolY ímpar)
malhas_cpu = [20, 40, 60, 80, 100]
# Na GPU cabem também malhas maiores; na CPU cada execução levaria horas
malhas_gpu = [20, 40, 60, 80, 100, 150, 200]

def mede(nome, script, malhas, repeticao):
  file_name = f"{nome}_convergencia_{repeticao}.txt"

  file = open(file_name, "xt")

  file.write("nVolY - Alocacao - Iteracoes - Numero de iteracoes - Estado\n")
  file.write("\n")

  file.close()

  for nVolY in malhas:
    tempo = script(nVolY=nVolY)

    # Execução que bate no limite de iterações não tem tempo válido
    estado = "convergiu" if tempo["converged"] else "nao_convergiu"

    linha = f"{nVolY} {tempo['allocation']:.6f} {tempo['iteration']:.6f} {tempo['iterations']} {estado}\n"

    print(linha)

    file = open(file_name, "a")

    file.write(linha)

    file.close()

    cp.get_default_memory_pool().free_all_blocks()

    gc.collect()

    time.sleep(0.1)

# Rodada de aquecimento, descartada: a primeira chamada à GPU no processo
# inclui a compilação dos kernels do CuPy
gpu_script(nVolY=20)

for i in range(5):
  print("Inicio com CPU")

  mede("cpu", cpu_script, malhas_cpu, i + 1)

  print("Inicio com GPU")

  mede("gpu", gpu_script, malhas_gpu, i + 1)
