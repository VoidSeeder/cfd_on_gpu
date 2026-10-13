import argparse
import ctypes
import datetime
import os
import subprocess
import sys
import threading
import time

# Executa um comando e mede a energia consumida enquanto ele roda: a do
# processador sempre e a da placa de vídeo com --gpu.
# O comando é um processo separado e não é alterado; este programa só lê os
# contadores e grava as leituras em um arquivo.
#
# Processador: contador de energia acumulada do pacote (RAPL), em
# /sys/class/powercap/intel-rapl:0/energy_uj. O usuário precisa de permissão
# de leitura nesse arquivo. É lido no início, a cada intervalo e no fim.
# Placa de vídeo: contador de energia acumulada da biblioteca da NVIDIA
# (libnvidia-ml). É lido só no início e no fim: ele perde energia quando é
# lido muitas vezes.
#
# Exemplo: python compare/energia.py --gpu -- build/jacobi_GPU 40

parser = argparse.ArgumentParser()
parser.add_argument("--intervalo", type=float, default=1.0, help="segundos entre as leituras")
parser.add_argument("--gpu", action="store_true", help="lê também a placa de vídeo")
parser.add_argument("--saida", default=".", help="pasta do arquivo de leituras")
parser.add_argument("comando", nargs="+", help="comando a executar, depois de --")
args = parser.parse_args()

pacote = "/sys/class/powercap/intel-rapl:0"

# O contador do processador volta a zero ao passar deste valor
with open(os.path.join(pacote, "max_energy_range_uj")) as file:
  faixa = int(file.read())

try:
  contador = open(os.path.join(pacote, "energy_uj"))
except PermissionError:
  sys.exit(f"sem permissão de leitura em {pacote}/energy_uj")

if args.gpu:
  nvml = ctypes.CDLL("libnvidia-ml.so.1")
  placa = ctypes.c_void_p()

  if nvml.nvmlInit_v2() != 0 or nvml.nvmlDeviceGetHandleByIndex_v2(0, ctypes.byref(placa)) != 0:
    sys.exit("não foi possível abrir a placa de vídeo pela libnvidia-ml")

  energia_placa = ctypes.c_ulonglong()

def leitura():
  # Tempo (s) e energia do processador (uJ)
  contador.seek(0)

  return [time.perf_counter(), int(contador.read())]

def leitura_placa():
  # Energia da placa (mJ)
  if nvml.nvmlDeviceGetTotalEnergyConsumption(placa, ctypes.byref(energia_placa)) != 0:
    sys.exit("a placa de vídeo não informou a energia acumulada")

  return energia_placa.value

inicio = datetime.datetime.now()

if args.gpu:
  placa_inicio = leitura_placa()

leituras = [leitura()]

processo = subprocess.Popen(args.comando)

terminou = threading.Event()

def leituras_no_intervalo():
  # Leituras em instantes fixos a partir do início, para o intervalo não
  # acumular atraso
  proxima = leituras[0][0] + args.intervalo

  while not terminou.wait(max(proxima - time.perf_counter(), 0)):
    leituras.append(leitura())

    proxima += args.intervalo

# As leituras do intervalo ficam em outra linha de execução, para esta ficar
# parada até o comando terminar e fazer a leitura do fim sem atraso (a espera
# com tempo limite confere o fim do comando só a cada 0,05 s)
intervalo = threading.Thread(target=leituras_no_intervalo)

intervalo.start()

processo.wait()

# Leitura do fim, logo depois de o comando terminar
fim = leitura()

terminou.set()

intervalo.join()

# Uma leitura do intervalo pode ter caído depois da do fim
leituras = [valores for valores in leituras if valores[0] < fim[0]] + [fim]

if args.gpu:
  placa_fim = leitura_placa()

  nvml.nvmlShutdown()

duracao = leituras[-1][0] - leituras[0][0]

# Soma das diferenças entre leituras seguidas, corrigindo a volta a zero
processador = sum((b[1] - a[1]) % faixa for a, b in zip(leituras, leituras[1:])) / 1e6

resumo = f"duracao={duracao:.6f} processador={processador:.3f}"

if args.gpu:
  resumo += f" placa={(placa_fim - placa_inicio) / 1e3:.3f}"

resumo += f" leituras={len(leituras)} retorno={processo.returncode}"

os.makedirs(args.saida, exist_ok=True)

file_name = os.path.join(args.saida, f"energia_{inicio:%Y%m%d_%H%M%S_%f}.txt")

# Modo "x": uma execução nova nunca sobrescreve um arquivo existente
file = open(file_name, "xt")

file.write(f"# Data: {inicio:%Y-%m-%d %H:%M:%S}\n")
file.write(f"# Comando: {' '.join(args.comando)}\n")
file.write(f"# Intervalo: {args.intervalo}\n")
file.write(f"# Faixa do contador do processador (uJ): {faixa}\n")

if args.gpu:
  file.write(f"# Contador da placa (mJ): inicio={placa_inicio} fim={placa_fim}\n")

file.write(f"# Energia (J): {resumo}\n")
file.write("\n")
file.write("Tempo (s) - Processador (uJ)\n")
file.write("\n")

for tempo, processador_uj in leituras:
  file.write(f"{tempo - leituras[0][0]:.6f} {processador_uj}\n")

file.close()

print(f"=> Energia: {resumo} arquivo={file_name}")

sys.exit(processo.returncode)
