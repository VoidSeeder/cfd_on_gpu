import argparse
import datetime
import importlib.metadata
import os
import platform
import subprocess
import sys

# Executa o script de um método, em CPU ou GPU, para valores crescentes de
# nVolY, com várias repetições por malha, e grava os tempos em um arquivo.
# Cada execução é um processo novo do script do método, que informa o próprio
# tempo na linha "=> Resultado: ..."
#
# Exemplo: python executa.py jacobi GPU --inicio 20 --fim 100 --passo 20

metodos = ["jacobi", "gauss_seidel_red_black", "successive_over_relaxation_red_black"]

parser = argparse.ArgumentParser()
parser.add_argument("metodo", choices=metodos)
parser.add_argument("plataforma", choices=["CPU", "GPU"])
parser.add_argument("--inicio", type=int, required=True, help="primeiro nVolY")
parser.add_argument("--fim", type=int, required=True, help="último nVolY")
parser.add_argument("--passo", type=int, required=True, help="incremento do nVolY")
parser.add_argument("--repeticoes", type=int, default=5)
parser.add_argument("--saida", default=".", help="pasta do arquivo de saída")
args = parser.parse_args()

raiz = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
script = os.path.join(raiz, args.metodo, f"{args.metodo}_{args.plataforma}.py")

def saida_de(comando):
  try:
    saida = subprocess.run(comando, capture_output=True, text=True, cwd=raiz).stdout

    # Em uma linha só, para caber no cabeçalho
    return ", ".join(linha.strip() for linha in saida.splitlines())
  except OSError:
    return "indisponivel"

inicio = datetime.datetime.now()

os.makedirs(args.saida, exist_ok=True)

file_name = os.path.join(args.saida, f"{args.metodo}_{args.plataforma}_{inicio:%Y%m%d_%H%M%S}.txt")

# Modo "x": uma execução nova nunca sobrescreve um arquivo existente
file = open(file_name, "xt")

file.write(f"# Data: {inicio:%Y-%m-%d %H:%M:%S}\n")
file.write(f"# Comando: {' '.join(sys.argv)}\n")
file.write(f"# Script: {os.path.relpath(script, raiz)}\n")
file.write(f"# Commit: {saida_de(['git', 'rev-parse', 'HEAD'])}\n")
file.write(f"# Alteracoes nao commitadas: {saida_de(['git', 'status', '--short']) or 'nenhuma'}\n")
file.write(f"# Sistema: {platform.platform()}\n")
file.write(f"# Python: {platform.python_version()}\n")
file.write(f"# Pacotes: {', '.join(sorted(f'{d.name}=={d.version}' for d in importlib.metadata.distributions()))}\n")
file.write(f"# GPU: {saida_de(['nvidia-smi', '--query-gpu=name,driver_version', '--format=csv,noheader'])}\n")

# Limite de threads das bibliotecas numéricas (vazio = padrão da biblioteca)
for variavel in ["OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"]:
  file.write(f"# {variavel}: {os.environ.get(variavel, '')}\n")

file.write("\n")
file.write("nVolY - Repeticao - Alocacao - Iteracoes - Numero de iteracoes - Residuo - Estado\n")
file.write("\n")

file.close()

for nVolY in range(args.inicio, args.fim + 1, args.passo):
  for repeticao in range(1, args.repeticoes + 1):
    execucao = subprocess.run([sys.executable, script, str(nVolY)], capture_output=True, text=True)

    resultado = [l for l in execucao.stdout.splitlines() if l.startswith("=> Resultado:")]

    if resultado:
      # "=> Resultado: nVolY=20 alocacao=0.01 ..." vira um dicionário
      campos = dict(campo.split("=") for campo in resultado[0].split()[2:])

      linha = f"{nVolY} {repeticao} {campos['alocacao']} {campos['iteracao']} {campos['numero_iteracao']} {campos['residuo']} {campos['estado']}\n"
    else:
      # Execução que termina sem resultado (por exemplo, falta de memória)
      linha = f"{nVolY} {repeticao} erro: {execucao.stderr.strip().splitlines()[-1:]}\n"

    print(linha, end="")

    file = open(file_name, "a")

    file.write(linha)

    file.close()

print(f"=> Arquivo de saída: {file_name}")
