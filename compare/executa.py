import argparse
import datetime
import importlib.metadata
import os
import platform
import subprocess
import sys
import time

# Executa o script de um método, em CPU ou GPU e em uma das abordagens, para
# valores crescentes de nVolY, com várias repetições por malha, e grava os
# tempos em um arquivo.
# Cada execução é um processo novo do script do método, que informa o próprio
# tempo na linha "=> Resultado: ..."
# As versões em C++ e em CUDA precisam estar compiladas (make, na pasta acima
# desta)
#
# Com --energia, cada execução passa pelo energia.py, que mede a energia do
# processo inteiro: a do processador sempre e a da placa de vídeo nas
# execuções de GPU. As leituras de cada execução ficam em uma pasta ao lado
# do arquivo de saída
#
# Antes de cada execução de GPU há uma espera (--espera, em segundos), para a
# placa voltar ao repouso depois da execução anterior. A espera fica fora do
# tempo e da energia medidos. As execuções de CPU não levam espera
#
# Com --iteracoes N, cada execução faz exatamente N iterações, sem parar pelo
# resíduo (o estado informado passa a ser "iteracoes_fixas")
#
# Com --omega W, só no SOR, o método usa o fator de relaxação W (sem a opção,
# o padrão do script, 1). Uma execução que diverge para cedo, com o estado
# "divergiu"
#
# Exemplo: python executa.py jacobi GPU cupy --inicio 20 --fim 100 --passo 20

metodos = ["jacobi", "gauss_seidel", "gauss_seidel_red_black", "successive_over_relaxation_red_black"]

parser = argparse.ArgumentParser()
parser.add_argument("metodo", choices=metodos)
parser.add_argument("plataforma", choices=["CPU", "GPU"])
parser.add_argument("abordagem", choices=["numpy", "cupy", "numba", "serial", "openmp", "cuda"])
parser.add_argument("--inicio", type=int, required=True, help="primeiro nVolY")
parser.add_argument("--fim", type=int, required=True, help="último nVolY")
parser.add_argument("--passo", type=int, required=True, help="incremento do nVolY")
parser.add_argument("--repeticoes", type=int, default=5)
parser.add_argument("--saida", default=".", help="pasta do arquivo de saída")
parser.add_argument("--energia", action="store_true", help="mede a energia de cada execução")
parser.add_argument("--espera", type=float, default=12.0, help="segundos de espera antes de cada execução de GPU")
parser.add_argument("--iteracoes", type=int, help="número fixo de iterações de cada execução (sem a opção, até convergir)")
parser.add_argument("--omega", type=float, help="fator de relaxação do SOR (sem a opção, 1)")
args = parser.parse_args()

sor = args.metodo == "successive_over_relaxation_red_black"

if args.omega is not None and not sor:
  parser.error("--omega só vale para o successive_over_relaxation_red_black")

# Só as execuções de GPU levam espera
espera = args.espera if args.plataforma == "GPU" else 0.0

raiz = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

if args.abordagem == "serial":
  # C++ serial: o arquivo não leva o nome da abordagem
  script = os.path.join(raiz, args.metodo, f"{args.metodo}_{args.plataforma}.cpp")
  executavel = os.path.join(raiz, "build", f"{args.metodo}_{args.plataforma}")
  comando = [executavel]
elif args.abordagem == "cuda":
  # CUDA C++: o arquivo também não leva o nome da abordagem
  script = os.path.join(raiz, args.metodo, f"{args.metodo}_{args.plataforma}.cu")
  executavel = os.path.join(raiz, "build", f"{args.metodo}_{args.plataforma}")
  comando = [executavel]
elif args.abordagem == "openmp":
  script = os.path.join(raiz, args.metodo, f"{args.metodo}_{args.plataforma}_{args.abordagem}.cpp")
  executavel = os.path.join(raiz, "build", f"{args.metodo}_{args.plataforma}_{args.abordagem}")
  comando = [executavel]
else:
  script = os.path.join(raiz, args.metodo, f"{args.metodo}_{args.plataforma}_{args.abordagem}.py")
  executavel = None
  comando = [sys.executable, script]

if not os.path.exists(script):
  parser.error(f"não existe o script {os.path.relpath(script, raiz)}")

# A compilação fica fora da medição: aqui só se confere se ela foi feita
if executavel and (not os.path.exists(executavel) or os.path.getmtime(executavel) < os.path.getmtime(script)):
  parser.error(f"{os.path.relpath(executavel, raiz)} não existe ou é mais antigo que {os.path.relpath(script, raiz)}; execute make")

def saida_de(comando):
  try:
    saida = subprocess.run(comando, capture_output=True, text=True, cwd=raiz).stdout

    # Em uma linha só, para caber no cabeçalho
    return ", ".join(linha.strip() for linha in saida.splitlines())
  except OSError:
    return "indisponivel"

inicio = datetime.datetime.now()

os.makedirs(args.saida, exist_ok=True)

file_name = os.path.join(args.saida, f"{args.metodo}_{args.plataforma}_{args.abordagem}_{inicio:%Y%m%d_%H%M%S}.txt")

if args.energia:
  # Pasta das leituras de energia, com o nome do arquivo de saída
  pasta_energia = file_name[:-len(".txt")] + "_energia"

  medidor = [sys.executable, os.path.join(raiz, "compare", "energia.py"), "--saida", pasta_energia]

  if args.plataforma == "GPU":
    medidor.append("--gpu")

  medidor.append("--")
else:
  medidor = []

# Opção repassada ao script do método
iteracoes = ["--iteracoes", str(args.iteracoes)] if args.iteracoes is not None else []

# repr() mantém todos os algarismos do valor pedido
omega = ["--omega", repr(args.omega)] if args.omega is not None else []

# Modo "x": uma execução nova nunca sobrescreve um arquivo existente
file = open(file_name, "xt")

file.write(f"# Data: {inicio:%Y-%m-%d %H:%M:%S}\n")
file.write(f"# Comando: {' '.join(sys.argv)}\n")
file.write(f"# Script: {os.path.relpath(script, raiz)}\n")
file.write(f"# Compilacao: {saida_de(comando + ['--compilacao']) if executavel else '-'}\n")
file.write(f"# Commit: {saida_de(['git', 'rev-parse', 'HEAD'])}\n")
file.write(f"# Alteracoes nao commitadas: {saida_de(['git', 'status', '--short']) or 'nenhuma'}\n")
file.write(f"# Sistema: {platform.platform()}\n")
file.write(f"# Python: {platform.python_version()}\n")
file.write(f"# Pacotes: {', '.join(sorted(f'{d.name}=={d.version}' for d in importlib.metadata.distributions()))}\n")
file.write(f"# GPU: {saida_de(['nvidia-smi', '--query-gpu=name,driver_version', '--format=csv,noheader'])}\n")

# Limite de threads das bibliotecas numéricas (vazio = padrão da biblioteca)
for variavel in ["OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"]:
  file.write(f"# {variavel}: {os.environ.get(variavel, '')}\n")

file.write(f"# Energia: {'sim' if args.energia else 'nao'}\n")
file.write(f"# Espera antes de cada execucao (s): {espera:g}\n")
file.write(f"# Iteracoes fixas: {args.iteracoes if args.iteracoes is not None else 'nao'}\n")
file.write(f"# Omega: {(repr(args.omega) if args.omega is not None else '1 (padrao)') if sor else '-'}\n")

file.write("\n")
file.write("nVolY - Repeticao - Alocacao - Iteracoes - Numero de iteracoes - Residuo - Estado - Aquecimento - Threads"
           + (" - Duracao do processo - Energia do processador (J) - Energia da placa (J)" if args.energia else "") + "\n")
file.write("\n")

file.close()

for nVolY in range(args.inicio, args.fim + 1, args.passo):
  for repeticao in range(1, args.repeticoes + 1):
    time.sleep(espera)

    execucao = subprocess.run(medidor + comando + [str(nVolY)] + iteracoes + omega, capture_output=True, text=True)

    resultado = [l for l in execucao.stdout.splitlines() if l.startswith("=> Resultado:")]
    energia = [l for l in execucao.stdout.splitlines() if l.startswith("=> Energia:")]

    if resultado and (energia or not args.energia):
      # "=> Resultado: nVolY=20 alocacao=0.01 ..." vira um dicionário
      campos = dict(campo.split("=") for campo in resultado[0].split()[2:])

      linha = f"{nVolY} {repeticao} {campos['alocacao']} {campos['iteracao']} {campos['numero_iteracao']} {campos['residuo']} {campos['estado']} {campos.get('aquecimento', '-')} {campos.get('threads', '-')}"

      if args.energia:
        campos = dict(campo.split("=") for campo in energia[0].split()[2:])

        linha += f" {campos['duracao']} {campos['processador']} {campos.get('placa', '-')}"

      linha += "\n"
    else:
      # Execução que termina sem resultado (por exemplo, falta de memória)
      linha = f"{nVolY} {repeticao} erro: {execucao.stderr.strip().splitlines()[-1:]}\n"

    print(linha, end="")

    file = open(file_name, "a")

    file.write(linha)

    file.close()

print(f"=> Arquivo de saída: {file_name}")
