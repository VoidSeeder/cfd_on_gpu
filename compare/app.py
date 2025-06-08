from gpu import gpu_script
from cpu import cpu_script
import gc
import cupy as cp
import time

for i in range(5):
  cpu_file_name = f"cpu_times_{i + 1}.txt"

  cpu_file = open(cpu_file_name, "xt")

  cpu_file.write("nVolY - Alocacao - Iteracoes\n")
  cpu_file.write("\n")

  cpu_file.close()

  print("Inicio com CPU")

  for nVolY in range(50, 610, 50):
    cpu_time = cpu_script(nVolY=nVolY)
    
    print(f"{nVolY} {cpu_time["allocation"]:.6f} {cpu_time["iteration"]:.6f}\n")
    
    cpu_file = open(cpu_file_name, "a")
    
    cpu_file.write(f"{nVolY} {cpu_time["allocation"]:.6f} {cpu_time["iteration"]:.6f}\n")
    
    cpu_file.close()
    
    gc.collect()

    time.sleep(0.1)

  gpu_file_name = f"gpu_times_{i + 1}.txt"
    
  gpu_file = open(gpu_file_name, "xt")

  gpu_file.write("nVolY - Alocacao - Iteraoees\n")
  gpu_file.write("\n")

  gpu_file.close()

  for nVolY in range(50, 610, 50):
    gpu_time = gpu_script(nVolY=nVolY)
    
    print(f"{nVolY} {gpu_time["allocation"]:.6f} {gpu_time["iteration"]:.6f}\n")
    
    gpu_file = open(gpu_file_name, "a")
    
    gpu_file.write(f"{nVolY} {gpu_time["allocation"]:.6f} {gpu_time["iteration"]:.6f}\n")
    
    gpu_file.close()
    
    cp.get_default_memory_pool().free_all_blocks()

    gc.collect()

    time.sleep(0.1)