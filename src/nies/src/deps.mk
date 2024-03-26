ies.o : definitions.o structure.o pot_mod.o dens_mod.o read_write.o maco_mod.o green_mod.o sgfva_mod.o 
definitions.o : 
inversion.o : definitions.o
matmul_blas.o : definitions.o
structure.o :definitions.o
quad3_mod.o :definitions.o
read_write.o :definitions.o
maco_mod.o : definitions.o
dens_mod.o : definitions.o
pot_mod.o :definitions.o structure.o
green_mod.o :definitions.o inversion.o matmul_blas.o structure.o
sgfva_mod.o : inversion.o green_mod.o
