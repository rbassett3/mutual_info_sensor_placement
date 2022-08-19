gfortran -Wall -funroll-loops -O3 -march=native -L$lapack_dir fortran_funcs.f90 -llapack -shared -fPIC -o libfortran_funcs.so
