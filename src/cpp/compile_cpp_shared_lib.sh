g++ -march=native -g -funroll-loops -L$fortran_dir -L$scs_libdir -L$lapack_dir -I$scs_incdir cpp_version.cpp -lfortran_funcs -lscsindir -llapack -shared -fPIC -o libcpp_version.so
