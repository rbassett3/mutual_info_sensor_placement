
#USER SHOULD EDIT THIS LINE WITH THE DIRECTORY CONTAINING LAPACK LIBRARY
lapack_dir=/path/to/system/libpack/directory
#USER SHOULD NOT EDIT ANYTHING BELOW THIS LINE

#get current directory
cur_dir=$(pwd)
scs_dir=$cur_dir/src/scs

#clone and compile scs
cd src
git clone https://github.com/rbassett3/scs.git
#checkout the for of scs used for this project
cd $scs_dir
export LDFLAGS=-L$lapack_dir
make

#set lib and include directories
scs_libdir=$cur_dir/src/scs/out
scs_incdir=$cur_dir/src/scs/include

#append lapack directory to ld_library_path
ld_library_path=$ld_library_path:$lapack_dir

#compile the fortran source
fortran_dir=$cur_dir/src/fortran/
cd $fortran_dir
source compile_linked_lib.sh

#append fortran and scs directory to ld_library_path
ld_library_path=$ld_library_path:$fortran_dir:$scs_libdir

#compile the cpp souce
cpp_dir=$cur_dir/src/cpp/
cd $cpp_dir
source compile_cpp_shared_lib.sh

#append cpp directory to ld_library_path
ld_library_path=$ld_library_path:$cpp_dir
export ld_library_path
export LD_LIBRARY_PATH=$ld_library_path
export cpp_dir

#compile the Cython source
python_dir=$cur_dir/src/python
pip install --global-option=build_ext --global-option="-R$cpp_dir:$fortran_dir:$scs_libdir" $python_dir

#change back to original directory
cd $cur_dir
