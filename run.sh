export NETCDF_DIR="/opt/netcdf_v4.9.2_openmpi"
export LD_LIBRARY_PATH="$NETCDF_DIR/lib"
export CMAKE_MODULE_PATH="/mnt/z/GitHub/cug-hydro/VegGroundWater.f90/cmake"

# export FCMP="/usr/bin/gfortran -fopenmp"
# export FF="mpif90 -fopenmp"
# export FC="/usr/bin/gfortran "

cmake -DCMAKE_MODULE_PATH="$CMAKE_MODULE_PATH" ..
