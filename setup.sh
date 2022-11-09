#!/bin/sh

echo "Starting Setup"
echo "Leaving insituBloodFlow"
cd ..
mv insituBloodFlow BloodFlow # preference step - changing name of insituBloodFlow to BloodFlow
mv InsituBloodFlow BloodFlow # preference step for parent fork
echo "Making Directories"
mkdir insituBloodFlow
mv BloodFlow insituBloodFlow # /insituBloodFlow/BloodFlow
cd insituBloodFlow
mkdir install
mkdir src
mkdir build 

INSTALLDIR="$PWD/install" # /insituBloodFlow/install
BUILDDIR="$PWD/build" # /insituBloodFlow/build
SRCDIR="$PWD/src" # /insituBloodFlow/src
BASEDIR=$PWD # /insituBloodFlow/
BFDIR="$PWD/BloodFlow" # /insituBloodFlow/src 

###########
# palabos #
###########

cd $SRCDIR
echo "downloading Palabos:"
git clone https://gitlab.com/unigespc/palabos.git

# echo "changing directory to src/Palabos"
# cd palabos
# echo "checking out the compatible version"
# git checkout e498e8ad7f24fd7ff87313670db7873703c1fd3f
# echo "changing directory back to src"
# cd ..

##########
# LAMMPS #
##########

echo "downloading LAMMPS:"
git clone https://github.com/lammps/lammps.git
# echo "changing directory to src/lammps"
# cd lammps
# echo "checking out the compatible version"
# git checkout e960674cea38515ae3749218c314a9e1a3c6c140
# cd ..
# echo "changing directory to insituBloodFlow/BloodFlow/rbc"
cd $BFDIR/rbc
echo "copying the necessary additions to lammps library:"
cp bond_wlc_pow.* $SRCDIR/lammps/src
cp angle_rbc.* $SRCDIR/lammps/src
cp dihedral_bend.* $SRCDIR/lammps/src
cp fix* $SRCDIR/lammps/src
echo "changing directory to src/lammps/src"
cd $SRCDIR/lammps/src
echo "adding MOLECULE package to LAMMPS:"
make yes-MOLECULE
echo "adding MC package to LAMMPS:"
make yes-MC
echo "compiling LAMMPS as a library:"
make mpi mode=lib -j 8
echo "done compiling LAMMPS"

############
# ParaView #
############

# echo "\nchanging directory to SRC"
# cd $SRCDIR 
# echo "\ncurrent directory is:" $PWD
# echo "\ninstalling ParaView:"
# git clone --recursive https://gitlab.kitware.com/paraview/paraview.git
# cd paraview
# git checkout v5.10.1
# git submodule update --init --recursive
# echo "\nsetting ParaView install options:"
# cmake -B $BUILDDIR/paraview -S $SRCDIR/paraview -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$INSTALLDIR/paraview -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc -DPARAVIEW_USE_PYTHON=ON -DPARAVIEW_USE_MPI=ON -DPARAVIEW_USE_QT=OFF -DVTK_SMP_IMPLEMENTATION_TYPE=TBB -DPARAVIEW_BUILD_EDITION=CATALYST_RENDERING -DVTK_USE_X=OFF -DVTK_OPENGL_HAS_OSMESA=ON -DOSMESA_INCLUDE_DIR=/usr/include/GL/ -DOSMESA_LIBRARY=/usr/lib/x86_64-linux-gnu/libOSMesa.so
# echo "\nchanging directory to build/paraview"
# cd $BUILDDIR/paraview
# echo "installing ParaView"
# make -j8
# make -j8 install
# echo "\ndone installing ParaView"

#########
# ADIOS #
#########

echo "\nchanging directory to src"
cd $SRCDIR
echo "\ndownloading ADIOS2:"
git clone https://github.com/ornladios/ADIOS2.git
cd ADIOS2
git checkout v2.7.1
echo "setting up ADIOS2 options:"
cmake -S $SRCDIR/ADIOS2 -B $BUILDDIR/ADIOS2 -DCMAKE_INSTALL_PREFIX=$INSTALLDIR/ADIOS2 -DADIOS2_USE_Fortran=OFF -DADIOS2_BUILD_EXAMPLES=OFF -DCMAKE_BUILD_TYPE=Release
echo "changing directory to build/ADIOS2:"
cd $BUILDDIR/ADIOS2
echo "installing ADIOS2:"
make -j8
make -j8 install
echo "done installing ADIOS2"

##########
# SENSEI #
##########

echo "changing directory to src"
cd $SRCDIR
echo "downloading SENSEI:"
git clone https://github.com/SENSEI-insitu/SENSEI.git
echo "changing to SENSEI directory"
cd SENSEI
echo "checking out a working version of SENSEI"
git checkout develop

cp -rf $BFDIR/examples/singleCell4.0 $SRCDIR/SENSEI/miniapps
cp $BFDIR/examples/singleCell4.0/SENSEICMakeLists/CMakeLists.txt $SRCDIR/SENSEI
cp $BFDIR/examples/singleCell4.0/miniappCMakeLists/CMakeLists.txt $SRCDIR/SENSEI/miniapps

cd $SRCDIR/SENSEI/miniapps

mv singleCell4.0 singleCell

PVDIR="/home/tishn/myFork/insituBloodFlow/install/paraview/lib/cmake/paraview-5.10"

cmake -S  $SRCDIR/SENSEI -B $BUILDDIR/sensei -DCMAKE_INSTALL_PREFIX=$INSTALLDIR/sensei -DParaView_DIR=$PVDIR -DENABLE_PYTHON=ON -DENABLE_CATALYST_PYTHON=ON -DENABLE_CATALYST=ON -DENABLE_VTK_IO=ON -DENABLE_LAMMPS=OFF -DENABLE_MANDELBROT=OFF -DENABLE_SINGLECELL=ON -DENABLE_OSCILLATORS=ON -DENABLE_ADIOS2=ON -DADIOS2_DIR=$INSTALLDIR/ADIOS2/lib/cmake/adios2 
echo "changing directory to build/SENSEI"
cd $BUILDDIR/sensei
echo "installing SENSEI"
make -j8

# cmake -S  $SRCDIR/SENSEI -B $BUILDDIR/sensei -DCMAKE_INSTALL_PREFIX=$INSTALLDIR/sensei -DParaView_DIR=$INSTALLDIR/paraview/lib/cmake/paraview-5.10 -DENABLE_PYTHON=ON -DENABLE_CATALYST_PYTHON=ON -DENABLE_CATALYST=ON -DENABLE_VTK_IO=ON -DENABLE_LAMMPS=OFF -DENABLE_MANDELBROT=OFF -DENABLE_SINGLECELL=ON -DENABLE_OSCILLATORS=ON -DENABLE_ADIOS2=ON -DADIOS2_DIR=$INSTALLDIR/ADIOS2/lib/cmake/adios2 
# echo "changing directory to build/SENSEI"
# cd $BUILDDIR/sensei
# echo "installing SENSEI"
# make -j8

# cmake -S  $SRCDIR/SENSEI -B $BUILDDIR/sensei -DCMAKE_INSTALL_PREFIX=$INSTALLDIR/sensei -DParaView_DIR=$PVIEWDIR -DENABLE_CATALYST=ON -DENABLE_VTK_IO=ON -DENABLE_LAMMPS=OFF -DENABLE_MANDELBROT=OFF -DENABLE_OSCILLATORS=OFF -DENABLE_ADIOS2=OFF -DADIOS2_DIR=$ADIOSDIR
# echo "changing directory to build/SENSEI"
# cd $BUILDDIR/sensei
# echo "installing SENSEI"
# make -j8 
# make -j8 install
echo "done installing SENSEI"

############################
# singleCell EXAMPLE setup #
############################

# echo "changing directory to BloodFlow/examples/singleCell"
# SINGCDIR="$BFDIR/examples/singleCell"
# cd $SINGCDIR
# echo "making build directory:"
# mkdir build
# echo "changing directory to base"
# cd $BASEDIR 
# cmake -S $SINGCDIR -B $SINGCDIR/build -DSENSEI_DIR=$INSTALLDIR/sensei/lib/cmake -DPALABOS_ROOT=$SRCDIR/palabos -DBLOODFLOW_ROOT=$BFDIR -DLAMMPS_DIR=$SRCDIR/lammps -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_STANDARD=11


# cd $SINGCDIR/build
# make -j8
# cd ..
# echo "running test simulation for single cell example:"
# mpirun -n 4 cellFlow in.lmp4cell 1000 10