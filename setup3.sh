#!/bin/sh

echo "Starting Setup"
echo "Leaving insituBloodFlow"
cd ..
# mv insituBloodFlow BloodFlow # preference step - changing name of insituBloodFlow to BloodFlow
# mv InsituBloodFlow BloodFlow # preference step for parent fork
# echo "Making Directories"
# mkdir insituBloodFlow
# mv BloodFlow insituBloodFlow # /insituBloodFlow/BloodFlow
# cd insituBloodFlow
# mkdir install
# mkdir src
# mkdir build 

INSTALLDIR="$PWD/install" # /insituBloodFlow/install
BUILDDIR="$PWD/build" # /insituBloodFlow/build
SRCDIR="$PWD/src" # /insituBloodFlow/src
BASEDIR=$PWD # /insituBloodFlow/
BFDIR="$PWD/BloodFlow" # /insituBloodFlow/src 
SINGCDIR="$BFDIR/examples/bidirectionalSingleCell"


##########
# SENSEI #
##########

# echo "changing directory to src"
# cd $SRCDIR
# echo "downloading SENSEI:"
# git clone https://github.com/SENSEI-insitu/SENSEI.git
# echo "changing to SENSEI directory"
# cd SENSEI
# echo "checking out a working version of SENSEI"
# git checkout develop

cd $SRCDIR/SENSEI/miniapps 
rm -rf singleCell

cd $BUILDDIR/SENSEI/miniapps 
rm -rf singleCell
rm -rf cellFlow

cd $SINGCDIR
rm cellFlow

cp -rf $BFDIR/examples/bidirectionalSingleCell $SRCDIR/SENSEI/miniapps
cp $BFDIR/examples/bidirectionalSingleCell/SENSEICMakeLists/CMakeLists.txt $SRCDIR/SENSEI
cp $BFDIR/examples/bidirectionalSingleCell/miniappCMakeLists/CMakeLists.txt $SRCDIR/SENSEI/miniapps

cd $SRCDIR/SENSEI/miniapps

mv bidirectionalSingleCell singleCell

# cp $BUILDDIR/sensei/miniapps/cellFlow $SINGCDIR

PVDIR="/home/tishn/myFork/insituBloodFlow/install/paraview/lib/cmake/paraview-5.10" #modify this.
# PVDIR="$INSTALLDIR/paraview/lib/cmake/paraview-5.10"
ADIOSDIR="/home/tishn/myFork/insituBloodFlow/install/ADIOS2/lib/cmake/adios2" #modify this.

cmake -S  $SRCDIR/SENSEI -B $BUILDDIR/sensei -DCMAKE_INSTALL_PREFIX=$INSTALLDIR/sensei -DParaView_DIR=$PVDIR -DENABLE_PYTHON=ON -DENABLE_CATALYST_PYTHON=ON -DENABLE_CATALYST=ON -DENABLE_VTK_IO=ON -DENABLE_LAMMPS=OFF -DENABLE_MANDELBROT=OFF -DENABLE_SINGLECELL=ON -DENABLE_OSCILLATORS=ON -DENABLE_ADIOS2=ON -DADIOS2_DIR=$ADIOSDIR 
echo "changing directory to build/SENSEI"
cd $BUILDDIR/sensei
echo "installing SENSEI"
make -j8
make -j8 install
echo "done installing SENSEI"

############################
# singleCell EXAMPLE setup #
############################

echo "changing directory to BloodFlow/examples/singleCell"
# SINGCDIR="$BFDIR/examples/bidirectionalSingleCell"
cp $BUILDDIR/sensei/miniapps/cellFlow $SINGCDIR
# cd $SINGCDIR
# mpirun -n 4 cellFlow in.lmp4cell 1000 10 $PWD/cellFlow.xml