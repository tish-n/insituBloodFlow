/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2015 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at 
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "palabos3D.h"
#include "palabos3D.hh"

#include "ibm3D.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "mpi.h"
#include "lammps.h"
#include "input.h"
#include "library.h"
#include "lammpsWrapper.h"

#include "latticeDecomposition.h"
// #include "nearestTwoNeighborLattices3D.h"

using namespace plb;
using namespace std;

typedef double T;
//#define DESCRIPTOR descriptors::ForcedN2D3Q19Descriptor
#define DESCRIPTOR descriptors::ForcedD3Q19Descriptor
//#define DYNAMICS BGKdynamics<T, DESCRIPTOR>(parameters.getOmega())
#define DYNAMICS GuoExternalForceBGKdynamics<T, DESCRIPTOR>(parameters.getOmega())

#define NMAX 150

const T pi = (T)4.*std::atan((T)1.);

static T poiseuillePressure(IncomprFlowParam<T> const &parameters, plint maxN)
{
    const T a = parameters.getNx()-1;
    const T b = parameters.getNy()-1;

    const T nu = parameters.getLatticeNu();
    const T uMax = parameters.getLatticeU();

    T sum = T();
    for (plint iN = 0; iN < maxN; iN += 2)
    {
        T twoNplusOne = (T)2*(T)iN+(T)1;
        sum += ((T)1 / (std::pow(twoNplusOne,(T)3)*std::cosh(twoNplusOne*pi*b/((T)2*a))));
    }
    for (plint iN = 1; iN < maxN; iN += 2)
    {
        T twoNplusOne = (T)2*(T)iN+(T)1;
        sum -= ((T)1 / (std::pow(twoNplusOne,(T)3)*std::cosh(twoNplusOne*pi*b/((T)2*a))));
    }

    T alpha = -(T)8 * uMax * pi * pi * pi / (a*a*(pi*pi*pi-(T)32*sum)); // alpha = -dp/dz / mu

    T deltaP = - (alpha * nu);

    return deltaP;
}

T poiseuilleVelocity(plint iX, plint iY, IncomprFlowParam<T> const& parameters, plint maxN)
{
    const T a = parameters.getNx()-1;
    const T b = parameters.getNy()-1;

    const T x = (T)iX - a / (T)2;
    const T y = (T)iY - b / (T)2;

    const T alpha = - poiseuillePressure(parameters,maxN) / parameters.getLatticeNu();

    T sum = T();

    for (plint iN = 0; iN < maxN; iN += 2)
    {
        T twoNplusOne = (T)2*(T)iN+(T)1;

        sum += (std::cos(twoNplusOne*pi*x/a)*std::cosh(twoNplusOne*pi*y/a)
             / ( std::pow(twoNplusOne,(T)3)*std::cosh(twoNplusOne*pi*b/((T)2*a)) ));
    }
    for (plint iN = 1; iN < maxN; iN += 2)
    {
        T twoNplusOne = (T)2*(T)iN+(T)1;

        sum -= (std::cos(twoNplusOne*pi*x/a)*std::cosh(twoNplusOne*pi*y/a)
             / ( std::pow(twoNplusOne,(T)3)*std::cosh(twoNplusOne*pi*b/((T)2*a)) ));
    }

    sum *= ((T)4 * alpha * a *a /std::pow(pi,(T)3));
    sum += (alpha / (T)2 * (x * x - a*a / (T)4));
    
    return sum;
}

T poiseuilleVelocity2D(plint iX, plint iY, IncomprFlowParam<T> const& parameters, plint maxN)
{
    const T a = parameters.getNx()-1;
    const T b = parameters.getNy()-1;

    const T x = (T)iX - a / (T)2;
    const T y = (T)iY - b / (T)2;

    const T uMax = parameters.getLatticeU();
    T sum = T();
    sum = uMax*(1-4*y*y/b/b); 
    return sum;
}

T poiseuilleVelocityHole2D(plint iX, plint iY, Array<T,3> center, plint r, Array<T,3> inletV)
{
    const T a = center[0]; // center of the inlet hole
    const T b = center[1];

    const T x = (T)iX - a ;
    const T y = (T)iY - b ;

    //const T uMax = parameters.getLatticeU();
    T sum = T();
    T R2 = x*x + y*y;
    if (R2 <= r*r){
        sum = inletV[2]*(1-R2/r/r); 
    }
    //pcout<<"velocity "<<sum<<endl;
    return sum;
}

template <typename T>
class SquarePoiseuilleVelocityHole {
public:
    SquarePoiseuilleVelocityHole(Array<T,3> center_, plint r_, Array<T,3> inletV_)
        : center(center_), r(r_), inletV(inletV_)
    { }
    void operator()(plint iX, plint iY, plint iZ, Array<T,3>& u) const  {
        u[0] = T();
        u[1] = T();
        u[2] = poiseuilleVelocityHole2D(iX, iY, center, r, inletV);
        //u[2] = poiseuilleVelocity2D(iX, iY, parameters, maxN);
    }
private:
    Array<T,3> inletV;
    Array<T,3> center;
    plint r;
};

template <typename T>
class SquarePoiseuilleDensityAndVelocity {
public:
    SquarePoiseuilleDensityAndVelocity(IncomprFlowParam<T> const& parameters_, plint maxN_)
        : parameters(parameters_),
          maxN(maxN_)
    { }
    void operator()(plint iX, plint iY, plint iZ, T &rho, Array<T,3>& u) const {
        rho = (T)1;
        u[0] = T();
        u[1] = T();
        u[2] = poiseuilleVelocity(iX, iY, parameters, maxN);
        //u[2] = poiseuilleVelocity2D(iX, iY, parameters, maxN);
    }
private:
    IncomprFlowParam<T> parameters;
    plint maxN;
};

template <typename T>
class SquarePoiseuilleVelocity {
public:
    SquarePoiseuilleVelocity(IncomprFlowParam<T> const& parameters_, plint maxN_)
        : parameters(parameters_),
          maxN(maxN_)
    { }
    void operator()(plint iX, plint iY, plint iZ, Array<T,3>& u) const  {
        u[0] = T();
        u[1] = T();
        u[2] = poiseuilleVelocity(iX, iY, parameters, maxN);
        //u[2] = poiseuilleVelocity2D(iX, iY, parameters, maxN);
    }
private:
    IncomprFlowParam<T> parameters;
    plint maxN;
};

template<typename T>
class WallDomain3D: public DomainFunctional3D{
  public:
    WallDomain3D(plint r_, plint r2_, plint n_):cx{230,125,125,187.5,125,125,125,125},cy{30,45,30, 30,15,30,15,10},
      cz{30,30,30,30,30,50,70,90},r(r_),r2(r2_),n(n_)
{
    /*plint cx[]={230,125,125,187.5,125,125,125,125};
    plint cy[]={30,45,30, 30,15,30,15,10};
    plint cz[]={30,30,30,30,30,50,70,90};*/
    }
  virtual bool operator ()(plint iX, plint iY, plint iZ) const {
    
    plint result = 1;
    if (iZ < 10 || iZ > 110) result=0;
    srand(time(NULL));
    int randnumber;
    randnumber = rand()%100;
    //if (randnumber > 50) result = 0;
    for (plint i=0;i<n;i++){
      T rSqr = (iX-cx[i])*(iX-cx[i]) + (iY-cy[i])*(iY-cy[i])+(iZ-cz[i])*(iZ-cz[i]);
      if (rSqr <= r*r)
        result = 0;
    }
    for (plint i=0;i<6;i++){
      T rSqr = util::sqr(iX-cx[i]) + util::sqr(iY-cy[i]);
      if (rSqr <= r2*r2)
        result = 0;
    }
    
   return result; 
  }
  virtual WallDomain3D<T> * clone() const {
    return new WallDomain3D<T>(*this);
  }
  private:
    plint cx[8],cy[8],cz[8];  
    plint r,r2,n;
};

Array<T,3> getVelocity(Array<T,3> targetValue, plint iT, plint periT)
{
    static T pi = std::acos((T)-1);

    return (targetValue * std::sin(2.0*pi*iT/periT));
}


void squarePoiseuilleSetup( MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
                            IncomprFlowParam<T> const& parameters,
                            OnLatticeBoundaryCondition3D<T,DESCRIPTOR>& boundaryCondition )
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    const plint nz = parameters.getNz();

    
    plint r(10),r2(6),n(8);

    Box3D top    = Box3D(0,    nx-1, ny-1, ny-1, 0, nz-1);
    Box3D bottom = Box3D(0,    nx-1, 0,    0,    0, nz-1);
    
    Box3D inlet    = Box3D(1,    nx-2, 1,    ny-2, 0,   0);
    //Box3D outlet = Box3D(0,    nx-1, 0,    ny-1, nz-1, nz-1);
    Box3D outlet = Box3D(1,    nx-2, 1,    ny-2, nz-1, nz-1);
    
    Box3D left   = Box3D(0,    0,    1,    ny-2, 0, nz-1);
    Box3D right  = Box3D(nx-1, nx-1, 1,    ny-2, 0, nz-1);
    
    //-- set velocity boundary--//
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, inlet );
    //boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, inlet, boundary::normalOutflow );
    //boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, outlet );
    //boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, outlet, boundary::normalOutflow );
    //boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, top, boundary::normalOutflow );
    //boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, bottom, boundary::normalOutflow );
    //boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, left, boundary::normalOutflow );
    //boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, right, boundary::normalOutflow );
    
    // set pressure boundary
    //boundaryCondition.setPressureConditionOnBlockBoundaries ( lattice, inlet );
    boundaryCondition.setPressureConditionOnBlockBoundaries ( lattice, outlet );
    boundaryCondition.setPressureConditionOnBlockBoundaries ( lattice, top );
    boundaryCondition.setPressureConditionOnBlockBoundaries ( lattice, bottom );
    boundaryCondition.setPressureConditionOnBlockBoundaries ( lattice, left );
    boundaryCondition.setPressureConditionOnBlockBoundaries ( lattice, right );

    //boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, top );
    //boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, bottom );
    
    //boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, left,boundary::outflow );
    //boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, right,boundary::outflow );
    //boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, left,boundary::freeslip );
    //boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, right,boundary::freeslip );
    //boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, left );
    //boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, right );
    
    //setBoundaryVelocity(lattice, inlet, SquarePoiseuilleVelocity<T>(parameters, NMAX));
    //setBoundaryVelocity(lattice, outlet, SquarePoiseuilleVelocity<T>(parameters, NMAX));
    
    //setBoundaryDensity(lattice, inlet, 1.004); // 10 mmHg over 1000 um->120 um
    setBoundaryDensity(lattice, outlet, 1./3.);
    setBoundaryDensity(lattice, top, 1./3.);
    setBoundaryDensity(lattice, bottom, 1./3.);
    setBoundaryDensity(lattice, left, 1./3.);
    setBoundaryDensity(lattice, right, 1./3.);
    
    //setBoundaryVelocity(lattice, inlet, Array<T,3>((T)0.0,(T)0.0,(T)0.0));
    /*setBoundaryVelocity(lattice, top, Array<T,3>((T)0.0,(T)0.0,(T)0.0));
    setBoundaryVelocity(lattice, bottom, Array<T,3>((T)0.0,(T)0.0,(T)0.0));
    setBoundaryVelocity(lattice, left, Array<T,3>((T)0.0,(T)0.0,(T)0.0));
    setBoundaryVelocity(lattice, right, Array<T,3>((T)0.0,(T)0.0,(T)0.0));
    setBoundaryVelocity(lattice, outlet, Array<T,3>((T)0.0,(T)0.0,(T)0.0));
*/
    //initializeAtEquilibrium(lattice, lattice.getBoundingBox(), SquarePoiseuilleDensityAndVelocity<T>(parameters, NMAX));
    //defineDynamics(lattice, lattice.getBoundingBox(), new WallDomain3D<T>(r,r2,n),new BounceBack<T,DESCRIPTOR>);
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), (T)1.0/3., Array<T,3>(0,0,0));

    lattice.initialize();
}

T computeRMSerror ( MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
                    IncomprFlowParam<T> const& parameters )
{
    MultiTensorField3D<T,3> analyticalVelocity(lattice);
    setToFunction( analyticalVelocity, analyticalVelocity.getBoundingBox(),
                   SquarePoiseuilleVelocity<T>(parameters, NMAX) );
    MultiTensorField3D<T,3> numericalVelocity(lattice);
    computeVelocity(lattice, numericalVelocity, lattice.getBoundingBox());

           // Divide by lattice velocity to normalize the error
    return 1./parameters.getLatticeU() *
           // Compute RMS difference between analytical and numerical solution
           std::sqrt( computeAverage( *computeNormSqr(
                          *subtract(analyticalVelocity, numericalVelocity)
                     ) ) );
}

void writeVTK(MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
              IncomprFlowParam<T> const& parameters, plint iter)
{
    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();
    VtkImageOutput3D<T> vtkOut(createFileName("vtk", iter, 6), dx);
    //vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", dx/dt);
    vtkOut.writeData<3,float>(*computeVelocity(lattice), "velocity", dx/dt);
    vtkOut.writeData<3,float>(*computeVorticity(*computeVelocity(lattice)), "vorticity", 1./dt);
    vtkOut.writeData<float>(*computeDensity(lattice), "density", 1.);
    //vtkOut.writeData<3,float>(*computeShearStress(lattice), "shearStress", 1);
}



int main(int argc, char* argv[]) {

    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");
/*
    if (argc != 2) {
        pcout << "Error the parameters are wrong. The structure must be :\n";
        pcout << "1 : N\n";
        exit(1);
    }*/

    //const plint N = atoi(argv[1]);
    const plint N = 1;// atoi(argv[1]);
    const T Re = 1e-2;
    const plint Nref = 50;
    const T uMaxRef = 0.01;
    const T uMax = 1.67e-3;//0.0000000225;//uMaxRef /(T)N * (T)Nref; // Needed to avoid compressibility errors.

    IncomprFlowParam<T> parameters(
            uMax,
            Re,
            N,
            60.,        // lx
            60.,        // ly
            80.         // lz
    );
    
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    const plint nz_s = 0;
    Box3D inlet    = Box3D(1,    nx-2, 1,    ny-2, nz_s,   nz_s);
    
    const T maxT     =400000;//6.6e4; //(T)0.01;
    plint iSave = 4000;//10;
    plint periT = 400000; // periodic time
    plint iniT = atoi(argv[2]);
    plint iCheck = 10*iSave;//10;
    plint r = 5;
    Array<T,3> center(31, 31, 1);
    writeLogFile(parameters, "3D square Poiseuille");

    LammpsWrapper wrapper(argv,global::mpi().getGlobalCommunicator());
    char * inlmp = argv[1];
    wrapper.execFile(inlmp);
   
    /*LammpsWrapper wrapper2(argv,global::mpi().getGlobalCommunicator());
    char * inlmp2 = argv[2];
    wrapper2.execFile(inlmp2);
    */

    //MultiTensorField3D<T,3> vel(parameters.getNx(),parameters.getNy(),parameters.getNz());
    plint Nx,Ny,Nz;
    Nx=parameters.getNx();Ny=parameters.getNy();Nz=parameters.getNz();
    pcout<<"Nx,Ny,Nz "<<parameters.getNx()<<" "<<parameters.getNy()<<" "<<parameters.getNz()<<endl;
    LatticeDecomposition lDec(parameters.getNx(),parameters.getNy(),parameters.getNz(),
                              wrapper.lmp);
    SparseBlockStructure3D blockStructure = lDec.getBlockDistribution();
    ExplicitThreadAttribution* threadAttribution = lDec.getThreadAttribution();
    plint envelopeWidth = 3;

    MultiBlockLattice3D<T, DESCRIPTOR> 
      lattice (MultiBlockManagement3D (blockStructure, threadAttribution, envelopeWidth ),
               defaultMultiBlockPolicy3D().getBlockCommunicator(),
               defaultMultiBlockPolicy3D().getCombinedStatistics(),
               defaultMultiBlockPolicy3D().getMultiCellAccess<T,DESCRIPTOR>(),
               new DYNAMICS );
    
    /*MultiScalarField3D<int> 
      geometry (MultiBlockManagement3D (blockStructure, threadAttribution, envelopeWidth ),
               defaultMultiBlockPolicy3D().getBlockCommunicator(),
               defaultMultiBlockPolicy3D().getCombinedStatistics(),
               defaultMultiBlockPolicy3D().getMultiScalarAccess<int>(),
               0);*/ // not working, not sure the reason, 4/20/2017
    
 
    //plint r(10),r2(6),n(8);
    //MultiScalarField3D<int> geometry(Nx,Ny,Nz);
    /*MultiScalarField3D<int> geometry(lattice);
    plb_ifstream geometryFile("in.geometry");
    if (!geometryFile.is_open()){
      pcout<<"Error: could not open file"<<endl;
      return -1;
    }
    geometryFile >>geometry;

    defineDynamics(lattice, geometry, new BounceBack<T,DESCRIPTOR>, 1);*/
    //Cell<T,DESCRIPTOR> &cell = lattice.get(550,5500,550);
    pcout<<"dx "<<parameters.getDeltaX()<<" dt  "<<parameters.getDeltaT()<<" tau "<<parameters.getTau()<<endl;
    pcout<<getMultiBlockInfo(lattice)<<std::endl;

/*
    MultiBlockLattice3D<T, DESCRIPTOR> lattice (
        parameters.getNx(), parameters.getNy(), parameters.getNz(), 
        new DYNAMICS );*/

    // Use periodic boundary conditions.
    //lattice.periodicity().toggle(2,true);
    //lattice.periodicity().toggle(0,true);

    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
        = createLocalBoundaryCondition3D<T,DESCRIPTOR>();

    squarePoiseuilleSetup(lattice, parameters, *boundaryCondition);

    // Loop over main time iteration.
    //util::ValueTracer<T> converge(parameters.getLatticeU(),parameters.getResolution(),1.0e-3);
      //coupling between lammps and palabos
   // v_max=1/(4*vis)*dp/dx*R^2; 
    Array<T,3> force(0,0.,0);
    //force[2]=1.5e-7;// velocity profile die out, 1e-7
    setExternalVector(lattice,lattice.getBoundingBox(),DESCRIPTOR<T>::ExternalField::forceBeginsAt,force);

    if (iniT > 0){
      loadBinaryBlock(lattice,"checkpoint.dat");
    }else{
      /*for (plint iT=0;iT<4e3;iT++){
          lattice.collideAndStream();
      }*/
    }
    T timeduration = T();
    writeVTK(lattice, parameters, 0);
    global::timer("mainloop").start();
    for (plint iT=0; iT<maxT; ++iT) {
    //for (plint iT=0; iT<2; ++iT) {
        if (iT%iSave ==0){
            pcout<<"Saving VTK file..."<<endl;
            writeVTK(lattice, parameters, iT);
        }
        if (iT%iCheck ==0 && iT >0){
            pcout<<"Timestep "<<iT<<" Saving checkPoint file..."<<endl;
            saveBinaryBlock(lattice,"checkpoint.dat");
        }

        Array<T,3> targetV(0,0,uMax);
        Array<T,3> inletV=getVelocity(targetV, iT, periT);
        //pcout<<"veloicty "<<inletV[0]<<" "<<inletV[1]<<" "<<inletV[2]<<endl;
        setBoundaryVelocity(lattice, inlet, SquarePoiseuilleVelocityHole<T>(center,r,inletV));// velocity is the value for v, inlet is box3d
        //setBoundaryVelocity(lattice, inlet, inletV);// velocity is the value for v, inlet is box3d
        // lammps to calculate force
        wrapper.execCommand("run 1 pre no post no");
        //wrapper2.execCommand("run 1 pre no post no");
        // Clear and spread fluid force
        //Array<T,3> force(0,0.,0);
        setExternalVector(lattice,lattice.getBoundingBox(),DESCRIPTOR<T>::ExternalField::forceBeginsAt,force);
        spreadForce3D(lattice,wrapper);
        //spreadForce3D(lattice,wrapper2);
        // Lattice Boltzmann iteration step.
        //forceCoupling3D(lattice,wrapper);
        //forceCoupling3D(lattice,wrapper2);
        lattice.collideAndStream();
        // Interpolate and update solid position
        interpolateVelocity3D(lattice,wrapper);
        //interpolateVelocity3D(lattice,wrapper2);
    }

    timeduration = global::timer("mainloop").stop();
    pcout<<"total execution time "<<timeduration<<endl;
    delete boundaryCondition;
}
