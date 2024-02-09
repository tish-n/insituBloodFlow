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

#include "mpi.h"
#include "lammps.h"
#include "input.h"
#include "library.h"
#include "lammpsWrapper.h"

#include "latticeDecomposition.h"
//#include "nearestTwoNeighborLattices3D.h"

using namespace plb;
using namespace std;

typedef double T;
//#define DESCRIPTOR descriptors::ForcedN2D3Q19Descriptor
#define DESCRIPTOR descriptors::ForcedD3Q19Descriptor
//#define DYNAMICS BGKdynamics<T, DESCRIPTOR>(parameters.getOmega())
#define DYNAMICS GuoExternalForceBGKdynamics<T, DESCRIPTOR>(parameters.getOmega())

#define NMAX 150

const T pi = (T)4.*std::atan((T)1.);

static T poiseuillePressure(IncomprFlowParam<T> const &parameters, plint maxN){
    const T a = parameters.getNx()-1;
    const T b = parameters.getNy()-1;

    const T nu = parameters.getLatticeNu();
    const T uMax = parameters.getLatticeU();

    T sum = T();
    for (plint iN = 0; iN < maxN; iN += 2){
        T twoNplusOne = (T)2*(T)iN+(T)1;
        sum += ((T)1 / (std::pow(twoNplusOne,(T)3)*std::cosh(twoNplusOne*pi*b/((T)2*a))));
    }
    for (plint iN = 1; iN < maxN; iN += 2){
        T twoNplusOne = (T)2*(T)iN+(T)1;
        sum -= ((T)1 / (std::pow(twoNplusOne,(T)3)*std::cosh(twoNplusOne*pi*b/((T)2*a))));
    }

    T alpha = -(T)8 * uMax * pi * pi * pi / (a*a*(pi*pi*pi-(T)32*sum)); // alpha = -dp/dz / mu
    T deltaP = - (alpha * nu);
    return deltaP;
}

T poiseuilleVelocity(plint iX, plint iY, IncomprFlowParam<T> const& parameters, plint maxN){
    const T a = parameters.getNx()-1;
    const T b = parameters.getNy()-1;

    const T x = (T)iX - a / (T)2;
    const T y = (T)iY - b / (T)2;

    const T alpha = - poiseuillePressure(parameters,maxN) / parameters.getLatticeNu();

    T sum = T();

    for (plint iN = 0; iN < maxN; iN += 2){
        T twoNplusOne = (T)2*(T)iN+(T)1;
        sum += (std::cos(twoNplusOne*pi*x/a)*std::cosh(twoNplusOne*pi*y/a)
             / ( std::pow(twoNplusOne,(T)3)*std::cosh(twoNplusOne*pi*b/((T)2*a)) ));
    }
    for (plint iN = 1; iN < maxN; iN += 2){
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

T poiseuilleVelocityHole2D(plint iX, plint iY, std::vector<Array<T,3>> centers, std::vector<plint> rs, vector<Array<T,3>> inletV)
{
    T sum = T();
    for (int i = 0; i < centers.size(); i++){
        const T a = centers[i][0]; // center of the inlet hole
        const T b = centers[i][1];

        const T x = (T)iX - a ;
        const T y = (T)iY - b ;

        //const T uMax = parameters.getLatticeU();
        T R2 = x*x + y*y;
        if (R2 <= rs[i]*rs[i]){
            sum += inletV[i][2]*(1-R2/rs[i]/rs[i]); 
        }
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

// SquarePoiseuilleVelocityHole for multiple centers and radii passed as vectors
template <typename T>
class MultipleSquarePoiseuilleVelocityHoles {
public:
    MultipleSquarePoiseuilleVelocityHoles(std::vector<Array<T,3>> centers_, std::vector<plint> rs_, vector<Array<T,3>> inletV_)
        : centers(centers_), rs(rs_), inletV(inletV_)
    { 
        assert(centers.size() == rs.size()); // Ensure the sizes match
    }
    void operator()(plint iX, plint iY, plint iZ, Array<T,3>& u) const  {
        u[0] = T();
        u[1] = T();
        u[2] = poiseuilleVelocityHole2D(iX, iY, centers, rs, inletV);

    }
private:
    std::vector<Array<T,3>> inletV;
    std::vector<Array<T,3>> centers;
    std::vector<plint> rs;
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
        int distanceFromCenter = iX*iX+iY*iY;
        int thisRadius = parameters.getNx()*parameters.getNx()/4;
        if (distanceFromCenter < thisRadius){
            u[2] = poiseuilleVelocity(iX, iY, parameters, maxN);
        }
        else{
            u[2] = T();
        }
        // u[2] = poiseuilleVelocity(iX, iY, parameters, maxN);
    }
private:
    IncomprFlowParam<T> parameters;
    plint maxN;
};

template <typename T>
class ShearTopVelocity {
public:
    ShearTopVelocity(IncomprFlowParam<T> const& parameters_, plint maxN_)
        : parameters(parameters_),
          maxN(maxN_)
    { }
    void operator()(plint iX, plint iY, plint iZ, Array<T,3>& u) const  {
        u[0] = T(); 
        u[1] = T(); 
        u[2] = parameters.getLatticeU();
    }
private:
    IncomprFlowParam<T> parameters;
    plint maxN;
};

template <typename T>
class ShearBottomVelocity {
public:
    ShearBottomVelocity(IncomprFlowParam<T> const& parameters_, plint maxN_)
        : parameters(parameters_),
          maxN(maxN_)
    { }
    void operator()(plint iX, plint iY, plint iZ, Array<T,3>& u) const  {
        u[0] = T();
        u[1] = T();
        u[2] = T();
    }
private:
    IncomprFlowParam<T> parameters;
    plint maxN;
};

// Array<T,3> getVelocity(Array<T,3> targetValue, plint iT, plint periT)
// {
//     static T pi = std::acos((T)-1);

//     return (targetValue);
// }


void squarePoiseuilleSetup( MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
                            IncomprFlowParam<T> const& parameters,
                            OnLatticeBoundaryCondition3D<T,DESCRIPTOR>& boundaryCondition )
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    const plint nz = parameters.getNz();
    // y 
    Box3D top    = Box3D(0,    nx-1, ny-1, ny-1, 0, nz-1); 
    Box3D bottom = Box3D(0,    nx-1, 0,    0,    0, nz-1); 
    
    // z
    Box3D inlet  = Box3D(0,    nx-1, 1,    ny-2, 0,    0);
    Box3D outlet = Box3D(0,    nx-1, 1,    ny-2, nz-1, nz-1);
    
    //x
    Box3D left   = Box3D(0,    0,    1,    ny-2, 1, nz-2);
    Box3D right  = Box3D(nx-1, nx-1, 1,    ny-2, 1, nz-2);

    // shear flow top bottom surface

    // boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, inlet, boundary::outflow );
    // boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, outlet, boundary::outflow );

    // boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, top );
    // boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, bottom );
    
    // boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, left, boundary::outflow );
    // boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, right, boundary::outflow );
    
    // setBoundaryVelocity(lattice, top, ShearTopVelocity<T>(parameters,NMAX));
    // setBoundaryVelocity(lattice, bottom, ShearBottomVelocity<T>(parameters,NMAX));
    
    // boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, inlet, boundary::outflow );
    // boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, outlet, boundary::outflow );

    // channel flow
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, inlet);
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, outlet);

    // turning off velocity bcs for now 7/31/23
    // boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, top );
    // boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, bottom );
    // boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, left );
    // boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, right );
    
    // setBoundaryVelocity(lattice, inlet, SquarePoiseuilleVelocity<T>(parameters, NMAX));
    // setBoundaryVelocity(lattice, outlet, SquarePoiseuilleVelocity<T>(parameters, NMAX));
    
    // setBoundaryVelocity(lattice, inlet, Array<T,3>((T)0.0,(T)0.0,(T)0.05)); // .5 in z direction XXX added by Nazariy 7/19
    // setBoundaryVelocity(lattice, outlet, Array<T,3>((T)0.0,(T)0.0,(T)0.05)); // .5 in z direction XXX added by Nazariy 7/19

    // turning off velocity BCs for now 7/31/23
    // setBoundaryVelocity(lattice, top, Array<T,3>((T)0.0,(T)0.0,(T)0.0));
    // setBoundaryVelocity(lattice, bottom, Array<T,3>((T)0.0,(T)0.0,(T)0.0));
    // setBoundaryVelocity(lattice, left, Array<T,3>((T)0.0,(T)0.0,(T)0.0));
    // setBoundaryVelocity(lattice, right, Array<T,3>((T)0.0,(T)0.0,(T)0.0));
    

    // initializeAtEquilibrium(lattice, lattice.getBoundingBox(), SquarePoiseuilleDensityAndVelocity<T>(parameters, NMAX));
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(),(T)1.0, Array<T,3>(0.0,0.0,0.0));

    lattice.initialize();
}



/// This functional defines a data processor for the instantiation
///   of bounce-back nodes following the half-circle geometry.
//XXXX ADDED by Nazariy: 
template <typename T, template <typename U> class Descriptor>
class DynamicBoundaryFunctional : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    DynamicBoundaryFunctional(plint xc_, plint yc_, plint radius_, T rn_, plint it_, plint zl_) : xc(xc_), yc(yc_), radius(radius_), rn(rn_), it(it_), zl(zl_) { }
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
    {
        //BounceBackNodes<T> bbDomain(N, radius);
        Dot3D relativeOffset = lattice.getLocation();
        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                // for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    T dist = (iY + relativeOffset.y - yc)*(iY + relativeOffset.y - yc) + 
                             (iX + relativeOffset.x - xc)*(iX + relativeOffset.x - xc); //XXXX change later - Nazariy
                    T xscale = .5 * rn * it;  // scaling factor for x (how wide) replace with a slider later.
                    T yscale = .5 * radius * sin(.0314*it);  // scaling factor for y (how close to center) dependent on iteration
                    // the representation is y = e ^ (-x^2), but the example is set up in Z direction so here 
                    // the x stands for domain in Z i.e. values from Z axis and y stands for value being subtracted from the radius
                    int cntZ = (iZ+relativeOffset.z-(zl / 2))/xscale; // position of Z with resepect to center of Z        
                    T modRad = yscale * exp(-cntZ*cntZ) + 1; // radius modification parameter
                    if (dist > (radius - modRad)*(radius - modRad)) {
                        lattice.attributeDynamics(iX, iY, iZ, new BounceBack<T, DESCRIPTOR>);
                    }
                    else{
                        lattice.attributeDynamics(iX, iY, iZ, new GuoExternalForceBGKdynamics<T, DESCRIPTOR>(1./0.95));
                    }
                }
            }
        }

        // //BounceBackNodes<T> bbDomain(N, radius);
        // Dot3D relativeOffset = lattice.getLocation();
        // for (plint iX = domain.x0; iX <= domain.x1; ++iX)
        // {
        //     for (plint iY = domain.y0; iY <= domain.y1; ++iY)
        //     {
        //         // just define bouncebacks for first and last Z layers 
        //         // first Z layer 
        //         plint iZ = domain.z0;
        //         T dist = (iY + relativeOffset.y - yc)*(iY + relativeOffset.y - yc) + 
        //                     (iX + relativeOffset.x - xc)*(iX + relativeOffset.x - xc); //XXXX change later - Nazariy
        //         T xscale = .5 * rn * it;  // scaling factor for x (how wide) replace with a slider later.
        //         T yscale = .5 * radius * sin(.0314*it);  // scaling factor for y (how close to center) dependent on iteration
        //         // the representation is y = e ^ (-x^2), but the example is set up in Z direction so here 
        //         // the x stands for domain in Z i.e. values from Z axis and y stands for value being subtracted from the radius
        //         int cntZ = (iZ+relativeOffset.z-(zl / 2))/xscale; // position of Z with resepect to center of Z        
        //         T modRad = yscale * exp(-cntZ*cntZ) + 1; // radius modification parameter
        //         if (dist > (radius - modRad)*(radius - modRad))
        //         {
        //             lattice.attributeDynamics(iX, iY, iZ, new BounceBack<T, DESCRIPTOR>);
        //         }
        //         else
        //         {
        //             lattice.attributeDynamics(iX, iY, iZ, new GuoExternalForceBGKdynamics<T, DESCRIPTOR>(1./0.95));
        //         }
        //         // last Z layer 
        //         iZ = domain.z1;
        //         dist = (iY + relativeOffset.y - yc)*(iY + relativeOffset.y - yc) + 
        //                     (iX + relativeOffset.x - xc)*(iX + relativeOffset.x - xc); //XXXX change later - Nazariy
        //         xscale = .5 * rn * it;  // scaling factor for x (how wide) replace with a slider later.
        //         yscale = .5 * radius * sin(.0314*it);  // scaling factor for y (how close to center) dependent on iteration
        //         // the representation is y = e ^ (-x^2), but the example is set up in Z direction so here 
        //         // the x stands for domain in Z i.e. values from Z axis and y stands for value being subtracted from the radius
        //         cntZ = (iZ+relativeOffset.z-(zl / 2))/xscale; // position of Z with resepect to center of Z        
        //         modRad = yscale * exp(-cntZ*cntZ) + 1; // radius modification parameter
        //         if (dist > (radius - modRad)*(radius - modRad))
        //         {
        //             lattice.attributeDynamics(iX, iY, iZ, new BounceBack<T, DESCRIPTOR>);
        //         }
        //         else
        //         {
        //             lattice.attributeDynamics(iX, iY, iZ, new GuoExternalForceBGKdynamics<T, DESCRIPTOR>(1./0.95));
        //         }
        //     }
        // }
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::dataStructure;
    }
    virtual DynamicBoundaryFunctional<T, Descriptor> *clone() const
    {
        return new DynamicBoundaryFunctional<T, Descriptor>(*this);
    }

private:
    plint xc,yc,it,zl; // zl added by Nazariy 7/19
    plint radius;
    T rn; //XXXX added by Nazariy 7/12
};

/// Automatic instantiation of the bounce-back nodes for the boundary,
///   using a data processor.
void createDynamicBoundaryFromDataProcessor(
    MultiBlockLattice3D<T, DESCRIPTOR> &lattice, plint xc, plint yc, plint radius, T rn, plint it, plint zl)
{
    applyProcessingFunctional(
        new DynamicBoundaryFunctional<T, DESCRIPTOR>(xc, yc, radius, rn, it, zl), lattice.getBoundingBox(),
        lattice);
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
    vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", dx/dt);
    vtkOut.writeData<3,float>(*computeVelocity(lattice), "velocity", dx/dt);
    vtkOut.writeData<3,float>(*computeVorticity(*computeVelocity(lattice)), "vorticity", 1./dt);
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
    T timeduration = T();
    global::timer("mainloop").start();

    //const plint N = atoi(argv[1]);
    const plint N = 1;// atoi(argv[1]);
    const T Re = 5e-3;
    const plint Nref = 50;
    //const T uMaxRef = 0.01;
    const T uMax = 0.00075;//uMaxRef /(T)N * (T)Nref; // Needed to avoid compressibility errors.

    IncomprFlowParam<T> parameters(
            uMax,
            Re,
            N,
            50,        // lx
            50,        // ly
            120         // lz
    );

    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    const plint nz = parameters.getNz();
    const T maxT = 2000;//6.6e4; //(T)0.01;
    plint iSave  = 20;//2000;//10;
    plint iCheck = 1000*iSave;
    // plint periT = 40000;
    writeLogFile(parameters, "3D square Poiseuille");
    
    // inlet and outlet boxes
    Box3D inlet    = Box3D(1,    nx-2, 1,    ny-2, 0,   0);
    Box3D outlet   = Box3D(1,    nx-2, 1,    ny-2, nz-1,   nz-1);

    // radii for the inlet and outlet holes
    plint r_inlet = 8; // edit this radius!!!!!!!!!!!!!
    plint r_smaller_outlet = 3;
    plint r_larger_outlet = 5;
    // vector of multiple radii
    std::vector<plint> outlet_rs = {r_smaller_outlet, r_larger_outlet};
    
    // centers for the inlet and outlet holes
    Array<T,3> center_inlet(21, 31, 2);
    Array<T,3> center_smaller_outlet(42.5, 17.5, 2);
    Array<T,3> center_larger_outlet(11, 13.5, 2);
    
    // vector of multiple outlet centers
    std::vector<Array<T,3>> outlet_centers = {center_smaller_outlet, center_larger_outlet};


    LammpsWrapper wrapper(argv,global::mpi().getGlobalCommunicator());
    char * inlmp = argv[1];
    wrapper.execFile(inlmp);
   
    //MultiTensorField3D<T,3> vel(parameters.getNx(),parameters.getNy(),parameters.getNz());
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
    
    //Cell<T,DESCRIPTOR> &cell = lattice.get(550,5500,550);
    pcout<<"dx "<<parameters.getDeltaX()<<" dt  "<<parameters.getDeltaT()<<" tau "<<parameters.getTau()<<endl;
    //pcout<<"51 works"<<endl;

/*
    MultiBlockLattice3D<T, DESCRIPTOR> lattice (
        parameters.getNx(), parameters.getNy(), parameters.getNz(), 
        new DYNAMICS );*/

    // Use periodic boundary conditions.
    // lattice.periodicity().toggle(2,true); // 2/9/2024 Nazariy: turning off periodicity in Z direction

    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
        = createLocalBoundaryCondition3D<T,DESCRIPTOR>();
    // OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
    //     = createInterpBoundaryCondition3D<T,DESCRIPTOR>();

    squarePoiseuilleSetup(lattice, parameters, *boundaryCondition);

    // Loop over main time iteration.
    util::ValueTracer<T> converge(parameters.getLatticeU(),parameters.getResolution(),1.0e-3);
    // uncomment to use force BC instead:  
        // Array<T,3> force(0,0,1e-6);
        // setExternalVector(lattice,lattice.getBoundingBox(),DESCRIPTOR<T>::ExternalField::forceBeginsAt,force);     
   
    plint xc,yc,radius, iterationCAS, zLength;

    xc = 20; // X center 
    yc = 20; // Y center
    radius = 15; // radius // not this radius, it's one for bounceback nodes edit the other one!!!!!!!!!!!!!!
    T radiusNorm = radius/maxT; // double radius/max iterations
    iterationCAS = 200; // iterations for collideAndStream in the loop 
    zLength = parameters.getNz(); // because domain.z0 gives local value pass this instead.  

    T targetV = 0.00075;
    // Array<T,3> targetV(0,0,uMax);
    // inlet target velocity 
    Array<T,3> inlet_targetV(0,0,targetV);

    // outlet target velocity
    // T outlet_scaled_targetV = targetV * (r_inlet*r_inlet)/((r_larger_outlet*r_larger_outlet)+(r_smaller_outlet*r_smaller_outlet));
    // T larger_outlet_scaled_targetV = outlet_scaled_targetV * (r_smaller_outlet*r_smaller_outlet)/((r_larger_outlet*r_larger_outlet)+(r_smaller_outlet*r_smaller_outlet));
    // T smaller_outlet_scaled_targetV = outlet_scaled_targetV * (r_larger_outlet*r_larger_outlet)/((r_larger_outlet*r_larger_outlet)+(r_smaller_outlet*r_smaller_outlet));
    
    // Array<T,3> small_outlet_targetV(0,0,smaller_outlet_scaled_targetV);
    // Array<T,3> large_outlet_targetV(0,0,larger_outlet_scaled_targetV);

    Array<T,3> small_outlet_targetV(0,0,targetV); // temp for commit purpose
    Array<T,3> large_outlet_targetV(0,0,targetV); // temp for commit purpose
    vector<Array<T,3>> outlet_targetV = {small_outlet_targetV, large_outlet_targetV};

    // setting the boundary conditions for the inlet and outlet holes:
    setBoundaryVelocity(lattice, inlet, SquarePoiseuilleVelocityHole<T>(center_inlet, r_inlet, inlet_targetV));  
    setBoundaryVelocity(lattice, outlet, MultipleSquarePoiseuilleVelocityHoles<T>(outlet_centers, outlet_rs, outlet_targetV));


    // createDynamicBoundaryFromDataProcessor(lattice, xc, yc, radius, radiusNorm, 0, zLength); // added by NT 7/18/2022
    plint prelim = 3e3;
    for (plint iT=0;iT<prelim;iT++){
        // progressbar
        plint progress = (plint) (iT * 100.0 / prelim);
        if(progress % 1 == 0) {
            int barWidth = 50;
            std::cout << "[";
            for (int j = 0; j < barWidth; ++j) {
                if (j < (barWidth*iT/prelim)) {
                    std::cout << "=";
                } else {
                    std::cout << " ";
                }
            }
            std::cout << "] " << int(progress) << "%" << "\r";
            std::cout.flush();
        }
        lattice.collideAndStream();
    }
 //           writeVTK(lattice, parameters, 4e3);
    // T timeduration = T();
    // global::timer("mainloop").start();
    for (plint iT=0; iT<=maxT; ++iT) {
    //for (plint iT=0; iT<2; ++iT) {
        if (iT%iSave ==0 && iT >0){
            pcout<<"Saving VTK file..."<<endl;
            writeVTK(lattice, parameters, iT);
        }
        // setBoundaryVelocity(lattice, inlet, SquarePoiseuilleVelocityHole<T>(center,r,targetV));  
        // setBoundaryVelocity(lattice, outlet, SquarePoiseuilleVelocityHole<T>(center,r,targetV));
        // lammps to calculate force
        wrapper.execCommand("run 1 pre no post no");
        // Clear and spread fluid force
        // createDynamicBoundaryFromDataProcessor(lattice, xc, yc, radius, radiusNorm, iT, zLength);
        ////-----classical ibm coupling-------------//
        spreadForce3D(lattice,wrapper);
        ////// Lattice Boltzmann iteration step.
        // for(int iteration=0; iteration<iterationCAS; iteration++){
        lattice.collideAndStream();
        // }
        // lattice.collideAndStream();
        ////// Interpolate and update solid position
        interpolateVelocity3D(lattice,wrapper);
        //-----force FSI ibm coupling-------------//
        //forceCoupling3D(lattice,wrapper);
        //lattice.collideAndStream();
    }

    timeduration = global::timer("mainloop").stop();
    pcout<<"total execution time "<<timeduration<<endl;
    delete boundaryCondition;
}
