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

//#include "opts.h"

#include "mpi.h"
#include "lammps.h"
#include "input.h"
#include "library.h"
#include "lammpsWrapper.h"
#include "update.h"
#include "neighbor.h"

#include "latticeDecomposition.h"
//#include "nearestTwoNeighborLattices3D.h"
// #include "senseiConfig.h"
//#ifdef ENABLE_SENSEI
#include <svtkVersion.h>
#include <svtkImageData.h>
#include <vtkXMLImageDataWriter.h>
// #include <vtkUnsignedCharArray.h>
#include <svtkUnsignedCharArray.h>
#include <svtkPointData.h>
#include <svtkDoubleArray.h>
#include <svtkSmartPointer.h>
#include <svtkUniformGrid.h> 
#include "Bridge.h"
//#endif

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
        u[2] = poiseuilleVelocity(iX, iY, parameters, maxN);
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



void squarePoiseuilleSetup( MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
                            IncomprFlowParam<T> const& parameters,
                            OnLatticeBoundaryCondition3D<T,DESCRIPTOR>& boundaryCondition )
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    const plint nz = parameters.getNz();
    Box3D top    = Box3D(0,    nx-1, ny-1, ny-1, 0, nz-1);
    Box3D bottom = Box3D(0,    nx-1, 0,    0,    0, nz-1);
    
    //Box3D inlet  = Box3D(0,    nx-1, 1,    ny-2, 0,    0);
    //Box3D outlet = Box3D(0,    nx-1, 1,    ny-2, nz-1, nz-1);
    
    Box3D left   = Box3D(0,    0,    1,    ny-2, 1, nz-2);
    Box3D right  = Box3D(nx-1, nx-1, 1,    ny-2, 1, nz-2);
    // shear flow top bottom surface
    /*
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, inlet, boundary::outflow );
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, outlet, boundary::outflow );

    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, top );
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, bottom );
    
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, left, boundary::outflow );
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, right, boundary::outflow );
    
    setBoundaryVelocity(lattice, top, ShearTopVelocity<T>(parameters,NMAX));
    setBoundaryVelocity(lattice, bottom, ShearBottomVelocity<T>(parameters,NMAX));
    
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, inlet, boundary::outflow );
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, outlet, boundary::outflow );
    */
    // channel flow
    //boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, inlet);
    //boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, outlet);
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, top );
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, bottom );
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, left );
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, right );
    
    //setBoundaryVelocity(lattice, inlet, SquarePoiseuilleVelocity<T>(parameters, NMAX));
    //setBoundaryVelocity(lattice, outlet, SquarePoiseuilleVelocity<T>(parameters, NMAX));
    
    setBoundaryVelocity(lattice, top, Array<T,3>((T)0.0,(T)0.0,(T)0.0));
    setBoundaryVelocity(lattice, bottom, Array<T,3>((T)0.0,(T)0.0,(T)0.0));
    setBoundaryVelocity(lattice, left, Array<T,3>((T)0.0,(T)0.0,(T)0.0));
    setBoundaryVelocity(lattice, right, Array<T,3>((T)0.0,(T)0.0,(T)0.0));
    

    //initializeAtEquilibrium(lattice, lattice.getBoundingBox(), SquarePoiseuilleDensityAndVelocity<T>(parameters, NMAX));
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(),(T)1.0, Array<T,3>(0.0,0.0,0.0));

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
    vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", dx/dt);
    vtkOut.writeData<3,float>(*computeVelocity(lattice), "velocity", dx/dt);
    vtkOut.writeData<3,float>(*computeVorticity(*computeVelocity(lattice)), "vorticity", 1./dt);
}
//**************************************Connor Changed 2/17/22
void writeVTK(MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
              Box3D domain, plint iter)
{
    //T dx = parameters.getDeltaX();
    //T dt = parameters.getDeltaT();
    VtkImageOutput3D<T> vtkOut(createFileName("vtk", iter, 6), 1);
    //vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", dx/dt);
    vtkOut.writeData<3,float>(*computeVelocity(lattice,domain), "velocity", 1);
    //vtkOut.writeData<3,float>(*computeVorticity(*computeVelocity(lattice)), "vorticity", 1./dt);
}
//**************************************

void VtkPalabos(MultiBlockLattice3D<T, DESCRIPTOR>& lattice,
    IncomprFlowParam<T> const& parameters,  plint iter)
{

	MultiTensorField3D<double, 3> velocityArray=*computeVelocity(lattice);
 	MultiTensorField3D<double, 3> VorticityArray=*computeVorticity(*computeVelocity(lattice));
	MultiScalarField3D<double> VelocityNormArray=*computeVelocityNorm(lattice);

 	pcout << VelocityNormArray << endl; 
	int  nx = parameters.getNx();  
	int  ny = parameters.getNy();  
	int  nz = parameters.getNz();
  
	svtkSmartPointer<svtkImageData> imageData =
        	svtkSmartPointer<svtkImageData>::New();

	imageData->SetDimensions(nx, ny, nz);

	//vtkSmartPointer<vtkDoubleArray> VelocityValues =
        //	vtkSmartPointer<vtkDoubleArray>::New();
	
	svtkDoubleArray *VelocityValues = svtkDoubleArray::New(); 	 
	
	VelocityValues->SetNumberOfComponents(3);
	VelocityValues->SetNumberOfTuples(nx * ny * nz); 

	//vtkSmartPointer<vtkDoubleArray> VorticityValues =
          //      vtkSmartPointer<vtkDoubleArray>::New();
      
	svtkDoubleArray *VorticityValues = svtkDoubleArray::New();
        
	VorticityValues->SetNumberOfComponents(3);
        VorticityValues->SetNumberOfTuples(nx * ny * nz);

//	vtkSmartPointer<vtkDoubleArray> VelocityNormValues =
  //              vtkSmartPointer<vtkDoubleArray>::New();

	svtkDoubleArray *VelocityNormValues = svtkDoubleArray::New();
	
        VelocityNormValues->SetNumberOfComponents(1);
        VelocityNormValues->SetNumberOfTuples(nx * ny * nz);
 
	for (int i=0; i<nz; i++)  
	{
 		for (int j=0; j<ny; j++)
        	{
                	for (int k=0; k<nx; k++) 
                	{	
			Array<double,3> vel = velocityArray.get(k,j,i); 
			Array<double,3> vor = VorticityArray.get(k,j,i); 
			double norm = VelocityNormArray.get(k,j,i);
  
			int index = j * nx + k + i * nx * ny; 
   
			VelocityValues->SetTuple3(index,vel[0],vel[1],vel[2]);
			VorticityValues->SetTuple3(index,vor[0],vor[1],vor[2]);
			VelocityNormValues->SetTuple1(index,norm);	  
			}
		}
	}


	imageData->GetPointData()->AddArray(VelocityValues);
		VelocityValues->SetName("Velocity");

	imageData->GetPointData()->AddArray(VorticityValues);
        VorticityValues->SetName("Vorticity");
	
    imageData->GetPointData()->AddArray(VelocityNormValues); // ass these lines to add Array pb_vel
        VelocityNormValues->SetName("Velocity Norm");

//	vtkSmartPointer<vtkXMLImageDataWriter> writer =
//		vtkSmartPointer<vtkXMLImageDataWriter>::New();

//	char filename[64];
//      sprintf (filename, "VtkDataStruc%d.vti", iter);

//	writer->SetInputData(imageData);
//	writer->SetFileName(filename);
//	writer->Write();
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
    const T Re = 5e-3;
    //const plint Nref = 50;
    //const T uMaxRef = 0.01;
    const T uMax = 0.00075;//uMaxRef /(T)N * (T)Nref; // Needed to avoid compressibility errors
    const int nx = 30;
    const int ny = 30;
    const int nz = 80;
    //using namespace opts;

    pcout<<"before configFile"<< endl;
    std::string config_file("/home/tishn/myFork/singleCell4.0/cellFlow.xml");//Configuration file to tell SENSEI what to do with data.
    // char * config_file = argv[4];
    // std::string  config_file = argv[4];
    pcout<<"before BRIDGE INITIALIZE"<< endl;
    Bridge::Initialize(global::mpi().getGlobalCommunicator(), config_file); // replaced with MPI_COMM_WORLD 
    // Bridge::Initialize(MPI_COMM_WORLD, config_file);
    /*Options ops(argc, argv);
    ops
    #ifdef ENABLE_SENSEI
    >> Option('f', "config", config_file, "Sensei analysis configuration xml (required)")
    #endif
    */
    IncomprFlowParam<T> parameters(
            uMax,
            Re,
            N,
            nx,        // lx
            ny,        // ly
            nz         // lz
    );
    const T maxT = atoi(argv[2]);//6.6e4; //(T)0.01;
    plint iSave = atoi(argv[3]);//10;//2000;//10;
    //plint iCheck = 10*iSave;
    writeLogFile(parameters, "3D square Poiseuille");

    LammpsWrapper wrapper(argv,global::mpi().getGlobalCommunicator()); // replaced with MPI_COMM_WORLD
    // LammpsWrapper wrapper(argv,MPI_COMM_WORLD);

    char * inlmp = argv[1];
    wrapper.execFile(inlmp);
   
    //MultiTensorField3D<T,3> vel(parameters.getNx(),parameters.getNy(),parameters.getNz());
    plint mysize = global::mpi().getSize();
    plint localdomain[mysize][6]; //First index: Rank value Second Index: extents (0: xlo 1: xhi 2: ylo 3: yhi 4: zlo 5: zhi)

    pcout<<"Nx,Ny,Nz "<<parameters.getNx()<<" "<<parameters.getNy()<<" "<<parameters.getNz()<<endl;
    LatticeDecomposition lDec(parameters.getNx(),parameters.getNy(),parameters.getNz(),
                              wrapper.lmp, localdomain); //XXX sending a double array to store extents of each processor for palabos data (Connor Murphy 2/19/22)

    SparseBlockStructure3D blockStructure = lDec.getBlockDistribution();
    ExplicitThreadAttribution* threadAttribution = lDec.getThreadAttribution();
    plint envelopeWidth = 3;


    MultiBlockManagement3D management = MultiBlockManagement3D(blockStructure, threadAttribution, envelopeWidth);//XXX New Change
    
    

    MultiBlockLattice3D<T, DESCRIPTOR> 
      lattice (management,
               defaultMultiBlockPolicy3D().getBlockCommunicator(),
               defaultMultiBlockPolicy3D().getCombinedStatistics(),
               defaultMultiBlockPolicy3D().getMultiCellAccess<T,DESCRIPTOR>(),
               new DYNAMICS );
    
    
    //********************************************* XXX Possible alternative for passing local extents to computeVelocity function
    /*
    ThreadAttribution const & orgThreadAttribution = management.getThreadAttribution();
    std::vector<plint> localBlocks = blockStructure.getLocalBlocks(orgThreadAttribution);
    std::map<plint, Box3D> bulksMap = blockStructure.getBulks();
    std::vector<Box3D> bulks;
    plint blockId;
    Box3D bulk;
    auto it=bulksMap.begin();
    cout << " bulkSmap size " << bulksMap.size() << endl;
    for(;it!=bulksMap.end();++it)
    {
        bulk = it->second;
        blockId = it->first;
        cout<<"block id= " << blockId<< " proc " << myrank <<" bulk: Nx " << bulk.getNx() << " Ny " << bulk.getNy() << " Nz " << bulk.getNz() << endl;
    }
    */
    //**********************************************
    
    //Cell<T,DESCRIPTOR> &cell = lattice.get(550,5500,550);
    pcout<<"dx "<<parameters.getDeltaX()<<" dt  "<<parameters.getDeltaT()<<" tau "<<parameters.getTau()<<endl;
    //pcout<<"51 works"<<endl;

    /*
    MultiBlockLattice3D<T, DESCRIPTOR> lattice (
        parameters.getNx(), parameters.getNy(), parameters.getNz(), 
        new DYNAMICS );*/

    // Use periodic boundary conditions.
    lattice.periodicity().toggle(2,true);
   

    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
        = createLocalBoundaryCondition3D<T,DESCRIPTOR>();

    squarePoiseuilleSetup(lattice, parameters, *boundaryCondition);

    // Loop over main time iteration.
    util::ValueTracer<T> converge(parameters.getLatticeU(),parameters.getResolution(),1.0e-3);
    //coupling between lammps and palabos
    Array<T,3> force(0,0,1e-6);
    setExternalVector(lattice,lattice.getBoundingBox(),DESCRIPTOR<T>::ExternalField::forceBeginsAt,force);
    //LAMMPS

    long time = 0; 
 
    for (plint iT=0;iT<4e3;iT++){
        lattice.collideAndStream();
    }
    T timeduration = T();
    global::timer("mainloop").start();
   plint nlocal;
   long ntimestep;
   int nanglelist;
   int nghost;
   double **x;
   
   int **anglelist;

   plint myrank = global::mpi().getRank();
   MultiTensorField3D<T,3> vel(lattice);
   MultiTensorField3D<double, 3> vort(lattice);
   MultiScalarField3D<double> velNorm(lattice);
   //TensorField3D<T,3> velocityArray = vel.getComponent(myrank);
   //TensorField3D<T,3> vorticityArray = vort.getComponent(myrank);
   //ScalarField3D<T> velocityNormArray = velNorm.getComponent(myrank);
   for (plint iT=0; iT<maxT; ++iT) {
        
        // lammps to calculate force
        wrapper.execCommand("run 1 pre no post no");

        //Some values are dynamically changing
        nlocal = wrapper.lmp->atom->nlocal;
        ntimestep = wrapper.lmp->update->ntimestep;
        nanglelist = wrapper.lmp->neighbor->nanglelist;
        nghost = wrapper.lmp->atom->nghost;
        x = wrapper.lmp->atom->x;
        anglelist = wrapper.lmp->neighbor->anglelist;
        

        //*************************************
        vel = *computeVelocity(lattice,lattice.getBoundingBox());
        vort = *computeVorticity(vel);
        velNorm = *computeVelocityNorm(lattice,lattice.getBoundingBox());

	TensorField3D<T,3> velocityArray = vel.getComponent(myrank);
   	TensorField3D<T,3> vorticityArray = vort.getComponent(myrank);
      	ScalarField3D<T> velocityNormArray = velNorm.getComponent(myrank);
      
        Box3D domain = Box3D(localdomain[myrank][0]-envelopeWidth,localdomain[myrank][1]+envelopeWidth,localdomain[myrank][2]-envelopeWidth,localdomain[myrank][3]+envelopeWidth,localdomain[myrank][4]-envelopeWidth,localdomain[myrank][5]+envelopeWidth);
        //*************************************
        
        //cout<<"Rank: " << myrank <<" Vorticity Extents: " <<vorticityArray.getNx() << " " << vorticityArray.getNy() << " " << vorticityArray.getNz()<<endl;
        //cout<<"Rank: " << myrank <<" Velocity Extents: " <<velocityArray.getNx() << " " << velocityArray.getNy() << " " << velocityArray.getNz()<<endl;
        //cout<<"Rank: " << myrank <<" Velocity Norm Extents: " <<velocityNormArray.getNx() << " " << velocityNormArray.getNy() << " " << velocityNormArray.getNz()<<endl;
        if (iT%iSave ==0 && iT >0){
        Bridge::SetData(x, ntimestep, nghost ,nlocal, anglelist, nanglelist,
			            velocityArray, vorticityArray, velocityNormArray, 
                        nx, ny, nz, domain, envelopeWidth);
        sensei::DataAdaptor *daOut = nullptr;
        Bridge::Analyze(time++, &daOut);
        }

        // Clear and spread fluid force
        setExternalVector(lattice,lattice.getBoundingBox(),DESCRIPTOR<T>::ExternalField::forceBeginsAt,force);
        ////-----classical ibm coupling-------------//
        spreadForce3D(lattice,wrapper);
        ////// Lattice Boltzmann iteration step.
        lattice.collideAndStream();
        ////// Interpolate and update solid position
        interpolateVelocity3D(lattice,wrapper);
        //-----force FSI ibm coupling-------------//
        //forceCoupling3D(lattice,wrapper);
        //lattice.collideAndStream();
        //writeVTK(lattice, domainBox, iT);
	
    }

    timeduration = global::timer("mainloop").stop();
    pcout<<"total execution time "<<timeduration<<endl;
    delete boundaryCondition;
    Bridge::Finalize();
}
