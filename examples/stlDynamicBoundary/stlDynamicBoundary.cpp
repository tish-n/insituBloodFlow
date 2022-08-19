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

using namespace plb;
using namespace std;

typedef double T;
//#define DESCRIPTOR descriptors::ForcedN2D3Q19Descriptor
#define DESCRIPTOR descriptors::ForcedD3Q19Descriptor
//#define DYNAMICS BGKdynamics<T, DESCRIPTOR>(parameters.getOmega())
#define DYNAMICS GuoExternalForceBGKdynamics<T, DESCRIPTOR>(parameters.getOmega())
#define NMAX 150

const T pi = (T)4.*std::atan((T)1.);

plint extraLayer = 0;  // Make the bounding box larger; for visualization purposes
                            //   only. For the simulation, it is OK to have extraLayer=0.
const plint blockSize = 20; // Zero means: no sparse representation.
const plint envelopeWidth = 1;  // For standard BGK dynamics.
const plint extendedEnvelopeWidth = 2;  // Because the Guo off lattice boundary condition
                                        //   needs 2-cell neighbor access.

bool performOutput = false;
bool doImages = false;
bool useAllDirections = false;
bool useRegularizedWall = false;
bool useIncompressible = false;
bool poiseuilleInlet = false;
bool convectiveScaling = false;

T kinematicViscosity       = 0.;
T averageInletVelocity     = 0.;
plint referenceResolution  = 0.;
T nuLB                     = 0.;
T fluidDensity             = 0.;
T volume                   = 0.;
T userDefinedInletDiameter = 0.;

plint referenceDirection = 0;
plint openingSortDirection = 0;

T simTime = 0;
plint startLevel = 0;
plint maxLevel   = 0;
T epsilon = 0.;

TriangleSet<T>* triangleSet = 0;
T currentTime = 0;

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

// This function assigns proper boundary conditions to the openings of the surface geometry
//   of the aneurysm. Which opening is inlet and which is outlet is defined by the user in
//   the input XML file. For the inlet, there is a choice between a Poiseuille velocity
//   profile and a simple plug velocity profile. At the outlets a Neumann boundary
//   condition with constant pressure is prescribed.
void setOpenings (
    std::vector<BoundaryProfile3D<T,Velocity>*>& inletOutlets,
    TriangleBoundary3D<T>& boundary, T uLB, T dx, T dt )
{
    for (pluint i=0; i<openings.size(); ++i) {
        Opening<T>& opening = openings[i];
        opening.center = computeBaryCenter (
                boundary.getMesh(),
                boundary.getInletOutlet(openingSortDirection)[i] );
        opening.innerRadius = computeInnerRadius (
                boundary.getMesh(),
                boundary.getInletOutlet(openingSortDirection)[i] );

        if (opening.inlet) {
            if (poiseuilleInlet) {
                inletOutlets.push_back (
                        new PoiseuilleProfile3D<T>(uLB) );
            }
            else {
                inletOutlets.push_back (
                        new VelocityPlugProfile3D<T>(uLB) );
            }
        }
        else {
            inletOutlets.push_back (
                    new DensityNeumannBoundaryProfile3D<T> );
        }
    }
}

/// This functional defines a data processor for the instantiation
///   of bounce-back nodes following the half-circle geometry. 
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
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    T dist = (iY + relativeOffset.y - yc)*(iY + relativeOffset.y - yc) + 
                             (iX + relativeOffset.x - xc)*(iX + relativeOffset.x - xc); //XXXX change later - Nazariy
                    T xscale = .5 * rn * it;  // scaling factor for x (how wide) replace with a slider later.
                    T yscale = .5 * radius * sin(.0314*it);  // scaling factor for y (how close to center) dependent on iteration
                    // the representation is y = e ^ (-x^2), but the example is set up in Z direction so here 
                    // the x stands for domain in Z i.e. values from Z axis and y stands for value being subtracted from the radius
                    int cntZ = (iZ+relativeOffset.z-(zl / 2))/xscale; // position of Z with resepect to center of Z        
                    T modRad = yscale * exp(-cntZ*cntZ) + 1; // radius modification parameter
                    if ((dist > (radius - modRad)*(radius - modRad)) && !lattice.get(iX, iY, iZ).isBoundary()) {
                        lattice.attributeDynamics(iX, iY, iZ, new BounceBack<T, DESCRIPTOR>);
                    }
                    else{
                        lattice.attributeDynamics(iX, iY, iZ, new GuoExternalForceBGKdynamics<T, DESCRIPTOR>(1./0.95));
                    }
                }
            }
        }
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

// T computeRMSerror ( MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
//                     IncomprFlowParam<T> const& parameters )
// {
//     MultiTensorField3D<T,3> analyticalVelocity(lattice);
//     setToFunction( analyticalVelocity, analyticalVelocity.getBoundingBox(),
//                    SquarePoiseuilleVelocity<T>(parameters, NMAX) );
//     MultiTensorField3D<T,3> numericalVelocity(lattice);
//     computeVelocity(lattice, numericalVelocity, lattice.getBoundingBox());

//            // Divide by lattice velocity to normalize the error
//     return 1./parameters.getLatticeU() *
//            // Compute RMS difference between analytical and numerical solution
//            std::sqrt( computeAverage( *computeNormSqr(
//                           *subtract(analyticalVelocity, numericalVelocity)
//                      ) ) );
// }

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


// This is the function that prepares and performs the actual simulation.
std::unique_ptr<MultiBlockLattice3D<T,DESCRIPTOR> > run (
        plint level, MultiBlockLattice3D<T,DESCRIPTOR>* iniVal=0 )
{
    plint margin = 3; // Extra margin of allocated cells around the obstacle. 
    plint borderWidth = 1; // Because the Guo boundary condition acts in a one-cell layer.
                           // Requirement: margin>=borderWidth.

    // The resolution is doubled at each coordinate direction with the increase of the
    //   resolution level by one. The parameter ``referenceResolution'' is by definition
    //   the resolution at grid refinement level 0.
    plint resolution = referenceResolution * util::twoToThePower(level);

    // The next few lines of code are typical. They transform the surface geometry of the
    //   aneurysm given by the user to more efficient data structures that are internally
    //   used by palabos. The TriangleBoundary3D structure will be later used to assign
    //   proper boundary conditions.
    DEFscaledMesh<T>* defMesh =
        new DEFscaledMesh<T>(*triangleSet, resolution, referenceDirection, margin, extraLayer);
    TriangleBoundary3D<T> boundary(*defMesh);
    delete defMesh;
    boundary.getMesh().inflate();

    // When convective scaling is used (relationship of dt with respect to dx as the grid is
    //   refined) the value of the kinematic viscosity must be also properly adjusted.
    T nuLB_ = nuLB;
    if (convectiveScaling) {
        nuLB_ = nuLB * util::twoToThePower(level);
    }
    T dx = boundary.getDx();
    T dt = nuLB_ / kinematicViscosity *dx*dx;
    T uAveLB = averageInletVelocity *dt/dx;
    T omega = 1./(3.*nuLB_+0.5);
    Array<T,3> location(boundary.getPhysicalLocation());


    pcout << "uLB=" << uAveLB << std::endl;
    pcout << "nuLB=" << nuLB_ << std::endl;
    pcout << "tau=" << 1./omega << std::endl;
    if (performOutput) {
        pcout << "dx=" << dx << std::endl;
        pcout << "dt=" << dt << std::endl;
    }

    // The aneurysm simulation is an interior (as opposed to exterior) flow problem. For
    //   this reason, the lattice nodes that lay inside the computational domain must
    //   be identified and distinguished from the ones that lay outside of it. This is
    //   handled by the following voxelization process.
    const int flowType = voxelFlag::inside;
    VoxelizedDomain3D<T> voxelizedDomain (
            boundary, flowType, extraLayer, borderWidth, extendedEnvelopeWidth, blockSize );
    if (performOutput) {
        pcout << getMultiBlockInfo(voxelizedDomain.getVoxelMatrix()) << std::endl;
    }

    MultiScalarField3D<int> flagMatrix((MultiBlock3D&)voxelizedDomain.getVoxelMatrix());
    setToConstant(flagMatrix, voxelizedDomain.getVoxelMatrix(),
                  voxelFlag::inside, flagMatrix.getBoundingBox(), 1);
    setToConstant(flagMatrix, voxelizedDomain.getVoxelMatrix(),
                  voxelFlag::innerBorder, flagMatrix.getBoundingBox(), 1);
    pcout << "Number of fluid cells: " << computeSum(flagMatrix) << std::endl;

    Dynamics<T,DESCRIPTOR>* dynamics = 0;
    if (useIncompressible) {
        dynamics = new IncBGKdynamics<T,DESCRIPTOR>(omega); // In this model velocity equals momentum.
    }
    else {
        dynamics = new BGKdynamics<T,DESCRIPTOR>(omega); // In this model velocity equals momentum
                                                         //   divided by density.
    }
    std::unique_ptr<MultiBlockLattice3D<T,DESCRIPTOR> > lattice 
        = generateMultiBlockLattice<T,DESCRIPTOR> (
                voxelizedDomain.getVoxelMatrix(), envelopeWidth, dynamics );
    lattice->toggleInternalStatistics(false);

    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
        = createLocalBoundaryCondition3D<T,DESCRIPTOR>();

    Array<T,3> inletRealPos(0.024, 0.018, 0.002);
    Array<T,3> outlet1RealPos(0.038, 0.06, 0.055);
    Array<T,3> outlet2RealPos(0.007,0.06,0.105);
    T diameterReal = 0.01;

    Array<T,3> inletPos((inletRealPos-location)/dx);
    Array<T,3> outlet1Pos((outlet1RealPos-location)/dx);
    Array<T,3> outlet2Pos((outlet2RealPos-location)/dx);
    plint diameter = util::roundToInt(diameterReal/dx);

    Box3D inletDomain(util::roundToInt(inletPos[0]-diameter), util::roundToInt(inletPos[0]+diameter),
                      util::roundToInt(inletPos[1]-diameter), util::roundToInt(inletPos[1]+diameter),
                      util::roundToInt(inletPos[2]), util::roundToInt(inletPos[2]));
    Box3D behindInlet(inletDomain.x0, inletDomain.x1,
                      inletDomain.y0, inletDomain.y1,
                      inletDomain.z0-diameter, inletDomain.z0-1);

    Box3D outlet1Domain(util::roundToInt(outlet1Pos[0]-diameter), util::roundToInt(outlet1Pos[0]+diameter),
                        util::roundToInt(outlet1Pos[1]-diameter), util::roundToInt(outlet1Pos[1]+diameter),
                        util::roundToInt(outlet1Pos[2]), util::roundToInt(outlet1Pos[2]));
    Box3D behindOutlet1(outlet1Domain.x0, outlet1Domain.x1,
                        outlet1Domain.y0, outlet1Domain.y1,
                        outlet1Domain.z0-diameter, outlet1Domain.z0-1);

    Box3D outlet2Domain(util::roundToInt(outlet2Pos[0]), util::roundToInt(outlet2Pos[0]),
                        util::roundToInt(outlet2Pos[1]-diameter), util::roundToInt(outlet2Pos[1]+diameter),
                        util::roundToInt(outlet2Pos[2]-diameter), util::roundToInt(outlet2Pos[2]+diameter));
    Box3D behindOutlet2(outlet2Domain.x0-diameter, outlet2Domain.x1-1,
                        outlet2Domain.y0, outlet2Domain.y1,
                        outlet2Domain.z0, outlet2Domain.z1);


    boundaryCondition->addVelocityBoundary2N(inletDomain, *lattice);
    setBoundaryVelocity(*lattice, inletDomain, Array<T,3>((T)0.,(T)0.,uAveLB));
    boundaryCondition->addPressureBoundary2N(outlet1Domain, *lattice);
    setBoundaryDensity(*lattice, outlet1Domain, (T)1.);
    boundaryCondition->addPressureBoundary0N(outlet2Domain, *lattice);
    setBoundaryDensity(*lattice, outlet2Domain, (T)1.);

    defineDynamics(*lattice, flagMatrix, lattice->getBoundingBox(), new BounceBack<T,DESCRIPTOR>(1.), 0);
    defineDynamics(*lattice, behindInlet, new BounceBack<T,DESCRIPTOR>(1.));
    defineDynamics(*lattice, behindOutlet1, new BounceBack<T,DESCRIPTOR>(1.));
    defineDynamics(*lattice, behindOutlet2, new BounceBack<T,DESCRIPTOR>(1.));

    iniLattice(*lattice, voxelizedDomain);
    if(iniVal) {
        Box3D toDomain(lattice->getBoundingBox());
        Box3D fromDomain(toDomain.shift(margin,margin,margin)); // During rescaling, the margin doubled in size,
                                                                //   an effect which is cancelled here through a shift.
        copy(*iniVal, fromDomain, *lattice, toDomain, modif::staticVariables);
    }

    // The ValueTracer is needed to check when a chosen quantity (in our case the average energy)
    //   has converged, so to conclude that steady state has been reached for the specific grid
    //   refinement level and stop the simulation.
    plint convergenceIter=20;
    util::ValueTracer<T> velocityTracer(0.05*convergenceIter, resolution, epsilon);
    global::timer("iteration").restart();
    plint i = util::roundToInt(currentTime/dt);
    lattice->resetTime(i);

    pcout << "Saving a " << lattice->getNx() << " by " << lattice->getNy()
          << " by " << lattice->getNz() << " lattice." << std::endl;
    global::timer("io").start();
    parallelIO::save(*lattice, "checkpoint", false);
    pcout << "Total time for i/o: " << global::timer("io").getTime() << std::endl;

    // Collision and streaming iterations.
    pcout << "Starting 100 iterations" << std::endl;
    global::timer("global").start();
    for (plint i=0; i<100; ++i) {

        lattice->collideAndStream();

        ++i;
        currentTime = i*dt;
    }
    pcout << "End of 100 iterations" << std::endl;
    pcout << "Total time of execution: " << global::timer("global").getTime() << std::endl;

    return lattice;
}

// Read the user input XML file provided at the command-line.
void readParameters(XMLreader const& document)
{
    std::string meshFileName;
    std::vector<std::string> openingType;
    document["geometry"]["mesh"].read(meshFileName);
    document["geometry"]["inletDiameter"].read(userDefinedInletDiameter);
    document["geometry"]["averageInletVelocity"].read(averageInletVelocity);
    document["geometry"]["openings"]["sortDirection"].read(openingSortDirection);
    document["geometry"]["openings"]["type"].read(openingType);

    document["fluid"]["kinematicViscosity"].read(kinematicViscosity);
    document["fluid"]["density"].read(fluidDensity);
    document["fluid"]["volume"].read(volume);

    document["numerics"]["referenceDirection"].read(referenceDirection);
    document["numerics"]["referenceResolution"].read(referenceResolution);
    document["numerics"]["nuLB"].read(nuLB);

    document["simulation"]["simTime"].read(simTime);
    document["simulation"]["maxLevel"].read(maxLevel);
    document["simulation"]["epsilon"].read(epsilon);

    document["simulation"]["performOutput"].read(performOutput);
    document["simulation"]["doImages"].read(doImages);
    document["simulation"]["useAllDirections"].read(useAllDirections);
    document["simulation"]["useRegularizedWall"].read(useRegularizedWall);
    document["simulation"]["useIncompressible"].read(useIncompressible);
    document["simulation"]["poiseuilleInlet"].read(poiseuilleInlet);
    document["simulation"]["convectiveScaling"].read(convectiveScaling);

    // At this part, the surface geometry of the aneurysm (as given by the user in
    //   the form of an ASCII or binary STL file) is read into a data structure
    //   comprised by a set of triangles. The DBL constant means that double
    //   precision accuracy will be used (generally the recommended choice).
    triangleSet = new TriangleSet<T>(meshFileName, DBL);
    pcout << "Reynolds number, based on provided inlet diameter: "
          << averageInletVelocity*userDefinedInletDiameter/kinematicViscosity
          << std::endl;
    plbIOError(openingSortDirection<0 || openingSortDirection>2,
               "Sort-direction of opening must be 0 (x), 1 (y), or 2 (z).");
    // The surface geometry, as provided by the STL file, must contain openings,
    //   namely inlets and outlets. On these openings, appropriate boundary conditions
    //   will be imposed by palabos. Which opening is inlet and which is outlet, is
    //   identified by the user in the input XML file.
    openings.resize(openingType.size());
    for (pluint i=0; i<openingType.size(); ++i) {
        std::string next_opening = util::tolower(openingType[i]);
        if (next_opening=="inlet") {
            openings[i].inlet = true;
        }
        else if (next_opening=="outlet") {
            openings[i].inlet = false;
        }
        else {
            plbIOError("Unknown opening type.");
        }
    }
}



int main(int argc, char* argv[]) {

    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");
    global::IOpolicy().activateParallelIO(true);
    
    string paramXmlFileName;
    try {
        global::argv(1).read(paramXmlFileName);
    }
    catch (PlbIOException& exception) {
        pcout << "Wrong parameters; the syntax is: " 
              << (std::string)global::argv(0) << " parameter-input-file.xml" << std::endl;
        return -1;
    }

    // Read the parameter XML input file. (Lots of comments are included there too).
    try {
        XMLreader document(paramXmlFileName);
        readParameters(paramXmlFileName);
    }
    catch (PlbIOException& exception) {
        pcout << "Error in input file " << paramXmlFileName
              << ": " << exception.what() << std::endl;
        return -1;
    }

    plint iniLevel=0;
    std::unique_ptr<MultiBlockLattice3D<T,DESCRIPTOR> > iniConditionLattice(nullptr);
    // This code incorporates the concept of smooth grid refinement until convergence is
    //   achieved. The word ``smooth'' indicates that as the refinement level increases
    //   by one, the whole grid doubles in each direction. When the grid is refined, both
    //   dx and dt have to change. Whether dt is changed as dx^2 (diffusive behavior)
    //   or as dx (convective behavior), is controlled by the input variable
    //   ``convectiveScaling'' (the recommended choice is not to use convective scaling).
    try {
        for (plint level=iniLevel; level<=maxLevel; ++level) {
            pcout << std::endl << "Running new simulation at level " << level << std::endl;
            std::unique_ptr<MultiBlockLattice3D<T,DESCRIPTOR> > convergedLattice (
                    run(level, iniConditionLattice.get()) );
            if (level != maxLevel) {
                plint dxScale = -1;
                plint dtScale = -2;
                if (convectiveScaling) {
                    dtScale = -1;
                }
                // The converged simulation of the previous grid level is used as the initial condition
                //   for the simulation at the next grid level (after appropriate interpolation has
                //   taken place).
                iniConditionLattice = std::unique_ptr<MultiBlockLattice3D<T,DESCRIPTOR> > (
                        refine(*convergedLattice, dxScale, dtScale, new BGKdynamics<T,DESCRIPTOR>(1.)) );
            }
        }
    }
    catch(PlbException& exception) {
        pcout << exception.what() << std::endl;
        return -1;
    }
//     /*
//     if (argc != 2) {
//         pcout << "Error the parameters are wrong. The structure must be :\n";
//         pcout << "1 : N\n";
//         exit(1);
//     }*/
//     T timeduration = T();
//     global::timer("mainloop").start();

//     //const plint N = atoi(argv[1]);
//     const plint N = 1;// atoi(argv[1]);
//     const T Re = 5e-3;
//     const plint Nref = 50;
//     //const T uMaxRef = 0.01;
//     const T uMax = 0.00075;//uMaxRef /(T)N * (T)Nref; // Needed to avoid compressibility errors.

//     IncomprFlowParam<T> parameters(
//             uMax,
//             Re,
//             N,
//             40.,        // lx
//             40.,        // ly
//             80.         // lz
//     );
//     const T maxT = 100;//6.6e4; //(T)0.01;
//     plint iSave  = 1;//2000;//10;
//     plint iCheck = 1000*iSave;
//     writeLogFile(parameters, "3D square Poiseuille");

//     LammpsWrapper wrapper(argv,global::mpi().getGlobalCommunicator());
//     char * inlmp = argv[1];
//     wrapper.execFile(inlmp);
   
//     //MultiTensorField3D<T,3> vel(parameters.getNx(),parameters.getNy(),parameters.getNz());
//     pcout<<"Nx,Ny,Nz "<<parameters.getNx()<<" "<<parameters.getNy()<<" "<<parameters.getNz()<<endl;
//     LatticeDecomposition lDec(parameters.getNx(),parameters.getNy(),parameters.getNz(),
//                               wrapper.lmp);
//     SparseBlockStructure3D blockStructure = lDec.getBlockDistribution();
//     ExplicitThreadAttribution* threadAttribution = lDec.getThreadAttribution();
//     plint envelopeWidth = 3;

//     MultiBlockLattice3D<T, DESCRIPTOR> 
//       lattice (MultiBlockManagement3D (blockStructure, threadAttribution, envelopeWidth ),
//                defaultMultiBlockPolicy3D().getBlockCommunicator(),
//                defaultMultiBlockPolicy3D().getCombinedStatistics(),
//                defaultMultiBlockPolicy3D().getMultiCellAccess<T,DESCRIPTOR>(),
//                new DYNAMICS );
    
//     //Cell<T,DESCRIPTOR> &cell = lattice.get(550,5500,550);
//     pcout<<"dx "<<parameters.getDeltaX()<<" dt  "<<parameters.getDeltaT()<<" tau "<<parameters.getTau()<<endl;
//     //pcout<<"51 works"<<endl;

// /*
//     MultiBlockLattice3D<T, DESCRIPTOR> lattice (
//         parameters.getNx(), parameters.getNy(), parameters.getNz(), 
//         new DYNAMICS );*/

//     // Use periodic boundary conditions.
//     lattice.periodicity().toggle(2,true);

//     OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
//         = createLocalBoundaryCondition3D<T,DESCRIPTOR>();

//     squarePoiseuilleSetup(lattice, parameters, *boundaryCondition);

//     // Loop over main time iteration.
//     util::ValueTracer<T> converge(parameters.getLatticeU(),parameters.getResolution(),1.0e-3);
//       //coupling between lammps and palabos
//         Array<T,3> force(0,0.,1e-6);
//         setExternalVector(lattice,lattice.getBoundingBox(),DESCRIPTOR<T>::ExternalField::forceBeginsAt,force);     
   
//     plint xc,yc,radius, iterationCAS, zLength;

//     xc = 20; // X center 
//     yc = 20; // Y center
//     radius = 20; // radius
//     T radiusNorm = radius/maxT; // double radius/max iterations
//     iterationCAS = 10; // iterations for collideAndStream in the loop 
//     zLength = parameters.getNz(); // because domain.z0 gives local value pass this instead.  
    
//     createDynamicBoundaryFromDataProcessor(lattice, xc, yc, radius, radiusNorm, 0, zLength); // added by NT 7/18/2022
//     for (plint iT=0;iT<6e3;iT++){
//         lattice.collideAndStream();
//     }
//  //           writeVTK(lattice, parameters, 4e3);
//     // T timeduration = T();
//     // global::timer("mainloop").start();
//     for (plint iT=0; iT<=maxT; ++iT) {
//     //for (plint iT=0; iT<2; ++iT) {
//         if (iT%iSave ==0 && iT >0){
//             pcout<<"Saving VTK file..."<<endl;
//             writeVTK(lattice, parameters, iT);
//         }
//         if (iT%iCheck ==0 && iT >0){
//             pcout<<"Timestep "<<iT<<" Saving checkPoint file..."<<endl;
//             saveBinaryBlock(lattice,"checkpoint.dat");
//         }
//         // lammps to calculate force
//         //wrapper.execCommand("run 1 pre no post no");
//         // Clear and spread fluid force
//         createDynamicBoundaryFromDataProcessor(lattice, xc, yc, radius, radiusNorm, iT, zLength);
//         ////-----classical ibm coupling-------------//
//         //spreadForce3D(lattice,wrapper);
//         ////// Lattice Boltzmann iteration step.
//         for(int iteration=0; iteration<iterationCAS; iteration++){
//             lattice.collideAndStream();
//         }

//         for(int iteration=0; iteration<iterationCAS; iteration++){
//             lattice.bulkCollideAndStream();
//         }
//         // lattice.collideAndStream();
//         ////// Interpolate and update solid position
//         //interpolateVelocity3D(lattice,wrapper);
//         //-----force FSI ibm coupling-------------//
//         //forceCoupling3D(lattice,wrapper);
//         //lattice.collideAndStream();
//     }

//     timeduration = global::timer("mainloop").stop();
//     pcout<<"total execution time "<<timeduration<<endl;
//     delete boundaryCondition;
// }
