#include "Bridge.h"
#include "LPdataAdaptor.h"
#include <svtkSmartPointer.h>
#include <ConfigurableAnalysis.h>
#include <iostream>

using namespace std;
using namespace plb; 
namespace Bridge
{
   static svtkSmartPointer<senseiLP::LPDataAdaptor>  GlobalDataAdaptor;
   static svtkSmartPointer<sensei::ConfigurableAnalysis> GlobalAnalysisAdaptor;
  
   svtkDoubleArray *velocityDoubleArray = svtkDoubleArray::New(); //jifu: 5/31/22, define it as global and release it later
   svtkDoubleArray *vorticityDoubleArray = svtkDoubleArray::New();
   svtkDoubleArray *velocityNormDoubleArray = svtkDoubleArray::New();

   // bidirectional_proxies arrays:
   svtkDoubleArray *center = svtkDoubleArray::New();

void Initialize(MPI_Comm world, const std::string& config_file){
   
   GlobalDataAdaptor = svtkSmartPointer<senseiLP::LPDataAdaptor>::New();
   GlobalDataAdaptor->Initialize();
   GlobalDataAdaptor->SetCommunicator(world);
   GlobalDataAdaptor->SetDataTimeStep(-1); //XXX Why -1?

   GlobalAnalysisAdaptor = svtkSmartPointer<sensei::ConfigurableAnalysis>::New();

   GlobalAnalysisAdaptor->Initialize(config_file);

}
void SetData(double **x, long ntimestep, int nghost, 
             int nlocal, int **anglelist, int nanglelist, 
	         TensorField3D<double, 3> velocityArray,
	         TensorField3D<double, 3> vorticityArray,
	         ScalarField3D<double> velocityNormArray,
           int nx, int ny, int nz, Box3D domainBox, plint envelopeWidth)
{
  GlobalDataAdaptor->AddLAMMPSData(x, ntimestep, nghost, nlocal, anglelist, nanglelist);

  /*vtkDoubleArray *velocityDoubleArray = vtkDoubleArray::New();
  vtkDoubleArray *vorticityDoubleArray = vtkDoubleArray::New();
  vtkDoubleArray *velocityNormDoubleArray = vtkDoubleArray::New();
  */

//XXXNew local values added with domainBox 2/23/22*****
  int nlx = velocityArray.getNx(); 
  int nly = velocityArray.getNy();
  int nlz = velocityArray.getNz();
  plint myrank = global::mpi().getRank();

//*****************************************************
  velocityDoubleArray->SetNumberOfComponents(3);
  velocityDoubleArray->SetNumberOfTuples((nlx) * (nly) * (nlz)); 

  vorticityDoubleArray->SetNumberOfComponents(3);
  vorticityDoubleArray->SetNumberOfTuples((nlx) * (nly) * (nlz));

  velocityNormDoubleArray->SetNumberOfComponents(1);
  velocityNormDoubleArray->SetNumberOfTuples((nlx) * (nly) * (nlz));

  center->SetNumberOfComponents(3);
  center->SetNumberOfTuples((nlx) * (nly) * (nlz));
  // center->SetNumberOfTuples(1);
  // center->SetTuple3(0,1.,2.,3.); // x, y, z

  // center->SetNumberOfComponents(3);
  // center->SetNumberOfTuples((nlx) * (nly) * (nlz));
  // center->SetTuple3(0,0.,0.,0.);

  // coordsArray->SetNumberOfComponents(1);
  // coordsArray->SetNumberOfTuples(3);
  // Array<int,3> center(0,0,0);
   //plint EW = envelopeWidth;

//XXX Need to convert this to zero copy: FUTURE WORK

  Array<double,3> vel(0.,0.,0.); //jifu 5/31/2022
  Array<double,3> vor(0.,0.,0.);
        
  for (int k=0; k<nlz; k++)
  {
    for (int j=0; j<nly; j++)
    {
     for (int i=0; i<nlx; i++)
      {
        vel = velocityArray.get(i,j,k); //jifu 5/31/2022 Array<double 3> vel ->vel
        vor = vorticityArray.get(i,j,k);
        double norm = velocityNormArray.get(i,j,k);
        int index = (j) * (nlx) + (i) + (k) * (nlx) * (nly);
        velocityDoubleArray->SetTuple3(index,vel[0],vel[1],vel[2]);
        vorticityDoubleArray->SetTuple3(index,vor[0],vor[1],vor[2]);
        velocityNormDoubleArray->SetTuple1(index,norm);
        
        // all internals associated with each datapoint ???? 2/6/23
        center->SetTuple3(index,vor[0],vor[1],vor[2]); // THIS IS AN UGLY HACK WILL FIX LATER. 
      }
    }
  }
 GlobalDataAdaptor->AddPalabosData(velocityDoubleArray, vorticityDoubleArray, velocityNormDoubleArray, nx, ny, nz, domainBox, envelopeWidth, center);

}
void Analyze(long ntimestep, sensei::DataAdaptor **dataOut)
{
   GlobalDataAdaptor->SetDataTimeStep(ntimestep);
   GlobalDataAdaptor->SetDataTime(ntimestep);
   GlobalAnalysisAdaptor->Execute(GlobalDataAdaptor.GetPointer(), dataOut);
   GlobalDataAdaptor->ReleaseData();
}
void Finalize()
   {
   GlobalAnalysisAdaptor->Finalize();
   GlobalAnalysisAdaptor = NULL;
   GlobalDataAdaptor = NULL;
   velocityDoubleArray->Delete(); 
   vorticityDoubleArray->Delete(); 
   velocityNormDoubleArray->Delete(); // jifu: 5/31/2022 working, memory leakage is fixed 
   center->Delete(); 
   }
}

