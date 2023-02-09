#pragma once
#include "palabos3D.h"
#include "palabos3D.hh"
#include <DataAdaptor.h>
#include <svtkDoubleArray.h>
#include <svtkIntArray.h>
#include <svtkImageData.h>
// #include <vtkUnsignedCharArray.h>
#include <svtkUnsignedCharArray.h>
using namespace plb; 
namespace senseiLP 
{
class LPDataAdaptor : public sensei::DataAdaptor
{
public:

  static LPDataAdaptor* New();

  senseiTypeMacro(LPDataAdaptor, sensei::DataAdaptor);

  void Initialize();

  void AddLAMMPSData(double **x, long ntimestep, int nghost, 
                     int nlocal, int **anglelist, int nanglelist);

  void AddPalabosData(svtkDoubleArray *velocityDoubleArray,
                      svtkDoubleArray *vorticityDoubleArray,
                      svtkDoubleArray *velocityNormDoubleArray,
                      int nx, int ny, int nz, Box3D domainBox, plint envelopeWidth, svtkDoubleArray *center);
// SENSEI API (Virtual functions overridden from sensei/DataAdaptor.h)
  int GetNumberOfMeshes(unsigned int &numMeshes) override;

  int GetMeshMetadata(unsigned int id, sensei::MeshMetadataPtr &metadata) override;

  int GetMesh(const std::string &meshName, bool structureOnly, svtkDataObject *&mesh); // override ;

  int GetMesh(const std::string &meshName, bool structureOnly, svtkCompositeDataSet *&mesh); // override ;

  int AddGhostNodesArray(svtkDataObject* mesh, const std::string &meshName) override;

  int AddGhostCellsArray(svtkDataObject* mesh, const std::string &meshName) override;

  int AddArray(svtkDataObject* mesh, const std::string &meshName, int association, const std::string &arrayName) override;

  int AddArrays(svtkDataObject* mesh, const std::string &meshName, int association, const std::vector<std::string> &arrayName) override;

  int ReleaseData() override;

protected:

LPDataAdaptor();

~LPDataAdaptor();

private:

  struct DInternals;
  DInternals* Internals;
};

}
