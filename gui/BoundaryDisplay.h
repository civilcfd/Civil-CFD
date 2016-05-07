/* 
 * BoundaryDisplay.h
 *
 * Descibe a class that shows a mesh and highlighed boundaries with a VTK file
 *
 * Inherits from RenderDisplay.h
 */

#include <vtkDataReader.h>
#include <vtkStructuredPointsReader.h>
#include <vtkPolyDataMapper.h>
#include <vtkGeometryFilter.h>
#include <vtkVolume.h>
#include <vtkPiecewiseFunction.h>
#include <vtkColorTransferFunction.h>
#include <vtkVolumeProperty.h>
#include <vtkSmartVolumeMapper.h>
#include <vtkFixedPointVolumeRayCastMapper.h>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPolyLine.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkHexahedron.h>
#include <vtkUnstructuredGrid.h>


#include <QString>
#include <QFile>

#include "RenderDisplay.h"

class BoundaryDisplay : public RenderDisplay
{
public:

  BoundaryDisplay(long int imax, long int jmax,
                  long int kmax, double delx,
                  double dely, double delz,
                  double o_x, double o_y, double o_z); 
  BoundaryDisplay(long int imax, long int jmax,
                  long int kmax); 


  void drawRectangle(double a_1, double a_2, double a_3,
                     double a_4, double a_5, double a_6, int normal);
  void clearRectangle();

  void update(double delx, double dely, double delz,
              double imax, double jmax, double kmax,
              double o_x,  double o_y,  double o_z);

private:
  vtkSmartPointer<vtkPoints> hexPoints;
  vtkSmartPointer<vtkHexahedron> hex;
  vtkSmartPointer<vtkPolyLine> polyLine;
  vtkSmartPointer<vtkCellArray> cells;
  vtkSmartPointer<vtkPolyData> polyData;
  vtkSmartPointer<vtkPolyDataMapper> polyDataMapper;
  vtkSmartPointer<vtkActor> rectangleActor;
  vtkSmartPointer<vtkUnstructuredGrid> hexGrid;
  vtkSmartPointer<vtkDataSetMapper> dataSetMapper;

  double dx, dy, dz;
};
