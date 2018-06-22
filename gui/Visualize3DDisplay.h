/* 
 * Visualize3DDisplay.h
 *
 * Descibe a class that shows a mesh with a VTK file
 *
 * Inherits from MeshDisplay.h
 */

#ifndef VISUALIZE3D_DISPLAY_H
#define VISUALIZE3D_DISPLAY_H

#include <vtkDataReader.h>
#include <vtkXMLImageDataReader.h>
#include <vtkStructuredPoints.h>
#include <vtkCutter.h>
#include <vtkPolyDataMapper.h>
#include <vtkStructuredGridGeometryFilter.h>
#include <vtkGeometryFilter.h>
#include <vtkVolume.h>
#include <vtkLookupTable.h>
#include <vtkPlane.h>
#include <vtkPiecewiseFunction.h>
#include <vtkColorTransferFunction.h>
#include <vtkVolumeProperty.h>
#include <vtkSmartVolumeMapper.h>
#include <vtkFixedPointVolumeRayCastMapper.h>
#include <vtkContourFilter.h>
#include <vtkArrowSource.h>
#include <vtkGlyph3D.h>
#include <vtkScalarBarActor.h>
#include <vtkScalarBarWidget.h>
#include <vtkVectorNorm.h>
#include <QString>
#include <QFile>

#include "VisualizeDisplay.h"
#include "GeometryDisplay.h"

class Visualize3DDisplay : public VisualizeDisplay //GeometryDisplay
{
public:

  //Visualize3DDisplay(long int imax, long int jmax,
  //                long int kmax, double delx,
  //                double dely, double delz, 
  //                double ox, double oy, double oz); 
  Visualize3DDisplay(long int imax, long int jmax,
                  long int kmax); 
                  
  void hideContour3d();
  void contour3d(QString vtkFile, int normal, double origin, double del);
  void block(QString vtkFile, int normal, double origin, double del);

private:

  vtkSmartPointer<vtkVolume> volume3d;
  vtkSmartPointer<vtkXMLImageDataReader> vol3dReader;
  vtkSmartPointer<vtkSmartVolumeMapper> vol3dMapper;
  vtkSmartPointer<vtkXMLImageDataReader> iso3dReader;
  vtkSmartPointer<vtkContourFilter> iso3dFilter;
  vtkSmartPointer<vtkPolyDataMapper> iso3dMapper;
  vtkSmartPointer<vtkActor> iso3dActor;

};

#endif
