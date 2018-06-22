/* 
 * VisualizeDisplay.h
 *
 * Descibe a class that shows a mesh with a VTK file
 *
 * Inherits from MeshDisplay.h
 */

#ifndef VISUALIZE_DISPLAY_H
#define VISUALIZE_DISPLAY_H

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

#include "GeometryDisplay.h"

class VisualizeDisplay : public MeshDisplay //GeometryDisplay
{
public:

  //VisualizeDisplay(long int imax, long int jmax,
  //                long int kmax, double delx,
  //                double dely, double delz, 
  //                double ox, double oy, double oz); 
  VisualizeDisplay(long int imax, long int jmax,
                  long int kmax); 


  virtual void block(QString vtkFile, int normal, double origin, double del);
  void clip(QString vtkFile, int normal, double origin);
  void clipVector(QString vtkFile, int normal, double origin, bool normalize);
  void vector(QString vtkFile, int normal, double origin);
  void hideVector();
  void hideBlock();
  void setRange(double a, double b);
  void getRange(double &a);
  void hideLegend();
  void showLegend();
  void clear();

protected:
  vtkSmartPointer<vtkXMLImageDataReader> reader;
  vtkSmartPointer<vtkPolyDataMapper> VTKmapper;
  vtkSmartPointer<vtkActor> VTKactor;
  vtkSmartPointer<vtkXMLImageDataReader> vectReader;
  vtkSmartPointer<vtkGeometryFilter> geometryFilter;
  vtkSmartPointer<vtkGeometryFilter> vectGeometryFilter;
  vtkSmartPointer<vtkPolyDataMapper> vectMapper;
  vtkSmartPointer<vtkActor> vectActor;
  vtkSmartPointer<vtkXMLImageDataReader> volReader;
  vtkSmartPointer<vtkGeometryFilter> volGeometryFilter;
	vtkSmartPointer<vtkScalarBarActor> legend;
//  vtkSmartPointer<vtkFixedPointVolumeRayCastMapper> volumeMapper;
  vtkSmartPointer<vtkSmartVolumeMapper> volMapper;

  vtkSmartPointer<vtkVolume> volume;
  //vtkSmartPointer<vtkActor> edgeLineActor;


};

#endif
