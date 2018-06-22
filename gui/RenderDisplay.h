/* 
 * RenderDisplay.h
 *
 * Descibe a class that shows a mesh with a VTK file
 *
 * Inherits from MeshDisplay.h
 */

#ifndef RENDER_DISPLAY_H
#define RENDER_DISPLAY_H

#include <vtkDataReader.h>
#include <vtkXMLImageDataReader.h>
#include <vtkPolyDataMapper.h>
#include <vtkGeometryFilter.h>
#include <vtkVolume.h>
#include <vtkPiecewiseFunction.h>
#include <vtkColorTransferFunction.h>
#include <vtkVolumeProperty.h>
#include <vtkSmartVolumeMapper.h>
#include <vtkFixedPointVolumeRayCastMapper.h>
#include <QString>
#include <QFile>

#include "GeometryDisplay.h"

class RenderDisplay : public GeometryDisplay
{
public:

  RenderDisplay(long int imax, long int jmax,
                  long int kmax, double delx,
                  double dely, double delz,
                  double ox, double oy, double oz); 
  RenderDisplay(long int imax, long int jmax,
                  long int kmax); 


  void connectVTK(QString vtkFile);


private:
  vtkSmartPointer<vtkXMLImageDataReader> reader;
  vtkSmartPointer<vtkPolyDataMapper> VTKmapper;

//  vtkSmartPointer<vtkFixedPointVolumeRayCastMapper> volumeMapper;
  vtkSmartPointer<vtkSmartVolumeMapper> volumeMapper;

  vtkSmartPointer<vtkVolume> volume;


};

#endif
