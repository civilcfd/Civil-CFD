/* 
 * GeometryDisplay.h
 *
 * Descibe a class that shows a mesh with an STL
 * cutting it
 *
 * Inherits from MeshDisplay.h
 */

#ifndef GEOMETRYDISPLAY_H
#define GEOMETRYDISPLAY_H

#include <vtkSTLReader.h>
#include <vtkPolyDataMapper.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <QString>
#include <QFile>

#include "MeshDisplay.h"

class GeometryDisplay : public MeshDisplay
{
public:

  GeometryDisplay(long int imax, long int jmax,
                  long int kmax); 

  void connectSTL(QString stlFile);
  void drawPoint(double x, double y, double z);
  void clearPoint();

private:
  vtkSmartPointer<vtkSTLReader> reader;
  vtkSmartPointer<vtkPolyDataMapper> STLmapper;
  vtkSmartPointer<vtkActor> STLactor;

  vtkSmartPointer<vtkActor> pointActor;
};

#endif
