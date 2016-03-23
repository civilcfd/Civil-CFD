/*
 * GeometryDisplay.cpp
 *
 * Implements the STL / mesh display
 */

#include "GeometryDisplay.h"

GeometryDisplay::GeometryDisplay(long int imax, long int jmax, long int kmax) : 
  MeshDisplay(imax,jmax,kmax) {

  reader = vtkSmartPointer<vtkSTLReader>::New();
 
   
}

void GeometryDisplay::drawPoint(double x, double y, double z) {
  // Create the geometry of a point (the coordinate)
  vtkSmartPointer<vtkPoints> points =
    vtkSmartPointer<vtkPoints>::New();
  const double p[3] = {x, y, z};
 
  // Create the topology of the point (a vertex)
  vtkSmartPointer<vtkCellArray> vertices =
    vtkSmartPointer<vtkCellArray>::New();
  vtkIdType pid[1];
  pid[0] = points->InsertNextPoint(p);
  vertices->InsertNextCell(1,pid);
 
  // Create a polydata object
  vtkSmartPointer<vtkPolyData> point =
    vtkSmartPointer<vtkPolyData>::New();
 
  // Set the points and vertices we created as the geometry and topology of the polydata
  point->SetPoints(points);
  point->SetVerts(vertices);
 
  // Visualize
  vtkSmartPointer<vtkPolyDataMapper> mapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInput(point);
 
  pointActor =
    vtkSmartPointer<vtkActor>::New();
  pointActor->SetMapper(mapper);
  pointActor->GetProperty()->SetPointSize(20);
  pointActor->GetProperty()->SetDiffuseColor(1,0,0);

  AddActor(pointActor);

  reset();
}

void GeometryDisplay::clearPoint() {
  RemoveActor(pointActor);
}

void GeometryDisplay::connectSTL(QString stlFile) {

	RemoveActor(STLactor);
	STLactor = NULL;
	
  if(!QFile::exists(stlFile)) return;

  reader->SetFileName(stlFile.toStdString().c_str());
  reader->Update();
 
  // Visualize

  STLmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  STLmapper->SetInputConnection(reader->GetOutputPort());
 
  STLactor = vtkSmartPointer<vtkActor>::New();
  STLactor->SetMapper(STLmapper);
  STLactor->GetProperty()->SetDiffuseColor(0,0,1);

  AddActor(STLactor);
  reset();
}
