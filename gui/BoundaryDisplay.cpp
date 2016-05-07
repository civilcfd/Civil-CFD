/*
 * BoundaryDisplay.cpp
 *
 * Implements VTK displayer for the redering window
 *
 */
#include "BoundaryDisplay.h"

void BoundaryDisplay::update(double delx, double dely, double delz,
              double imax, double jmax, double kmax,
              double o_x,  double o_y,  double o_z) {

  dx = delx;
  dy = dely;
  dz = delz;
  

  MeshDisplay::update(delx, dely, delz, imax, jmax, kmax, o_x, o_y, o_z);
}
              
BoundaryDisplay::BoundaryDisplay(long int imax, long int jmax, long int kmax, double delx, double dely, double delz, double o_x, double o_y, double o_z) : 
  RenderDisplay(imax,jmax,kmax,delx,dely,delz,o_x,o_y,o_z) {

  dx = delx;
  dy = dely;
  dz = delz;
  

}

BoundaryDisplay::BoundaryDisplay(long int imax, long int jmax, long int kmax) : 
  RenderDisplay(imax,jmax,kmax) {

}

void BoundaryDisplay::clearRectangle() {

  RemoveActor(rectangleActor);
  hex = NULL;
  hexPoints = NULL;
  hexGrid = NULL;
  dataSetMapper = NULL;
  rectangleActor = NULL;
}

#define set_vector(a,b,c,d) { a[0] = b; a[1] = c; a[2] = d; }

void BoundaryDisplay::drawRectangle(double a_1, double a_2, double a_3,
                     double b_1, double b_2, double b_3, int normal) {

  clearRectangle();

  double norm[3][3] = { { dx, 0, 0 },
                        { 0, dy, 0 },
                        { 0, 0, dz } };

  // Create five points.
  double p0[3];
  double p1[3];
  double p2[3];
  double p3[3];
  double p4[3];
  double p5[3];
  double p6[3];
  double p7[3];
 
  switch(normal) {
  case 2: // z-normal
    set_vector(p0, a_1, a_2, a_3);
    set_vector(p1, b_1, a_2, a_3);
    set_vector(p2, b_1, b_2, a_3);
    set_vector(p3, a_1, b_2, a_3);
    break;
  case 1: // y-normal
    set_vector(p0, a_1, a_2, a_3);
    set_vector(p1, b_1, a_2, a_3);
    set_vector(p2, b_1, a_2, b_3);
    set_vector(p3, a_1, a_2, b_3);
    break;
  case 0: // x-normal
    set_vector(p0, a_1, a_2, a_3);
    set_vector(p1, a_1, b_2, a_3);
    set_vector(p2, a_1, b_2, b_3);
    set_vector(p3, a_1, a_2, b_3);
    break;
  }
  
  set_vector(p4, p0[0] + norm[normal][0], 
                 p0[1] + norm[normal][1], p0[2] + norm[normal][2]);
  set_vector(p5, p1[0] + norm[normal][0], 
                 p1[1] + norm[normal][1], p1[2] + norm[normal][2]);
  set_vector(p6, p2[0] + norm[normal][0], 
                 p2[1] + norm[normal][1], p2[2] + norm[normal][2]);
  set_vector(p7, p3[0] + norm[normal][0],
                 p3[1] + norm[normal][1], p3[2] + norm[normal][2]);  
 
  // Create a vtkPoints object and store the points in it
  hexPoints =vtkSmartPointer<vtkPoints>::New();
  hexPoints->InsertPoint(0, p0);
  hexPoints->InsertPoint(1, p1);
  hexPoints->InsertPoint(2, p2);
  hexPoints->InsertPoint(3, p3);
  hexPoints->InsertPoint(4, p4);
  hexPoints->InsertPoint(5, p5);
  hexPoints->InsertPoint(6, p6);
  hexPoints->InsertPoint(7, p7);
  
  /*points->InsertNextPoint(origin);
  points->InsertNextPoint(p0);
  points->InsertNextPoint(p1);
  points->InsertNextPoint(p2);
  points->InsertNextPoint(p3);*/
  
  hex = vtkSmartPointer<vtkHexahedron>::New();
  hex->GetPointIds()->SetId(0, 0);
  hex->GetPointIds()->SetId(1, 1);
  hex->GetPointIds()->SetId(2, 2);
  hex->GetPointIds()->SetId(3, 3);
  hex->GetPointIds()->SetId(4, 4);
  hex->GetPointIds()->SetId(5, 5);
  hex->GetPointIds()->SetId(6, 6);
  hex->GetPointIds()->SetId(7, 7);
  
  hexGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  hexGrid->Allocate(1, 1);
  hexGrid->InsertNextCell(hex->GetCellType(),
                          hex->GetPointIds());
  hexGrid->SetPoints(hexPoints);  
 
  dataSetMapper = vtkSmartPointer<vtkDataSetMapper>::New();
  dataSetMapper->SetInput(hexGrid); 
 
//  polyLine = vtkSmartPointer<vtkPolyLine>::New();
//  polyLine->GetPointIds()->SetNumberOfIds(5);
//  for(unsignedg int i = 0; i < 5; i++)
//    {
//    polyLine->GetPointIds()->SetId(i,i);
//    }
 
  // Create a cell array to store the lines in and add the lines to it
//  cells = vtkSmartPointer<vtkCellArray>::New();
//  cells->InsertNextCell(polyLine);
 
  // Create a polydata to store everything in
//  polyData = vtkSmartPointer<vtkPolyData>::New();
 
  // Add the points to the dataset
//  polyData->SetPoints(points);
 
  // Add the lines to the dataset
//  polyData->SetLines(cells);
 
  // Setup actor and mapper
//  polyDataMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//  polyDataMapper->SetInput(polyData);

 
  rectangleActor = vtkSmartPointer<vtkActor>::New();
  rectangleActor->SetMapper(dataSetMapper);
  
  rectangleActor->GetProperty()->SetDiffuseColor(1,1,0);
  rectangleActor->GetProperty()->SetOpacity(1.0);
 
  AddActor(rectangleActor);
 
}
