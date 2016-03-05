/*
 * MeshDisplay.h
 *
 * Use VTK to display the mesh
 */

#ifndef MESHDISPLAY_H
#define MESHDISPLAY_H

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkMath.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkStructuredGridOutlineFilter.h>
#include <vtkStructuredGridGeometryFilter.h>
#include <vtkProperty.h>
#include <vtkAxesActor.h>
#include <vtkTransform.h>
#include <vtkCaptionActor2D.h>
#include <vtkTextProperty.h>
#include <vtkTextActor.h>
#include <vtkCamera.h>
#include <vtkScalarBarActor.h>

class MeshDisplay {

public:
  MeshDisplay(long int numi, long int numj, long int numk);
  vtkSmartPointer<vtkRenderWindow> getRenderWindow(); 
  vtkSmartPointer<vtkRenderWindowInteractor> getRenderWindowInteractor(); 

  void update(double delx, double dely, double delz,
              double imax, double jmax, double kmax,
              double o_x,  double o_y,  double o_z);

  void AddActor(vtkSmartPointer<vtkActor> &new_actor);
  void RemoveActor(vtkSmartPointer<vtkActor> &new_actor);
  void reset();

  void AddVolume(vtkSmartPointer<vtkVolume> &volume);
  void RemoveVolume(vtkSmartPointer<vtkVolume> &volume);

  void AddScalarBar(vtkSmartPointer<vtkScalarBarActor> &new_actor);
  void RemoveScalarBar(vtkSmartPointer<vtkScalarBarActor> &new_actor);
  
  //vtkSmartPointer<vtkRenderWindowInteractor> &GetInteractor() { return renderWindowInteractor; }
  
  void HideMesh();
  void ShowMesh();

  void HideAxis();
  void ShowAxis();

private:
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;
  vtkSmartPointer<vtkStructuredGrid> structuredGrid;
  vtkSmartPointer<vtkPoints> points;
  vtkSmartPointer<vtkRenderer> renderer;
  vtkSmartPointer<vtkRenderWindow> renderWindow;
  vtkSmartPointer<vtkDataSetMapper> mapper;
  vtkSmartPointer<vtkStructuredGridGeometryFilter> outlineFilter;
  vtkSmartPointer<vtkActor> actor;
  vtkSmartPointer<vtkAxesActor> axes;
  vtkSmartPointer<vtkTransform> transform;

};

#endif
