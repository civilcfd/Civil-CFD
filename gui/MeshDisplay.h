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
#include <vtkGenericOpenGLRenderWindow.h>
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
#include <QVTKInteractor.h>

class MeshDisplay {

public:
  MeshDisplay(long int numi, long int numj, long int numk);
  vtkSmartPointer<vtkGenericOpenGLRenderWindow> getRenderWindow(); 
  vtkSmartPointer<QVTKInteractor> getRenderWindowInteractor(); 

  void completeSetup(); // call this after assigning the render window
  void completeSetup(double delx, double dely, double delz,
                     double imax, double jmax, double kmax, 
                     double o_x,  double o_y,  double o_z); // call this after assigning the render window if an update is also needed

  void update(double delx, double dely, double delz,
              double imax, double jmax, double kmax, 
              double o_x,  double o_y,  double o_z);

  void AddActor(vtkSmartPointer<vtkActor> &new_actor);
  void RemoveActor(vtkSmartPointer<vtkActor> &new_actor);
  void reset();
	void normalizeCamera(int x, double y, double z);

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
  vtkSmartPointer<QVTKInteractor> renderWindowInteractor;
  vtkSmartPointer<vtkStructuredGrid> structuredGrid;
  vtkSmartPointer<vtkPoints> points;
  vtkSmartPointer<vtkDataSetMapper> mapper;
  vtkSmartPointer<vtkStructuredGridGeometryFilter> outlineFilter;
  vtkSmartPointer<vtkActor> actor;
  vtkSmartPointer<vtkAxesActor> axes;
  vtkSmartPointer<vtkTransform> transform;

  vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow;
  vtkSmartPointer<vtkRenderer> renderer;
};

#endif
