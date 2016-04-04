/*
 * MeshDisplay.cpp
 *
 * Implements mesh renderWindow using VTK
 */

#include "MeshDisplay.h"

void MeshDisplay::update(double delx, double dely, double delz,
                         double imax, double jmax, double kmax,
                         double o_x,  double o_y,  double o_z) {
  //structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
 
  //structuredGrid->SetDimensions(imax * delx,jmax * dely,kmax * delz);  
  structuredGrid->SetDimensions(imax,jmax,kmax);
  structuredGrid->SetUpdateExtent(0,imax-1,0,jmax-1,0,kmax-1);

  points = NULL;
  points = vtkSmartPointer<vtkPoints>::New();


  //for(double k = 0; k < kmax; k++)  {
  //  for(double  j = 0; j < jmax; j++)  {
  //    for(double  i = 0; i < imax; i++)  {
  //      points->InsertNextPoint(i * delx, j * dely, k * delz);
  //    }
  //  }
  //}
  for(unsigned int k = 0; k < kmax; k++)  {
    for(unsigned int j = 0; j < jmax; j++)  {
      for(unsigned int i = 0; i < imax; i++)  {
        points->InsertNextPoint(i*delx + o_x, j*dely + o_y, k*delz + o_z);
      }
    }
  }


  structuredGrid->SetPoints(points);
  structuredGrid->Update();
  // Create a mapper and actor
  /* mapper = vtkSmartPointer<vtkDataSetMapper>::New();
  mapper->SetInputConnection(structuredGrid->GetProducerPort());
  actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  actor->GetProperty()->SetRepresentationToWireframe();

  // Visualize
  renderer = vtkSmartPointer<vtkRenderer>::New();
  renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  //renderWindowInteractor->SetRenderWindow(renderWindow); 
  renderer->AddActor(actor);
  //renderer->SetBackground(1,1,1); // Background color white
*/
  reset();
}

void MeshDisplay::reset() {
  renderer->ResetCamera();
}

vtkSmartPointer<vtkRenderWindow> MeshDisplay::getRenderWindow() {
  return renderWindow;
}

vtkSmartPointer<vtkRenderWindowInteractor> MeshDisplay::getRenderWindowInteractor() {
  return renderWindowInteractor;
}

void MeshDisplay::AddActor(vtkSmartPointer<vtkActor> &new_actor) {
  renderer->AddActor(new_actor);
}

void MeshDisplay::AddScalarBar(vtkSmartPointer<vtkScalarBarActor> &new_actor) {
  renderer->AddActor2D(new_actor);
}

void MeshDisplay::RemoveScalarBar(vtkSmartPointer<vtkScalarBarActor> &new_actor) {
  renderer->RemoveActor2D(new_actor);
}

void MeshDisplay::RemoveActor(vtkSmartPointer<vtkActor> &new_actor) {
  renderer->RemoveActor(new_actor);
}

void MeshDisplay::AddVolume(vtkSmartPointer<vtkVolume> &new_volume) {
  renderer->AddVolume(new_volume);
}

void MeshDisplay::HideMesh() {
  RemoveActor(actor);
}

void MeshDisplay::ShowMesh() {
  AddActor(actor);
}
void MeshDisplay::HideAxis() {
  renderer->RemoveActor(axes);
}

void MeshDisplay::ShowAxis() {
  renderer->AddActor(axes);
}
void MeshDisplay::RemoveVolume(vtkSmartPointer<vtkVolume> &new_volume) {
  renderer->RemoveVolume(new_volume);
}
/*
MeshDisplay::MeshDisplay(long int imax, long int jmax, long int kmax, double delx, double dely, double delz) {

  MeshDisplay(imax,jmax,kmax);
  update(imax,jmax,kmax,delz,dely,delz);

}*/

MeshDisplay::MeshDisplay(long int numi, long int numj, long int numk) {

  structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
 
  points = vtkSmartPointer<vtkPoints>::New();
 
  for(unsigned int k = 0; k < numk; k++)  {
    for(unsigned int j = 0; j < numj; j++)  {
      for(unsigned int i = 0; i < numi; i++)  {
        points->InsertNextPoint(i, j, k);
      }
    }
  }

  //specify the dimensions of the grid
  structuredGrid->SetDimensions(numi, numj, numk);
  structuredGrid->SetPoints(points);
 
  std::cout << "There are " << structuredGrid->GetNumberOfPoints() << " points." << std::endl;
  std::cout << "There are " << structuredGrid->GetNumberOfCells() << " cells." << std::endl;
 
  outlineFilter = vtkSmartPointer<vtkStructuredGridGeometryFilter>::New();
  outlineFilter->SetInputConnection(structuredGrid->GetProducerPort());
  outlineFilter->Update();
 
  // Create a mapper and actor
  mapper = vtkSmartPointer<vtkDataSetMapper>::New();
  mapper->SetInputConnection(structuredGrid->GetProducerPort());
  actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  actor->GetProperty()->SetRepresentationToWireframe();

  transform = vtkSmartPointer<vtkTransform>::New();
  transform->Scale(0.8, 0.8, 0.8);
 
  axes = vtkSmartPointer<vtkAxesActor>::New();
 
  // The axes are positioned with a user transform
  axes->SetUserTransform(transform);
 
  // properties of the axes labels can be set as follows
  // this sets the x axis label to red
  //axes->GetXAxisCaptionActor2D()->GetCaptionTextProperty()->SetFontSize(2);
  axes->GetXAxisCaptionActor2D()->GetTextActor()->SetTextScaleModeToViewport();   
  axes->GetYAxisCaptionActor2D()->GetTextActor()->SetTextScaleModeToViewport();   
  axes->GetZAxisCaptionActor2D()->GetTextActor()->SetTextScaleModeToViewport();   

  // the actual text of the axis label can be changed:
  // axes->SetXAxisLabelText("test");

  // Visualize
  renderer = vtkSmartPointer<vtkRenderer>::New();
  renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  //renderWindowInteractor->SetRenderWindow(renderWindow); 
  renderer->AddActor(actor);
  renderer->AddActor(axes);

  renderer->GetActiveCamera()->SetPosition(1,-1,1);
  renderer->GetActiveCamera()->SetFocalPoint(0,0,0);
  renderer->GetActiveCamera()->SetViewUp(0,0,1);
  renderer->ResetCamera();

 renderer->SetBackground(0.4,0.4,0.4); 

  //renderWindow->Render();
}

void MeshDisplay::normalizeCamera(int x, double y, double z) {
	switch(x) {
	case 1:
		renderer->GetActiveCamera()->SetPosition(10,y/2,z/2);
		renderer->GetActiveCamera()->SetViewUp(0,0,1);
		renderer->GetActiveCamera()->SetFocalPoint(0,y/2,z/2);
		break;
		
	case 2:
		renderer->GetActiveCamera()->SetPosition(y/2,-10,z/2);
		renderer->GetActiveCamera()->SetViewUp(0,0,1);
		renderer->GetActiveCamera()->SetFocalPoint(y/2,0,z/2);
		break;
		
	case 3:
		renderer->GetActiveCamera()->SetPosition(y/2,z/2,10);
		renderer->GetActiveCamera()->SetViewUp(0,1,0);
		renderer->GetActiveCamera()->SetFocalPoint(y/2,z/2,0);
		break;
	}
	
	renderer->ResetCamera();
}