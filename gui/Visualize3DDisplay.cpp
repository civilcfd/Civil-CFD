/*
 * Visualize3DDisplay.cpp
 *
 * Implements VTK displayer for the redering window
 *
 */
#include "Visualize3DDisplay.h"
/*
Visualize3DDisplay::Visualize3DDisplay(long int imax, long int jmax, long int kmax, double delx, double dely, double delz, double ox, double oy, double oz) : 
  VisualizeDisplay (imax,jmax,kmax,delx,dely,delz,ox,oy,oz) {


  //update(delx,dely,delz,imax,jmax,kmax,ox,oy,oz);
  //VTKmapper = NULL;
  //VTKactor = NULL;
  iso3dActor = NULL;
}*/

Visualize3DDisplay::Visualize3DDisplay(long int imax, long int jmax, long int kmax) : 
  VisualizeDisplay(imax,jmax,kmax) {
  //VTKmapper = NULL;
  //VTKactor = NULL;
  iso3dActor = NULL;
}
 
void Visualize3DDisplay::block(QString vtkFile, int normal, double origin, double del) {
  
  RemoveVolume(volume);
  
  if(!QFile::exists(vtkFile)) return;

	volReader = NULL;
  volReader = vtkSmartPointer<vtkXMLImageDataReader>::New();
  volReader->SetFileName(vtkFile.toStdString().c_str());
  volReader->Update();

  vtkSmartPointer<vtkPlane> planePos = vtkSmartPointer<vtkPlane>::New();
  vtkSmartPointer<vtkPlane> planeNeg = vtkSmartPointer<vtkPlane>::New();

  switch(normal) {
  case 0:
    planePos->SetOrigin(origin,0,0);
    planePos->SetNormal(1,0,0);
    planeNeg->SetOrigin(origin+del,0,0);
    planeNeg->SetNormal(-1,0,0);
    break;
  case 1:
    planePos->SetOrigin(0,origin,0);
    planePos->SetNormal(0,1,0);
    planeNeg->SetOrigin(0,origin+del,0);
    planeNeg->SetNormal(0,-1,0);
    break;
  case 2:
    planePos->SetOrigin(0,0,origin);
    planePos->SetNormal(0,0,1);
    planeNeg->SetOrigin(0,0,origin+del);
    planeNeg->SetNormal(0,0,-1);
    break;
  }

  vtkSmartPointer<vtkVolumeProperty> volumeProperty = 
    vtkSmartPointer<vtkVolumeProperty>::New();
//  volumeProperty->ShadeOff();
  volumeProperty->SetInterpolationType(VTK_LINEAR_INTERPOLATION);

  vtkSmartPointer<vtkPiecewiseFunction> compositeOpacity = 
    vtkSmartPointer<vtkPiecewiseFunction>::New();
  compositeOpacity->AddPoint(0.0,1.0);
  compositeOpacity->AddPoint(0.5,1.0);
  compositeOpacity->AddPoint(1.0,0.0);
  volumeProperty->SetScalarOpacity(compositeOpacity); // composite first.
 
  vtkSmartPointer<vtkColorTransferFunction> color = 
    vtkSmartPointer<vtkColorTransferFunction>::New();
  color->AddRGBPoint(0.0,0.5,0.5,0.5);
  color->AddRGBPoint(0.5,0.5,0.5,0.5);
  color->AddRGBPoint(1.0,1.0,1.0,1.0);
  volumeProperty->SetColor(color);
 
  volume = 
    vtkSmartPointer<vtkVolume>::New();
  volume->SetProperty(volumeProperty);

//  volumeMapper = vtkSmartPointer<vtkFixedPointVolumeRayCastMapper>::New();
  volMapper = vtkSmartPointer<vtkSmartVolumeMapper>::New();
  volMapper->SetInputConnection(volReader->GetOutputPort());
  volMapper->AddClippingPlane(planePos);
  volMapper->AddClippingPlane(planeNeg);
  //volMapper->SetRequestedRenderModeToRayCastAndTexture();
//  volumeMapper->SetBlendModeToComposite();
  volume->SetMapper(volMapper);
  volume->Update();

//  vtkSmartPointer<vtkCutter> cutter = vtkSmartPointer<vtkCutter>::New();
//  cutter->SetInputConnection(volReader->GetOutputPort());
//  cutter->SetCutFunction(plane);
 
  AddVolume(volume);

  //reset();
}

void Visualize3DDisplay::contour3d(QString vtkFile, int normal, double origin, double del) {
  
  if(iso3dActor != NULL) RemoveActor(iso3dActor);
  
  if(!QFile::exists(vtkFile)) return;

	iso3dReader = NULL;
  iso3dReader = vtkSmartPointer<vtkXMLImageDataReader>::New();
  iso3dReader->SetFileName(vtkFile.toStdString().c_str());
  iso3dReader->Update();

  vtkSmartPointer<vtkPlane> planePos = vtkSmartPointer<vtkPlane>::New();
  vtkSmartPointer<vtkPlane> planeNeg = vtkSmartPointer<vtkPlane>::New();

  switch(normal) {
  case 0:
    planePos->SetOrigin(origin,0,0);
    planePos->SetNormal(1,0,0);
    planeNeg->SetOrigin(origin+del,0,0);
    planeNeg->SetNormal(-1,0,0);
    break;
  case 1:
    planePos->SetOrigin(0,origin,0);
    planePos->SetNormal(0,1,0);
    planeNeg->SetOrigin(0,origin+del,0);
    planeNeg->SetNormal(0,-1,0);
    break;
  case 2:
    planePos->SetOrigin(0,0,origin);
    planePos->SetNormal(0,0,1);
    planeNeg->SetOrigin(0,0,origin+del);
    planeNeg->SetNormal(0,0,-1);
    break;
  }
 
 	iso3dFilter = NULL;
 	iso3dFilter = vtkSmartPointer<vtkContourFilter>::New();
 	iso3dFilter->SetInputConnection(iso3dReader->GetOutputPort());
 	iso3dFilter->SetValue(0, 0.5);
 	
 	iso3dMapper = NULL;
  iso3dMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  iso3dMapper->SetInputConnection(iso3dFilter->GetOutputPort());
  iso3dMapper->AddClippingPlane(planePos);
  iso3dMapper->AddClippingPlane(planeNeg);
  iso3dMapper->ScalarVisibilityOff();
  
  iso3dActor = NULL;
  iso3dActor = vtkSmartPointer<vtkActor>::New();
  iso3dActor->SetMapper(iso3dMapper);
  iso3dActor->GetProperty()->SetColor(0,0,1); 
  
  AddActor(iso3dActor);

}

void Visualize3DDisplay::hideContour3d() {
  RemoveVolume(volume3d);
}


