/*
 * RenderDisplay.cpp
 *
 * Implements VTK displayer for the redering window
 *
 */
#include "RenderDisplay.h"

RenderDisplay::RenderDisplay(long int imax, long int jmax, long int kmax, double delx, double dely, double delz, double ox, double oy, double oz) : 
  GeometryDisplay(imax,jmax,kmax) {


  //update(delx,dely,delz,imax,jmax,kmax,ox,oy,oz);
}

RenderDisplay::RenderDisplay(long int imax, long int jmax, long int kmax) : 
  GeometryDisplay(imax,jmax,kmax) {

}

void RenderDisplay::connectVTK(QString vtkFile) {
  
  RemoveVolume(volume);
  volume = NULL;

  if(!QFile::exists(vtkFile)) {
    getRenderWindow()->Render();
    return;
  }

	reader = NULL;
  reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
  reader->SetFileName(vtkFile.toStdString().c_str());
  reader->Update();

  vtkSmartPointer<vtkVolumeProperty> volumeProperty = 
    vtkSmartPointer<vtkVolumeProperty>::New();
  //volumeProperty->ShadeOff();
  volumeProperty->SetInterpolationType(VTK_LINEAR_INTERPOLATION);
  vtkSmartPointer<vtkPiecewiseFunction> compositeOpacity = 
    vtkSmartPointer<vtkPiecewiseFunction>::New();
  compositeOpacity->AddPoint(0.0,0.0);
  compositeOpacity->AddPoint(1.0,0.99);
  volumeProperty->SetScalarOpacity(compositeOpacity); // composite first.
 
  vtkSmartPointer<vtkColorTransferFunction> color = 
    vtkSmartPointer<vtkColorTransferFunction>::New();
  color->AddRGBPoint(0.0,0.0,0.0,1.0);
  color->AddRGBPoint(0.5,1.0,1.0,1.0);
  color->AddRGBPoint(1.0,1.0,0.0,0.0);
  volumeProperty->SetColor(color);
 
  volume = 
    vtkSmartPointer<vtkVolume>::New();
  volume->SetProperty(volumeProperty);

//  volumeMapper = vtkSmartPointer<vtkFixedPointVolumeRayCastMapper>::New();
  volumeMapper = vtkSmartPointer<vtkSmartVolumeMapper>::New();
  volumeMapper->SetInputConnection(reader->GetOutputPort());
  volume->SetMapper(volumeMapper);
  volume->Update();

/*
  vtkSmartPointer<vtkOpenGLExtensionManager> extManager = vtkSmartPointer<vtkOpenGLExtensionManager>::New();
*/

  AddVolume(volume);
/*
  getRenderWindow()->Render();

  volumeMapper->SetRequestedRenderModeToRayCastAndTexture();
  getRenderWindow()->Render();
*/
  //volumeMapper->SetRequestedRenderModeToRayCast();
/*
  extManager->SetRenderWindow(getRenderWindow());
  extManager->LoadExtension("GL_VERSION_4_0");
  extManager->LoadExtension("GL_ARB_compatibility");
  extManager->LoadExtension("GL_ARB_multitexture");
  */
  reset();
}
