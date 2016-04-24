/*
 * VisualizeDisplay.cpp
 *
 * Implements VTK displayer for the redering window
 *
 */
#include "VisualizeDisplay.h"

VisualizeDisplay::VisualizeDisplay(long int imax, long int jmax, long int kmax, double delx, double dely, double delz, double ox, double oy, double oz) : 
  MeshDisplay (imax,jmax,kmax) {


  update(delx,dely,delz,imax,jmax,kmax,ox,oy,oz);
  VTKmapper = NULL;
  VTKactor = NULL;
}

VisualizeDisplay::VisualizeDisplay(long int imax, long int jmax, long int kmax) : 
  MeshDisplay(imax,jmax,kmax) {
  VTKmapper = NULL;
  VTKactor = NULL;
} 

void VisualizeDisplay::block(QString vtkFile, int normal, double origin, double del) {
  
  RemoveVolume(volume);
//  RemoveActor(edgeLineActor);
  
  if(!QFile::exists(vtkFile)) return;

	volReader = NULL;
  volReader = vtkSmartPointer<vtkStructuredPointsReader>::New();
  volReader->SetFileName(vtkFile.toStdString().c_str());
  volReader->Update();

  vtkSmartPointer<vtkPlane> planePos = vtkSmartPointer<vtkPlane>::New();
  vtkSmartPointer<vtkPlane> planeNeg = vtkSmartPointer<vtkPlane>::New();

  switch(normal) {
  case 0:
    planePos->SetOrigin(origin-del,0,0);
    planePos->SetNormal(1,0,0);
    planeNeg->SetOrigin(origin+del,0,0);
    planeNeg->SetNormal(-1,0,0);
    break;
  case 1:
    planePos->SetOrigin(0,origin-del,0);
    planePos->SetNormal(0,1,0);
    planeNeg->SetOrigin(0,origin+del,0);
    planeNeg->SetNormal(0,-1,0);
    break;
  case 2:
    planePos->SetOrigin(0,0,origin-del);
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
  color->AddRGBPoint(0.0  ,0.5,0.5,0.5);
  color->AddRGBPoint(0.9	,0.5,0.5,0.5);
  color->AddRGBPoint(1.0	,0,0,0);
  volumeProperty->SetColor(color);
 
  volume = 
    vtkSmartPointer<vtkVolume>::New();
  volume->SetProperty(volumeProperty);
	 

//  volumeMapper = vtkSmartPointer<vtkFixedPointVolumeRayCastMapper>::New();
  volMapper = vtkSmartPointer<vtkSmartVolumeMapper>::New();
  volMapper->SetInputConnection(volReader->GetOutputPort());
  volMapper->AddClippingPlane(planePos);
  volMapper->AddClippingPlane(planeNeg);
//  volumeMapper->SetBlendModeToComposite();
  volume->SetMapper(volMapper);
  volume->Update();

//  vtkSmartPointer<vtkCutter> cutter = vtkSmartPointer<vtkCutter>::New();
//  cutter->SetInputConnection(volReader->GetOutputPort());
//  cutter->SetCutFunction(plane);
 
  AddVolume(volume);


	/*vtkSmartPointer<vtkContourFilter> edgeLine = vtkSmartPointer<vtkContourFilter>::New();
	edgeLine->SetInputConnection(volReader->GetOutputPort());
	edgeLine->SetValue(0,0.01);
	
  vtkSmartPointer<vtkPolyDataMapper> edgeLineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	edgeLineActor = NULL;
	edgeLineActor = vtkSmartPointer<vtkActor>::New();
 
  edgeLineMapper->SetInputConnection(edgeLine->GetOutputPort());
  edgeLineMapper->AddClippingPlane(planePos);
  edgeLineMapper->AddClippingPlane(planeNeg);
  edgeLineActor->SetMapper(edgeLineMapper);
  edgeLineActor->GetProperty()->SetLineWidth(4);
  edgeLineActor->GetProperty()->SetColor(0,0,0); 
  
  AddActor(edgeLineActor); */

	
  //reset();
}

void VisualizeDisplay::hideBlock() {
  RemoveVolume(volume);
//  RemoveActor(edgeLineActor);
}
void VisualizeDisplay::hideLegend() {
  	RemoveScalarBar(legend);
}
void VisualizeDisplay::showLegend() {
  	AddScalarBar(legend);
}

void VisualizeDisplay::clip(QString vtkFile, int normal, double origin) {
	
	reader = NULL;
 	reader = vtkSmartPointer<vtkStructuredPointsReader>::New();
  reader->SetFileName(vtkFile.toStdString().c_str());
  reader->Update();

	geometryFilter = NULL;
  geometryFilter =
    vtkSmartPointer<vtkGeometryFilter>::New();
  //geometryFilter->SetInputConnection(reader->GetOutputPort());

  vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();

  switch(normal) {
  case 0:
    plane->SetOrigin(origin,0,0);
    plane->SetNormal(1,0,0);
    break;
  case 1:
    plane->SetOrigin(0,origin,0);
    plane->SetNormal(0,1,0);
    break;
  case 2:
    plane->SetOrigin(0,0,origin);
    plane->SetNormal(0,0,1);
    break;
  }

  vtkSmartPointer<vtkCutter> cutter = vtkSmartPointer<vtkCutter>::New();
  cutter->SetInputConnection(reader->GetOutputPort());
  cutter->SetCutFunction(plane);

  vtkSmartPointer<vtkContourFilter> iso = vtkSmartPointer<vtkContourFilter>::New();
  iso->SetInputConnection(reader->GetOutputPort());
  iso->SetValue(0, 0.5);

  geometryFilter->SetInputConnection(cutter->GetOutputPort());

  vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
  /* lut->SetNumberOfColors(64); */
  lut->SetHueRange(0.667, 0.0);
  lut->Build();

  RemoveActor(VTKactor);

	VTKmapper = NULL;
	VTKactor = NULL;
  VTKmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  VTKactor = vtkSmartPointer<vtkActor>::New();
 
  VTKmapper->SetInputConnection(geometryFilter->GetOutputPort());
  VTKmapper->SetLookupTable(lut);
  //VTKmapper->UseLookupTableScalarRangeOff();
  VTKmapper->SetScalarRange(reader->GetOutput()->GetScalarRange());

  VTKactor->SetMapper(VTKmapper);

  AddActor(VTKactor);

	RemoveScalarBar(legend);
	legend = NULL;
	legend = vtkSmartPointer<vtkScalarBarActor>::New();
	legend->SetLookupTable(VTKmapper->GetLookupTable());
	legend->SetWidth(0.07);
	legend->GetLabelTextProperty()->SetFontSize(9);
	AddScalarBar(legend);

}
void VisualizeDisplay::clipVector(QString vtkFile, int normal, double origin, bool normalize) {
	
	reader = NULL;
 	reader = vtkSmartPointer<vtkStructuredPointsReader>::New();
  reader->SetFileName(vtkFile.toStdString().c_str());
  reader->Update();

	geometryFilter = NULL;
  geometryFilter =
    vtkSmartPointer<vtkGeometryFilter>::New();
  //geometryFilter->SetInputConnection(reader->GetOutputPort());

  vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();

  switch(normal) {
  case 0:
    plane->SetOrigin(origin,0,0);
    plane->SetNormal(1,0,0);
    break;
  case 1:
    plane->SetOrigin(0,origin,0);
    plane->SetNormal(0,1,0);
    break;
  case 2:
    plane->SetOrigin(0,0,origin);
    plane->SetNormal(0,0,1);
    break;
  }

  vtkSmartPointer<vtkVectorNorm> vecScalar = vtkSmartPointer<vtkVectorNorm>::New();
  vecScalar->SetInputConnection(reader->GetOutputPort());
  if(!normalize) vecScalar->NormalizeOff();
  else {
    vecScalar->NormalizeOn();
    vecScalar->SetNormalize(normal);
  }

  vtkSmartPointer<vtkCutter> cutter = vtkSmartPointer<vtkCutter>::New();
  cutter->SetInputConnection(vecScalar->GetOutputPort());
  cutter->SetCutFunction(plane);

  vtkSmartPointer<vtkContourFilter> iso = vtkSmartPointer<vtkContourFilter>::New();
  iso->SetInputConnection(vecScalar->GetOutputPort());
  iso->SetValue(0, 0.5);

  geometryFilter->SetInputConnection(cutter->GetOutputPort());

  vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
  /* lut->SetNumberOfColors(64);*/
  lut->SetHueRange(0.667, 0.0);
  lut->Build();

  RemoveActor(VTKactor);

	VTKmapper = NULL;
	VTKactor = NULL;
  VTKmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  VTKactor = vtkSmartPointer<vtkActor>::New();
 
  VTKmapper->SetInputConnection(geometryFilter->GetOutputPort());
  VTKmapper->SetLookupTable(lut);
  //VTKmapper->UseLookupTableScalarRangeOff();
  VTKmapper->SetScalarRange(geometryFilter->GetOutput()->GetScalarRange());

  VTKactor->SetMapper(VTKmapper);

  AddActor(VTKactor);

	RemoveScalarBar(legend);
	legend = NULL;
	legend = vtkSmartPointer<vtkScalarBarActor>::New();
	legend->SetLookupTable(VTKmapper->GetLookupTable());
	legend->SetWidth(0.1);
	legend->GetLabelTextProperty()->SetFontSize(9);
	AddScalarBar(legend);

}
void VisualizeDisplay::clear() {
  if(VTKactor != NULL)    RemoveActor(VTKactor);
}
void VisualizeDisplay::setRange(double a, double b) {
  if (VTKmapper != NULL) VTKmapper->SetScalarRange(a,b);
}
void VisualizeDisplay::getRange(double &a) {
  VTKmapper->GetScalarRange(&a);
}
void VisualizeDisplay::vector(QString vtkFile, int normal, double origin) {

	vectReader = NULL;
 	vectReader = vtkSmartPointer<vtkStructuredPointsReader>::New();
  vectReader->SetFileName(vtkFile.toStdString().c_str());
  vectReader->Update();

  vectGeometryFilter =
    vtkSmartPointer<vtkGeometryFilter>::New();
  //geometryFilter->SetInputConnection(reader->GetOutputPort());

  vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();

  switch(normal) {
  case 0:
    plane->SetOrigin(origin,0,0);
    plane->SetNormal(1,0,0);
    break;
  case 1:
    plane->SetOrigin(0,origin,0);
    plane->SetNormal(0,1,0);
    break;
  case 2:
    plane->SetOrigin(0,0,origin);
    plane->SetNormal(0,0,1);
    break;
  }

  vtkSmartPointer<vtkCutter> cutter = vtkSmartPointer<vtkCutter>::New();
  cutter->SetInputConnection(vectReader->GetOutputPort());
  cutter->SetCutFunction(plane);
  
  vectGeometryFilter->SetInputConnection(cutter->GetOutputPort());

  vtkSmartPointer<vtkArrowSource> arrow = vtkSmartPointer<vtkArrowSource>::New();
  arrow->SetTipResolution(6);
  arrow->SetTipRadius(0.1);
  arrow->SetTipLength(0.35);
  arrow->SetShaftResolution(6);
  arrow->SetShaftRadius(0.03);

  vtkSmartPointer<vtkGlyph3D> glyph = vtkSmartPointer<vtkGlyph3D>::New();
  glyph->SetInputConnection(cutter->GetOutputPort());
  glyph->SetSourceConnection(arrow->GetOutputPort());
  glyph->SetVectorModeToUseVector();
//  glyph->SetColorModeToColorByScalar();
  glyph->SetScaleModeToScaleByVector();
  glyph->SetOrient(normal);
  glyph->OrientOn();
  glyph->SetScaleFactor(0.3);


  RemoveActor(vectActor);

  vectMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  vectActor = vtkSmartPointer<vtkActor>::New();
 
  vectMapper->SetInput(glyph->GetOutput());
  //VTKmapper->SetLookupTable(lut);
  vectMapper->ScalarVisibilityOn();
  vectMapper->SetScalarRange(vectReader->GetOutput()->GetScalarRange());

  vectActor->SetMapper(vectMapper);

  AddActor(vectActor);

}

void VisualizeDisplay::hideVector() {

  RemoveActor(vectActor);
}
