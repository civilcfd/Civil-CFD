/*
 * Visualize3D.cpp: implements the 3d isosurface tab
 * To maximize code re-use, this mostly rehashes Visualize 
 * The code to block obstacles is very similar to that needed to create the fluid isosurface
 */

#include <QDebug>
#include <QMessageBox>
#include <QPixmap>
#include <QFileDialog>

#include <vtkSmartPointer.h>
#include <vtkStructuredPoints.h>
#include <vtkStructuredPointsReader.h>
#include <vtkCubeSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkPlane.h>
#include <vtkCutter.h>
#include <vtkProperty.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkOBJReader.h>
#include <vtkLookupTable.h>

#include "MainWindow.h"

void MainWindow::visualize3dUpdate() {

  if(build3dTimesteps() < 1) {
    visualize3dDisplay->clear();
    return;
  }
  
  on_origin3d_valueChanged();
    
  visualize3dRender();    
  
  
  visualize3dDisplay->update(sim.getDelx().toDouble(), sim.getDely().toDouble(),
                      sim.getDelz().toDouble(),
                      sim.getImax().toDouble(), sim.getJmax().toDouble(),
                      sim.getKmax().toDouble(),
                      sim.getOrigin(0).toDouble(), sim.getOrigin(1).toDouble(), sim.getOrigin(2).toDouble());

  //visualizeDisplay->connectSTL(sim.getStlFilename());
 

}

void MainWindow::on_timesteps3d_currentItemChanged() {
  visualize3dRender();
}
void MainWindow::on_xNormal3d_toggled() {
  updateSlider3d();
  visualize3dRender();
}
void MainWindow::on_yNormal3d_toggled() {
  updateSlider3d();
  visualize3dRender();
}
void MainWindow::on_zNormal3d_toggled() {
  updateSlider3d();
  visualize3dRender();
}
/*
void MainWindow::on_contourVOF_toggled() {
  visualizeRender();
}
void MainWindow::on_contourP_toggled() {
  visualizeRender();
}
void MainWindow::on_contourK_toggled() {
  visualizeRender();
}
void MainWindow::on_showAxis_toggled() {
  visualizeRender();
}
void MainWindow::on_showLegend_toggled() {
  visualizeRender();
}
void MainWindow::on_contourVorticity_toggled() {
  visualizeRender();
} 
void MainWindow::on_showVectors_toggled() {
  visualizeRender();
} */
void MainWindow::on_blockObstacles3d_toggled() {
  if(!ui.blockObstacles3d->isChecked())
    visualize3dDisplay->hideBlock();

  visualize3dRender();
}
void MainWindow::on_showMesh3d_toggled() {


  visualize3dRender();
}
void MainWindow::on_showAxis3d_toggled() {
  visualize3dRender();
}

void MainWindow::on_extent3d_valueChanged() {
	updateSlider3d();
	
	visualize3dRender();
}

void MainWindow::updateSlider3d() {
  
  double extent, maxExtent, position, extentPos;
  bool ok;

  if(ui.xNormal3d->isChecked()) {
    maxExtent = sim.getDelx().toDouble(&ok) * sim.getImax().toDouble(&ok);
  } else if(ui.yNormal3d->isChecked()) {
    maxExtent = sim.getDely().toDouble(&ok) * sim.getJmax().toDouble(&ok); 
  } else {
    maxExtent = sim.getDelz().toDouble(&ok) * sim.getKmax().toDouble(&ok);
  }

  position = (double) ui.origin3d->sliderPosition();
  position /= 99;
  position *= maxExtent;
  
  extentPos = (double) ui.extent3d->sliderPosition();
  extent = (extentPos / 99) * (maxExtent - position);

  ui.origin3dText->setText(QString::number(position, 'f', 2));
  
  ui.extent3dText->setText(QString::number(extent, 'f', 2));
  
}

void MainWindow::on_origin3d_valueChanged() {
  updateSlider3d();

  visualize3dRender();
}

int MainWindow::build3dTimesteps() {
  int count = 0;

  ui.timesteps3d->clear();  
  
  sim.trackRewind();

  while(sim.getTrackNext() >= 0) {
    ui.timesteps3d->addItem(sim.getTrackT());
    count++;
  }
  
  return count;
}
/*
void MainWindow::on_updateRange_clicked() {
  double a, b;
  bool ok_a, ok_b;
  
  a = ui.from->toPlainText().toDouble(&ok_a);
  b = ui.to->toPlainText().toDouble(&ok_b);
  
  if(b>a && ok_a && ok_b) {
    visualizeDisplay->setRange(a,b);
    ui.vis->update();
  }
}
*/
void MainWindow::on_saveJPEG3d_clicked() {
  QString file = QFileDialog::getSaveFileName(this, tr("Save Image"), 
    QDir::homePath(), tr("png files (*.png)"));
  QPixmap p = QPixmap::grabWidget(ui.vis3d);
  p.save(QString(file),0,100);
  
}

void MainWindow::visualize3dRender() {

    double origin, del;
    int normal;

    if(ui.timesteps3d->currentItem() == NULL) return;
    QListWidgetItem *item = ui.timesteps3d->currentItem();

    QString vtkFile;
    vtkFile = sim.getTrackN(item->text());
    
    
    QString volFile;
    bool vtkFile_d = false;
    
    vtkFile.append(".vtk");

  	vtkFile.prepend("vof_");

    volFile = "vtk/fv_0.vtk";
    vtkFile.prepend("vtk/");

  if(!QFile::exists(vtkFile)) {
    sim.decompressFile(vtkFile);
    if(!QFile::exists(vtkFile)) {
      QMessageBox msgBox;
      vtkFile.prepend("Failed to open VTK file: ");
      msgBox.setText(vtkFile);
      msgBox.exec();
    }
    else {
      vtkFile_d = true;
    }
  }

  origin = ui.origin3dText->text().toDouble();

  if(ui.xNormal3d->isChecked()) {
    normal = 0;
  } else if(ui.yNormal3d->isChecked()) {
    normal = 1;
  } else {
    normal = 2;
  }
  del = ui.extent3dText->text().toDouble();

	visualize3dDisplay->contour3d(vtkFile, normal, origin, del);


  if(ui.blockObstacles3d->isChecked())
    visualize3dDisplay->block(volFile, normal, origin, del);
    
  if(ui.showMesh3d->isChecked())
    visualize3dDisplay->ShowMesh();
  else
    visualize3dDisplay->HideMesh();
    
  if(ui.showAxis3d->isChecked()) 
    visualize3dDisplay->ShowAxis();
  else
    visualize3dDisplay->HideAxis();
    /*
  if(ui.showLegend->isChecked()) 
    visualizeDisplay->showLegend();
  else
    visualizeDisplay->hideLegend();*/
  
  ui.vis3d->update();
  
  if(vtkFile_d) QFile::remove(vtkFile);
}
