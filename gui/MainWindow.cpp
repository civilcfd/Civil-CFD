/*
 * MainWindow.cpp: implements routines
 */

#include <QDebug>
#include <QFileDialog>
#include <QMessageBox>
#include <QInputDialog>

#include <cmath>

#include "MainWindow.h"

MainWindow::MainWindow() : QMainWindow() {
  meshDisplay = new MeshDisplay(10,10,10);
  geometryDisplay = new GeometryDisplay(10,10,10);
  boundaryDisplay = new BoundaryDisplay(10,10,10);
  baffleDisplay = new BoundaryDisplay(10,10,10);
  visualizeDisplay = new VisualizeDisplay(10,10,10);

  appPath = QCoreApplication::applicationDirPath();

  ui.setupUi( this );

  ui.VTKMesh->SetRenderWindow(meshDisplay->getRenderWindow());
  //ui.VTKMesh->GetRenderWindow()->SetInteractor(meshDisplay->getRenderWindowInteractor());
  //ui.VTKMesh->GetRenderWindow()->GetInteractor()->Start();

  ui.VTKGeometry->SetRenderWindow(geometryDisplay->getRenderWindow());

  ui.VTKBoundaries->SetRenderWindow(boundaryDisplay->getRenderWindow());  
  
  ui.VTKBaffles->SetRenderWindow(baffleDisplay->getRenderWindow());  

  ui.vis->SetRenderWindow(visualizeDisplay->getRenderWindow());
}

bool MainWindow::saveNotify() {
  QString name;
  
  if(QMessageBox::question(this, "Save required", "This will save your simulation.  Proceed?",
       QMessageBox::Ok|QMessageBox::Cancel) == QMessageBox::Cancel)
      return false;

  name = ui.Name->toPlainText();
  if(name.isEmpty()) {
    QMessageBox msgBox;
    msgBox.setText("Please enter a project name");
    msgBox.exec();    
    return false;
  }

  qDebug() << "saving:" << name;

  sim.setName(name);
  write();
  sim.save();
  
  return true;
}

void MainWindow::on_RunSimulation_clicked() {
  if(!saveNotify()) return;

  solverDialog = new SolverDialog(sim, appPath, ui.t->currentText());
  solverDialog->exec();

  delete solverDialog;

  sim.loadSimulation(ui.Path->toPlainText());
  update();
}

void MainWindow::on_STLRender_clicked() {
  if(!saveNotify()) return;

  renderDialog = new RenderDialog(sim, appPath);
  renderDialog->exec();

  delete renderDialog;

  update();
}

void MainWindow::on_STLOpen_clicked() {
   QString file = QFileDialog::getOpenFileName(this, tr("Open STL File"), QDir::homePath(), tr("STL Files (*.stl)"));

  if(file.isEmpty()) {
    qDebug() << "no file chosen";
    return;
  }

  qDebug() << "opening file:" << file;

   if(!sim.loadSTL(file)) {
    QMessageBox msgBox;
    msgBox.setText("Failed to open STL file");
    msgBox.exec();    
    return;
  }

  ui.STLFile->setPlainText(file);
  ui.STLFile_2->setPlainText(file);

  geometryDisplay->connectSTL(file);
  ui.VTKGeometry->update();

  return;
}

void MainWindow::on_MeshUndo_clicked() {
  QTreeWidgetItem *item = ui.MeshParameters->topLevelItem(0);
  item->child(0)->setText(1,sim.getDelx());
  item->child(1)->setText(1,sim.getDely());
  item->child(2)->setText(1,sim.getDelz());
  item = ui.MeshParameters->topLevelItem(1);
  item->child(0)->setText(1,sim.getImax());
  item->child(1)->setText(1,sim.getJmax());
  item->child(2)->setText(1,sim.getKmax());
  item = ui.MeshParameters->topLevelItem(2);
  item->child(0)->setText(1,sim.getOrigin(0));
  item->child(1)->setText(1,sim.getOrigin(1));
  item->child(2)->setText(1,sim.getOrigin(2));

}

void MainWindow::on_MeshUpdate_clicked() {
  meshUpdate();
  update();
}

void MainWindow::on_kEpsilon_toggled() {
  toggle();
}

void MainWindow::on_Laminar_toggled() {
  toggle();
}


void MainWindow::on_actionSave_triggered() {
  on_Save_clicked();
}

void MainWindow::on_Save_clicked() {
  QString name;

  name = ui.Name->toPlainText();
  if(name.isEmpty()) {
    QMessageBox msgBox;
    msgBox.setText("Please enter a project name");
    msgBox.exec();    
    return;
  }

  qDebug() << "saving:" << name;

  sim.setName(name);
  write();
  sim.save();
  
  return;
}
void MainWindow::on_actionNew_triggered() {
  on_New_clicked();
}
void MainWindow::on_New_clicked() {
  QString path;
  QString path_file;
  QString name;

  name = ui.Name->toPlainText();
  if(name.isEmpty()) {
    QMessageBox msgBox;
    msgBox.setText("Please enter a project name");
    msgBox.exec();    
    return;
  }

  path = QFileDialog::getExistingDirectory(this, tr("Choose Project Directory"), QDir::homePath());

  if(path.isEmpty()) {
    qDebug() << "no path chosen";
    return;
  }

  path_file = path + "/casefile";
  QFile ftest(path_file);
  if(ftest.exists()) {
    if(QMessageBox::question(this, "Case exists", "Overwrite existing case?",
       QMessageBox::Yes|QMessageBox::No) == QMessageBox::No)
      return;
  }

  qDebug() << "new project:" << name;
  qDebug() << "path:" << path;

  if(!sim.newSimulation(name,path)) {
    QMessageBox msgBox;
    msgBox.setText("Failed to create simulation");
    msgBox.exec();    
    return;
  }

  ui.Path->setPlainText(path);

  enableAll();
 
  update();

  return;
}

void MainWindow::on_actionQuit_triggered() {
  QApplication::quit();
}

void MainWindow::on_actionOpen_triggered() {
   QString path = QFileDialog::getExistingDirectory(this, tr("Choose Project Directory"), QDir::homePath());

  if(path.isEmpty()) {
    qDebug() << "no path chosen";
    return;
  }

  QString path_file = path + "/casefile";
  QFile ftest(path_file);
  if(!ftest.exists()) {
    QMessageBox::warning(this, "Case does not exist", "No case found in this directory",  QMessageBox::Ok);
    return;
  }

  qDebug() << "opening path:" << path;

   if(!sim.loadSimulation(path)) {
    QMessageBox msgBox;
    msgBox.setText("Failed to create simulation");
    msgBox.exec();    
    return;
  }

  ui.Path->setPlainText(path);
  enableAll();

  update();

  return;
}

void MainWindow::enableAll() {
  setEnabled(true);
  ui.Save->setEnabled(true);
  ui.Models->setEnabled(true);
  ui.Mesh->setEnabled(true);
  ui.MeshParameters->expandAll();
  ui.Geometry->setEnabled(true);  
  ui.Boundaries->setEnabled(true);
  ui.Baffles->setEnabled(true);
  ui.BoundaryTree->expandAll();
  ui.BaffleTree->expandAll();
  ui.InitialConditions->setEnabled(true);
  ui.Simulate->setEnabled(true);
  ui.Visualize->setEnabled(true);
}

void MainWindow::disableAll() {
  setEnabled(false);
}

void MainWindow::on_MeshParameters_itemDoubleClicked(QTreeWidgetItem *item, int column) {
 if(column == 1) {
   ui.MeshParameters->editItem(item, column);
 }

}

void MainWindow::meshUpdate() {
  /* write mesh values to the solver
   * and then update the QVTKWidget */

  QTreeWidgetItem *item = ui.MeshParameters->topLevelItem(0);
  sim.setDelx(item->child(0)->text(1));
  sim.setDely(item->child(1)->text(1));
  sim.setDelz(item->child(2)->text(1));
  item = ui.MeshParameters->topLevelItem(1);
  sim.setImax(item->child(0)->text(1));
  sim.setJmax(item->child(1)->text(1));
  sim.setKmax(item->child(2)->text(1));
  item = ui.MeshParameters->topLevelItem(2);
  sim.setOrigin(0,item->child(0)->text(1));
  sim.setOrigin(1,item->child(1)->text(1));
  sim.setOrigin(2,item->child(2)->text(1));
 
  meshDisplay->update(sim.getDelx().toDouble(), sim.getDely().toDouble(),
                      sim.getDelz().toDouble(),
                      sim.getImax().toDouble(), sim.getJmax().toDouble(),
                      sim.getKmax().toDouble(),
                      sim.getOrigin(0).toDouble(), sim.getOrigin(1).toDouble(),
                      sim.getOrigin(2).toDouble());
  geometryDisplay->update(sim.getDelx().toDouble(), sim.getDely().toDouble(),
                      sim.getDelz().toDouble(),
                      sim.getImax().toDouble(), sim.getJmax().toDouble(),
                      sim.getKmax().toDouble(),
                      sim.getOrigin(0).toDouble(), sim.getOrigin(1).toDouble(),
                      sim.getOrigin(2).toDouble());

  movePoint();
  ui.VTKMesh->update();
  geometryDisplay->connectSTL(sim.getStlFilename());
  ui.VTKGeometry->update();
}
 
void MainWindow::update() {
  /* update all constants displayed in the application */
  ui.Name->setPlainText(sim.getName());
  ui.Path->setPlainText(sim.getPath());
 
  QString Dims;
  Dims = QString::number(sim.getImax().toLong() * sim.getDelx().toDouble()) + " x " + QString::number(sim.getJmax().toLong() * sim.getDely().toDouble()) + " x " + QString::number(sim.getKmax().toLong() * sim.getDelz().toDouble()) + " m";
  ui.Dims->setPlainText(Dims);

  QString ElementSize;
  ElementSize = sim.getDelx() + " x " + sim.getDely() + " x " + sim.getDelz() + " m";
  ui.ElementSize->setPlainText(ElementSize);

  QString ElementNumber;
  ElementNumber = sim.getImax() + " x " + sim.getJmax() + " x " + sim.getKmax() + " cells";
  ui.ElementNumber->setPlainText(ElementNumber);
 
  QString OutputSize;
  OutputSize = QString::number((sim.getImax().toLong() * sim.getJmax().toLong() * sim.getKmax().toLong()) * 0.125) + " kB per saved timestep";
  ui.OutputSize->setPlainText(OutputSize);
  
  ui.ResultList->clear();
  ui.t->clear();
  buildResultList();
  

  /* Models Page */
  ui.gx->setPlainText(sim.getGx());
  ui.gy->setPlainText(sim.getGy());
  ui.gz->setPlainText(sim.getGz());
  
  ui.nu->setPlainText(sim.getNu());
  ui.rho->setPlainText(sim.getRho());

  if(sim.laminar()) {
    ui.Laminar->setChecked(true);
    ui.kEpsilon->setChecked(false);
  }
  else if(sim.kEpsilon()) {
    ui.Laminar->setChecked(false);
    ui.kEpsilon->setChecked(true);
    ui.rough->setPlainText(sim.getRough());
    ui.length_scale->setPlainText(sim.getLength_scale());
  }
  toggle();

  QTreeWidgetItem *item = ui.MeshParameters->topLevelItem(0);
  item->child(0)->setText(1,sim.getDelx());
  item->child(1)->setText(1,sim.getDely());
  item->child(2)->setText(1,sim.getDelz());
  item = ui.MeshParameters->topLevelItem(1);
  item->child(0)->setText(1,sim.getImax());
  item->child(1)->setText(1,sim.getJmax());
  item->child(2)->setText(1,sim.getKmax());
  item = ui.MeshParameters->topLevelItem(2);
  item->child(0)->setText(1,sim.getOrigin(0));
  item->child(1)->setText(1,sim.getOrigin(1));
  item->child(2)->setText(1,sim.getOrigin(2));

  meshUpdate();

  // Geometry tab

  ui.inside_x->setPlainText(sim.getInside(0));
  ui.inside_y->setPlainText(sim.getInside(1));
  ui.inside_z->setPlainText(sim.getInside(2));

  ui.STLFile->setPlainText(sim.getStlFilename());
  ui.STLFile_2->setPlainText(sim.getStlFilename());

  boundariesUpdate();
  bafflesUpdate();

  // Initial Values Tab
  double vector[3];
  sim.getInitialVector("velocity", vector);  
  
  ui.initialU->setPlainText(QString::number(vector[0]));
  ui.initialV->setPlainText(QString::number(vector[1]));
  ui.initialW->setPlainText(QString::number(vector[2]));
  ui.height->setPlainText(sim.getInitialScalar("vof_height"));
  ui.initialk->setPlainText(sim.getInitialScalar("kE_k"));

  if(sim.getInitialScalar("hydrostatic").toDouble()>0) ui.hydrostatic->setChecked(true);
  else ui.hydrostatic->setChecked(false);

  ui.fillPoints->clearContents();
  int fillPointIndex = 0;
  while(sim.getInitialVector("inside", vector) != 1); // cycle to beginning of list
  while(sim.getInitialVector("inside", vector) != 1) { // now populate fillPoints table 
    ui.fillPoints->setItem(fillPointIndex, 0, new QTableWidgetItem(QString::number(vector[0])));
    ui.fillPoints->setItem(fillPointIndex, 1, new QTableWidgetItem(QString::number(vector[1])));
    ui.fillPoints->setItem(fillPointIndex, 2, new QTableWidgetItem(QString::number(vector[2])));
    fillPointIndex++;
  }
   
  // Simulate Tab

  ui.t->setCurrentIndex(ui.t->findText(sim.getT(), Qt::MatchContains));
  ui.endt->setPlainText(sim.getEndt());
  ui.writet->setPlainText(sim.getWritet());
  ui.delt->setPlainText(sim.getDelt());
  ui.autot->setChecked(sim.getAutot());


  visualizeUpdate();
}

void MainWindow::write() {
  /* write all constants into the solver */
  sim.setGx(ui.gx->toPlainText());
  sim.setGy(ui.gy->toPlainText());
  sim.setGz(ui.gz->toPlainText());

  sim.setNu(ui.nu->toPlainText());
  sim.setRho(ui.rho->toPlainText());
  if(ui.kEpsilon->isChecked()) {
    sim.setTurbulence("kEpsilon");
    sim.setRough(ui.rough->toPlainText());
    sim.setLength_scale(ui.length_scale->toPlainText());
  }
  else {
    sim.setTurbulence("Laminar");
  }
  toggle();  
  meshUpdate();  

  sim.setInside(0, ui.inside_x->toPlainText());
  sim.setInside(1, ui.inside_y->toPlainText());
  sim.setInside(2, ui.inside_z->toPlainText());

  // Initial Values Tab

  sim.setInitial("velocity", ui.initialU->toPlainText(),
                             ui.initialV->toPlainText(),
                             ui.initialW->toPlainText());

  sim.setInitial("vof_height", ui.height->toPlainText());
  sim.setInitial("kE_k", ui.initialk->toPlainText());

  if(ui.hydrostatic->isChecked()) sim.setInitial("hydrostatic", "1");
  else sim.setInitial("hydrostatic", "0");

  QString fp1, fp2, fp3;
  for(int fillPointIndex = 0; fillPointIndex < 8; fillPointIndex++) {
    
      if(!ui.fillPoints->item(fillPointIndex, 0) ||
         !ui.fillPoints->item(fillPointIndex, 1) ||
         !ui.fillPoints->item(fillPointIndex, 2)) continue;
    
      fp1 = ui.fillPoints->item(fillPointIndex, 0)->text();
      fp2 = ui.fillPoints->item(fillPointIndex, 1)->text();
      fp3 = ui.fillPoints->item(fillPointIndex, 2)->text();
    
      if(!fp1.isEmpty() && !fp2.isEmpty() && !fp3.isEmpty())    
        sim.setInitial("inside", fp1, fp2, fp3);
  }
  
  // This must be called to flush out initial values
  sim.setInitial("end");

  // Simulate tab

  sim.setT(ui.t->currentText());
  sim.setEndt(ui.endt->toPlainText());
  sim.setWritet(ui.writet->toPlainText());
  sim.setDelt(ui.delt->toPlainText());
  sim.setAutot(ui.autot->isChecked());

}

void MainWindow::toggle() {
  if(ui.kEpsilon->isChecked()) {
    ui.length_scale->setEnabled(true);
    ui.initialk->setEnabled(true);
    ui.rough->setEnabled(true);
    
    if(ui.length_scale->toPlainText().isEmpty() && 
       ui.rough->toPlainText().isEmpty()) {
      sim.setTurbulence("kEpsilon");
      ui.rough->setPlainText(sim.getRough());
      ui.length_scale->setPlainText(sim.getLength_scale());
    }
    
    if(ui.initialk->toPlainText() == "-1") ui.initialk->setPlainText("0.001");
  }
  else {
    ui.length_scale->setEnabled(false);
    ui.rough->setEnabled(false);
    ui.initialk->setEnabled(false);
  }
}

void MainWindow::on_inside_x_textChanged() {
   sim.setInside(0, ui.inside_x->toPlainText()); 
   movePoint();
}

void MainWindow::on_inside_y_textChanged() {
 sim.setInside(1, ui.inside_y->toPlainText()); 
 movePoint();
}

void MainWindow::on_inside_z_textChanged() {
 sim.setInside(2, ui.inside_z->toPlainText()); 
 movePoint();
}

void MainWindow::movePoint() {
  geometryDisplay->clearPoint();
  geometryDisplay->drawPoint(sim.getInside(0).toDouble(),
                             sim.getInside(1).toDouble(),
                             sim.getInside(2).toDouble());
  ui.VTKGeometry->update();
}

void MainWindow::on_earthGravity_clicked() {
  ui.gx->setPlainText("0");
  ui.gy->setPlainText("0");
  ui.gz->setPlainText("-9.81");
}

void MainWindow::on_water20C_clicked() {
  ui.nu->setPlainText("1.004e-06");
  ui.rho->setPlainText("998.2");
}

void MainWindow::on_defaultLength_clicked() {
  ui.length_scale->setPlainText("0.038");
}

void MainWindow::on_calcRough_clicked() {

  bool ok; 

  QString text =  QInputDialog::getText(this ,"Civil CFD",
                                          "Manning's n:", QLineEdit::Normal,
                                          "0.013", &ok);
          
  if (ok && !text.isEmpty()) {
    double ks, n;
    n = text.toDouble(&ok);

    if(ok) {
      ks = n * 8.41 * std::sqrt(std::abs(sim.getGz().toDouble()));
      ks = std::pow(ks, 6);
      ui.rough->setPlainText(QString::number(ks,'g',3));
    }
  }
}


