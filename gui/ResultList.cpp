/*
 * ResultList.cpp: maintain the result list
 */

#include <QDebug>
#include <QFileDialog>
#include <QMessageBox>
#include <QFile>
#include <QProcess>
#include <QSettings>

#include "MainWindow.h"

void MainWindow::buildResultList() {
  sim.trackRewind();

  ui.t->clear();
  ui.t->addItem("0");
  
  bool ok;
  QString t_last;
  
  while(sim.getTrackNext() >= 0) {
    ui.ResultList->addItem(sim.getTrackT());
    if(sim.getTrackT().toDouble(&ok) > 0.00000001)
      if(ok) t_last = sim.getTrackT();
  }

  if(t_last.toDouble(&ok) > 0.0000001)
  {
    if(ok)  {
      ui.t->addItem(t_last);
      ui.t->setCurrentIndex(1);
    } else {
      ui.t->setCurrentIndex(0);
    }
  }
  

}

void MainWindow::on_SelectAll_clicked() {
	ui.ResultList->selectAll();
}

void MainWindow::on_Clear_clicked() {
	ui.ResultList->clearSelection();
}


void MainWindow::on_Delete_clicked() {
	
	int index;
	int flg=0;
	
	foreach(QListWidgetItem *item, ui.ResultList->selectedItems()) {
		if(sim.deleteTrack(item->text())) {
			index = ui.t->findText(item->text());
			if(index != -1) {
				flg = 1;
			}
			
			delete item;
		}
	}
	
  buildTimesteps();
  build3dTimesteps();
  
  if(flg == 1) {
    sim.trackRewind();

    ui.t->clear();
    ui.t->addItem("0");
  
    bool ok;
    QString t_last;
  
    while(sim.getTrackNext() >= 0) {
      if(sim.getTrackT().toDouble(&ok) > 0.00000001)
        if(ok) t_last = sim.getTrackT();
    }

    if(t_last.toDouble(&ok) > 0.0000001)
    {
      if(ok)  {
        ui.t->addItem(t_last);
        ui.t->setCurrentIndex(1);
      } else {
        ui.t->setCurrentIndex(0);
      }
    } 
  }
}

void MainWindow::on_Paraview_clicked() {
	
	int flg=0;
  QString index;
  QString timestep;
  QProcess *process;

  if(ui.ResultList->selectedItems().count() < 1) {
    QMessageBox msgBox;
    msgBox.setText("Please select some timesteps to open");
    msgBox.exec();    
    return;
  }

  QFile file("fv.pvd");
	if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {
    QTextStream out(&file);
    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" >\n ";
    out << "<Collection>\n";
    foreach(QListWidgetItem *item, ui.ResultList->selectedItems()) {
      index = sim.getTrackN(item->text());
      timestep = item->text();
      out << "<DataSet timestep=\"" << timestep << "\" file=\"vtk/fv_0.vti\" />\n";
    }
    out << "</Collection>\n";
    out << "</VTKFile>\n";
    file.close();
  }

  file.setFileName("U.pvd");
	if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {
    QTextStream out(&file);
    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" >\n ";
    out << "<Collection>\n";
    foreach(QListWidgetItem *item, ui.ResultList->selectedItems()) {
      index = sim.getTrackN(item->text());
      timestep = item->text();
      out << "<DataSet timestep=\"" << timestep << "\" file=\"vtk/U_" << index << ".vti\" />\n";
    }
    out << "</Collection>\n";
    out << "</VTKFile>\n";
    file.close();
  }

  file.setFileName("P.pvd");
	if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {
    QTextStream out(&file);
    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" >\n ";
    out << "<Collection>\n";
    foreach(QListWidgetItem *item, ui.ResultList->selectedItems()) {
      index = sim.getTrackN(item->text());
      timestep = item->text();
      out << "<DataSet timestep=\"" << timestep << "\" file=\"vtk/P_" << index << ".vti\" />\n";
    }
    out << "</Collection>\n";
    out << "</VTKFile>\n";
    file.close();
  }
 
  file.setFileName("vof.pvd");
	if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {
    QTextStream out(&file);
    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" >\n ";
    out << "<Collection>\n";
    foreach(QListWidgetItem *item, ui.ResultList->selectedItems()) {
      index = sim.getTrackN(item->text());
      timestep = item->text();
      out << "<DataSet timestep=\"" << timestep << "\" file=\"vtk/vof_" << index << ".vti\" />\n";
    }
    out << "</Collection>\n";
    out << "</VTKFile>\n";
    file.close();
  }

  file.setFileName("vorticity.pvd");
	if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {
    QTextStream out(&file);
    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" >\n ";
    out << "<Collection>\n";
    foreach(QListWidgetItem *item, ui.ResultList->selectedItems()) {
      index = sim.getTrackN(item->text());
      timestep = item->text();
      out << "<DataSet timestep=\"" << timestep << "\" file=\"vtk/vorticity_" << index << ".vti\" />\n";
    }
    out << "</Collection>\n";
    out << "</VTKFile>\n";
    file.close();
  }

  if(sim.kEpsilon()) {

    file.setFileName("k.pvd");
    if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {
      QTextStream out(&file);
      out << "<?xml version=\"1.0\"?>\n";
      out << "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" >\n ";
      out << "<Collection>\n";
      foreach(QListWidgetItem *item, ui.ResultList->selectedItems()) {
        index = sim.getTrackN(item->text());
        timestep = item->text();
        out << "<DataSet timestep=\"" << timestep << "\" file=\"vtk/k_" << index << ".vti\" />\n";
      }
      out << "</Collection>\n";
      out << "</VTKFile>\n";
      file.close();
    }

    file.setFileName("E.pvd");
    if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {
      QTextStream out(&file);
      out << "<?xml version=\"1.0\"?>\n";
      out << "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" >\n ";
      out << "<Collection>\n";
      foreach(QListWidgetItem *item, ui.ResultList->selectedItems()) {
        index = sim.getTrackN(item->text());
        timestep = item->text();
        out << "<DataSet timestep=\"" << timestep << "\" file=\"vtk/E_" << index << ".vti\" />\n";
      }
      out << "</Collection>\n";
      out << "</VTKFile>\n";
      file.close();
    }

  }

  file.setFileName("results.py");
  if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {
    QTextStream out(&file);
    out << "from paraview.simple import *\n";
    out << "paraview.simple._DisableFirstRenderCameraReset() \n";

    out << "fvpvd = PVDReader(FileName='fv.pvd') \n";
    out << "RenameSource('fv', fvpvd) \n";
    out << "ppvd = PVDReader(FileName='P.pvd') \n";
    out << "RenameSource('P', ppvd) \n";
    out << "Upvd = PVDReader(FileName='U.pvd') \n";
    out << "RenameSource('U', Upvd) \n";
    out << "vofpvd = PVDReader(FileName='vof.pvd') \n";
    out << "RenameSource('vof', vofpvd) \n";
    out << "vorticitypvd = PVDReader(FileName='vorticity.pvd') \n";
    out << "RenameSource('vorticity', vorticitypvd) \n";

    if(sim.kEpsilon()) {
      out << "kpvd = PVDReader(FileName='k.pvd') \n";
      out << "RenameSource('k', kpvd) \n";
      out << "Epvd = PVDReader(FileName='E.pvd') \n";
      out << "RenameSource('E', Epvd) \n";
    }
    out << "";
    file.close();
  }

  QString cmd;

  #ifdef _WIN32

  QSettings settings(appPath + "/Paraview.ini", QSettings::IniFormat);
  cmd = settings.value("paraview_path", "paraview.exe").toString();

  if(!QFile::exists(cmd)) {
    cmd = QFileDialog::getOpenFileName(this, tr("Where is the Paraview Executable"), QDir::homePath(), tr("EXE Files (*.exe)"));

    if(cmd.isEmpty()) {
      qDebug() << "no file chosen";
      return;
    }
    settings.setValue("paraview_path", cmd);

  }

  cmd = "\"" + cmd + "\"";
  #else
  cmd = "paraview";
  #endif

  cmd = cmd + " --script=results.py";

  process = new QProcess(this);
  process->start(cmd);
}
