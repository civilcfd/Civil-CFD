/*
 * Render.cpp
 *
 * Implementation of render dialog
 */

#include "Render.h"

RenderDialog::RenderDialog(Simulation &sim, QString appPath) {
  QString cmd;

  stopped = false;

  ui.setupUi(this);

  ui.status->setText("Starting...");
  ui.Return->setEnabled(false);

  process = new QProcess(this);

  connect(process,
            SIGNAL(error(QProcess::ProcessError)),
            SLOT(error(QProcess::ProcessError)));
  connect(process,
            SIGNAL(finished(int, QProcess::ExitStatus)),
            SLOT(finished(int, QProcess::ExitStatus)));
  connect(process,
            SIGNAL(stateChanged(QProcess::ProcessState)),
            SLOT(stateChanged(QProcess::ProcessState)));
  connect(process,
            SIGNAL(readyReadStandardOutput()),
            SLOT(readyReadStandardOutput()));
			
#ifdef _WIN32
  cmd = appPath + "/mesh3d.exe";
#else
  cmd = appPath + "/mesh3d";
#endif

  if(!QFile::exists(cmd)) {
    ui.status->setText("Could not find executable: " + cmd);
    ui.Return->setEnabled(true);
  }
  else {
    cmd = "\"" + cmd + "\"";

    cmd = cmd + " solver.xml " + sim.getStlFilename(); 

    process->start(cmd);
  }

  renderDisplay = new RenderDisplay(sim.getImax().toLong(), 
                                    sim.getJmax().toLong(),
                                    sim.getKmax().toLong(),
                                    sim.getDelx().toDouble(),
                                    sim.getDely().toDouble(),
                                    sim.getDelz().toDouble(),
                                    sim.getOrigin(0).toDouble(),
                                    sim.getOrigin(1).toDouble(),
                                    sim.getOrigin(2).toDouble());

  ui.VTKRender->SetRenderWindow(renderDisplay->getRenderWindow());
  renderDisplay->completeSetup( sim.getDelx().toDouble(),
                                    sim.getDely().toDouble(),
                                    sim.getDelz().toDouble(),
                                    
                                    sim.getImax().toLong(), 
                                    sim.getJmax().toLong(),
                                    sim.getKmax().toLong(),
                                   
                                    sim.getOrigin(0).toDouble(),
                                    sim.getOrigin(1).toDouble(),
                                    sim.getOrigin(2).toDouble());
}

void RenderDialog::on_Return_clicked() {
  close();
}

void RenderDialog::on_Stop_clicked() {

  stopped = true;
  
#ifdef _WIN32
  process->kill();
#else
  process->terminate();
#endif

}

void RenderDialog::readyReadStandardOutput() {
  QString str = process->readAllStandardOutput(); 
  ui.output->appendPlainText(str);

  if(str.contains("Marking cells")) {
    ui.progressBar->setValue(5);
    ui.status->setText("Marking cells with no intersections");
  }
  if(str.contains("area fraction")) {
    ui.progressBar->setValue(30);
    ui.status->setText("Calculating area fractions");
  }
  if(str.contains("volume fractions")) {
    ui.progressBar->setValue(60);
    ui.status->setText("Calculating volume fractions");
  }
  if(str.contains("around obstacles")) {
    ui.progressBar->setValue(85);
    ui.status->setText("Filling mesh cells around obstacles");
  }
  if(str.contains("small mesh cells")) {
    ui.progressBar->setValue(90);
    ui.status->setText("Eliminating very small mesh cells");
  }
  if(str.contains("mesh to file")) {
    ui.progressBar->setValue(95);
    ui.status->setText("Writing mesh to file");
  }
}

void RenderDialog::error(QProcess::ProcessError error) {

  if(stopped == true) {
    ui.status->setText("Stopped by user");
  } else if(error == QProcess::FailedToStart) {
    ui.status->setText("Error: mesh3d could not start");
  } else if(error == QProcess::Crashed) {
    ui.status->setText("Error: mesh3d has crashed");
  }

  ui.Return->setEnabled(true);
}

void RenderDialog::finished(int exitCode, QProcess::ExitStatus status) {
  QString str = process->readAllStandardOutput(); 
  ui.output->appendPlainText(str);
  
  ui.progressBar->setValue(100);

  if(exitCode == 0) {
    if(stopped == false) ui.status->setText("Finished successfully");
    else ui.status->setText("Stopped by user");

    renderDisplay->connectVTK("vtk/fv_0.vti");
    ui.VTKRender->update();
  }
  
  if(exitCode == 1) {
    ui.status->setText("Process failed");
  }

  ui.Return->setEnabled(true);
}

void RenderDialog::stateChanged(QProcess::ProcessState state) {
}
