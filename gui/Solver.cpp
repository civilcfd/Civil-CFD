/*
 * Solver.cpp
 *
 * Implementation of solver dialog
 */

#include "SolverDialog.h"

SolverDialog::SolverDialog(Simulation &sim, QString appPath, QString t) {
  QString cmd;
  bool ok;

  stopT = sim.getEndt().toDouble();

  ui.setupUi(this);

  ui.status->setText("Simulation Running");

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

  cmd = appPath + "/solver3d";

  if(t.toDouble(&ok) < 0.00001) t = "";
  if(!ok) t="";
  
  if(t != "") cmd = cmd + " " + t;
  
  if(!QFile::exists(cmd)) {
    ui.status->setText("Could not find executable: " + cmd);
    ui.Return->setEnabled(true);
  }
  else {
    if(t != "") cmd = cmd + " " + t;
    process->start(cmd);
  }

}

void SolverDialog::on_Return_clicked() {
  close(); 
}

void SolverDialog::on_Stop_clicked() {
  process->terminate();

  stopped = true;
}

void SolverDialog::readyReadStandardOutput() {
  QString str = process->readAllStandardOutput(); 
  ui.output->appendPlainText(str);

  if(str.contains("timestep")) {
    int progressVal = 0;

    QStringList list = str.split(" ");
    for(int i=0; i < list.count()-1; i++) {
      if(list[i].contains("timestep")) {
        progressVal = list[i+1].toDouble() * 100;
        break;
      }
    }

    if(progressVal > 0) {
      progressVal = progressVal / stopT;
      ui.progressBar->setValue(progressVal);
    }
  }

}

void SolverDialog::error(QProcess::ProcessError error) {

  if(error == QProcess::FailedToStart) {
    ui.status->setText("Error: solver3d could not start");
  } else if(error == QProcess::Crashed) {
    ui.status->setText("Error: solver3d has crashed");
  }

  ui.Return->setEnabled(true);
}

void SolverDialog::finished(int exitCode, QProcess::ExitStatus status) {
  ui.progressBar->setValue(100);

  if(exitCode == 0) {
    if(stopped == false) ui.status->setText("Finished successfully");
    else ui.status->setText("Stopped by user");
  }
  
  if(exitCode == 1) {
    ui.status->setText("Process failed");
  }

  ui.Return->setEnabled(true);
}

void SolverDialog::stateChanged(QProcess::ProcessState state) {
}
