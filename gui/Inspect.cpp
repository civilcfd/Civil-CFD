/*
 * Inspect.cpp
 *
 * Implementation of inspect dialog
 */

#include <QMessageBox>

#include "Inspect.h"

InspectDialog::InspectDialog(Simulation &sim, QString path, QString t) {
  QString cmd;
  
  appPath = path;
  timestep = t;
  
  ui.maxX->text() = "of " + sim.getImax();
  ui.maxY->text() = "of " + sim.getJmax();
  ui.maxZ->text() = "of " + sim.getKmax();
  
  ui.setupUi(this);
}

void InspectDialog::on_Return_clicked() {
  close();
}

void InspectDialog::on_InspectCell_clicked() {  

  process = new QProcess(this);
  bool ok_i, ok_j, ok_k;
  QString i,j,k, cmd;

  i = ui.x->toPlainText();
  j = ui.y->toPlainText();
  k = ui.z->toPlainText();
  
  i.toLong(&ok_i);
  j.toLong(&ok_j);
  k.toLong(&ok_k);

  if(!ok_i || !ok_j || !ok_k) {
    QMessageBox msgBox;
    msgBox.setText("Please enter integer value for x, y, z");
    msgBox.exec();
    return;
  }

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
  cmd = appPath + "/inspect_cell.exe";
#else
  cmd = appPath + "/inspect_cell";
#endif


  if(!QFile::exists(cmd)) {
    ui.status->setText("Could not find executable: " + cmd);
  }
  else {
    cmd = cmd + " " + timestep + " " + i + " " + j + " " + k;

    process->start(cmd);
  }
  
  ui.Return->setEnabled(false);

}

void InspectDialog::readyReadStandardOutput() {
  QString str = process->readAllStandardOutput(); 
  ui.output->appendPlainText(str);
}

void InspectDialog::error(QProcess::ProcessError error) {

  ui.Return->setEnabled(true);
}

void InspectDialog::finished(int exitCode, QProcess::ExitStatus status) {

  ui.Return->setEnabled(true);
}

void InspectDialog::stateChanged(QProcess::ProcessState state) {
}
