/*
 * Inspect.h 
 *
 * User interface for render dialog
 */

#ifndef INSPECT_H
#define INSPECT_H

#include <QDialog>
#include <QProcess>

#include "ui_inspect.h"

#include "Simulation.h"

class InspectDialog : public QDialog {

  Q_OBJECT

public:

  InspectDialog(Simulation &sim, QString path, QString t);

private slots:
  void error(QProcess::ProcessError error);
  void finished(int exitCode, QProcess::ExitStatus status);
  void stateChanged(QProcess::ProcessState state);
  void readyReadStandardOutput();
  void on_Return_clicked();
  void on_InspectCell_clicked();

private:
  Ui::Inspect ui;
  
  QString timestep;
  QString appPath;

  QProcess *process;


};



#endif
