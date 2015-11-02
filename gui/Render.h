/*
 * Render.h 
 *
 * User interface for render dialog
 */

#ifndef RENDER_H
#define RENDER_H

#include <QDialog>
#include <QProcess>

#include "ui_render.h"

#include "RenderDisplay.h"
#include "Simulation.h"

class RenderDialog : public QDialog {

  Q_OBJECT

public:

  RenderDialog(Simulation &sim, QString appPath);

private slots:
  void error(QProcess::ProcessError error);
  void finished(int exitCode, QProcess::ExitStatus status);
  void stateChanged(QProcess::ProcessState state);
  void readyReadStandardOutput();
  void on_Return_clicked();
  void on_Stop_clicked();

private:
  Ui::Render ui;

  RenderDisplay *renderDisplay;

  QProcess *process;

  bool stopped;

};



#endif
