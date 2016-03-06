/*
 * Boundary.h 
 *
 * User interface for boundary select dialogs
 */

#ifndef BOUNDARY_H 
#define BOUNDARY_H

#include <QDialog>

#include "ui_boundary.h"
#include "ui_sboundary.h"

#include "Simulation.h"

class BoundaryDialog : public QDialog {

  Q_OBJECT

public:

  BoundaryDialog(QString text);

  QString getBoundaryText();

private slots:
  void accept();

private:
  QString boundaryText;
  void setBoundaryText(QString text);

  Ui::BoundaryDialog ui;

};

class SBoundaryDialog : public QDialog {

  Q_OBJECT

public:

  SBoundaryDialog(int wall);
  SBoundaryDialog(QString type, int wall, long int extentA1, long int extentA2, long int extentB1, long int extentB2, double value, double turbulence); 
 
  QString getBoundaryText();
  double getValue();
  double getTurbulence();
  long int getExtentA1();
  long int getExtentA2();
  long int getExtentB1();
  long int getExtentB2();

private slots:
  void accept();
  void on_calculate_clicked();
  void on_select_currentIndexChanged(QString str);

private:
  void setLabels(int wall);

  Ui::SBoundaryDialog ui;
  double value;
  double turbulence;

  QString boundaryText;
  long int extent_a[2];
  long int extent_b[2]; 
};

#endif
