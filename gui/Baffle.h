/*
 * Baffle.h 
 *
 * User interface for baffle select dialogs
 */

#ifndef BAFFLE_H 
#define BAFFLE_H

#include <QDialog>

#include "ui_baffle.h"

#include "Simulation.h"

class BaffleDialog : public QDialog {

  Q_OBJECT

public:

  BaffleDialog(int wall,
				QString maxA, QString maxB, QString maxPos);
  BaffleDialog(QString type, int wall, long int extentA1, long int extentA2, long int extentB1, long int extentB2, double value, long int pos,
				QString maxA, QString maxB, QString maxPos); 
 
  QString getBaffleText();
  double getValue();
  long int getPos();
  long int getExtentA1();
  long int getExtentA2();
  long int getExtentB1();
  long int getExtentB2();

private slots:
  void accept();
  void on_select_currentIndexChanged(QString str);

private:
  void setLabels(int wall);

  Ui::BaffleDialog ui;
  double value;
  long int pos;

  QString baffleText;
  long int extent_a[2];
  long int extent_b[2]; 
};

#endif
