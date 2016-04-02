/*
 * About.h 
 *
 * User interface for about dialog
 */

#ifndef ABOUT_H 
#define ABOUT_H

#include <QDialog>

#include "ui_about.h"

#include "Simulation.h"

class AboutDialog : public QDialog {

  Q_OBJECT

public:

  AboutDialog();

private slots:
  void on_OK_clicked();

private:
  Ui::AboutDialog ui;

};

#endif