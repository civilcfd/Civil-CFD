/*
 * About.cpp: implementation
 */
 
#include "About.h"

AboutDialog::AboutDialog() {
  ui.setupUi(this);
}

void AboutDialog::on_OK_clicked() {
	
  QDialog::accept();
}