/*
 * Boundary.cpp
 *
 * Implements the boundary select dialog
 *
 */

#include <QDebug>
#include <QMessageBox>
#include <QInputDialog>

#include "Boundary.h"

#include <math.h>


BoundaryDialog::BoundaryDialog(QString text) {

  ui.setupUi(this);

  ui.select->setCurrentIndex(ui.select->findText(text, Qt::MatchFixedString));

}

void BoundaryDialog::accept() {

  setBoundaryText(ui.select->currentText());
  QDialog::accept();

}

void BoundaryDialog::setBoundaryText(QString text) {
  boundaryText = text;
}

QString BoundaryDialog::getBoundaryText() {
  return boundaryText;
}

SBoundaryDialog::SBoundaryDialog(QString type, int wall, long int extentA1, long int extentA2, long int extentB1, long int extentB2, double value, double turbulence) {
  ui.setupUi(this);

  ui.extentA1->setPlainText(QString::number(extentA1));
  ui.extentA2->setPlainText(QString::number(extentA2));
  ui.extentB1->setPlainText(QString::number(extentB1));
  ui.extentB2->setPlainText(QString::number(extentB2));
 
  ui.value->setPlainText(QString::number(value));
  ui.turbulence->setPlainText(QString::number(turbulence));

	ui.select->setCurrentIndex(ui.select->findText(type));

  setLabels(wall);
}

SBoundaryDialog::SBoundaryDialog(int wall) {
  ui.setupUi(this);

  setLabels(wall);
}

void SBoundaryDialog::setLabels(int wall) {

  switch(wall) {
  case 0: // west
    ui.A1->setText("min Y");
    ui.B1->setText("max Y");
    ui.A2->setText("min Z");
    ui.B2->setText("max Z");
    break;
  case 1: // east
    ui.A1->setText("min Y");
    ui.B1->setText("max Y");
    ui.A2->setText("min Z");
    ui.B2->setText("max Z");
    break;
  case 2: // south
    ui.A1->setText("min X");
    ui.B1->setText("max X");
    ui.A2->setText("min Z");
    ui.B2->setText("max Z");
    break;
  case 3: // north
    ui.A1->setText("min X");
    ui.B1->setText("max X");
    ui.A2->setText("min Z");
    ui.B2->setText("max Z");
    break;
  case 4: // bottom
    ui.A1->setText("min X");
    ui.B1->setText("max X");
    ui.A2->setText("min Y");
    ui.B2->setText("max Y");
    break;
  case 5: // top
    ui.A1->setText("min X");
    ui.B1->setText("max X");
    ui.B1->setText("min Y");
    ui.B2->setText("max Y");
    break;
  }
}

void SBoundaryDialog::accept() {
  bool ok_value, ok_turbulence, ok_extent_a[2], ok_extent_b[2];

  boundaryText = ui.select->currentText();
  value = ui.value->toPlainText().toDouble(&ok_value);
  turbulence = ui.turbulence->toPlainText().toDouble(&ok_turbulence);
  extent_a[0] = ui.extentA1->toPlainText().toLong(&ok_extent_a[0]);
  extent_a[1] = ui.extentA2->toPlainText().toLong(&ok_extent_a[1]);
  extent_b[0] = ui.extentB1->toPlainText().toLong(&ok_extent_b[0]);
  extent_b[1] = ui.extentB2->toPlainText().toLong(&ok_extent_b[1]);

  if(!ok_value || !ok_turbulence || !ok_extent_a[0] ||
     !ok_extent_a[1] || !ok_extent_b[0] || !ok_extent_b[1]) {
    QMessageBox::warning(this, "Civil CFD", "Type checking error: Boundary extents must be integers.  All fields must be numeric.", QMessageBox::Ok, QMessageBox::Ok);
    return;
  }

  QDialog::accept();
}

void SBoundaryDialog::on_calculate_clicked() {

	double velocity;
	bool ok_value;
	
	if(ui.select->currentText() == "fixed velocity") {
		velocity = ui.value->toPlainText().toDouble(&ok_value);
		if(!ok_value) return;
	}
	else {
		bool ok; 

		QString text =  QInputDialog::getText(this ,"Civil CFD",
																						"Freestream velocity:", QLineEdit::Normal,
																						"0", &ok);
					
		if (ok && !text.isEmpty()) {
			velocity = text.toDouble(&ok);
		}
		else return;
  }

	turbulence = 1.5 * pow(0.02,2) * pow(velocity, 2);
	ui.turbulence->setPlainText(QString::number(turbulence));

}

double SBoundaryDialog::getValue() {
  return value;
}

double SBoundaryDialog::getTurbulence() {
  return turbulence;
}

long int SBoundaryDialog::getExtentA1() {
  return extent_a[0];
}

long int SBoundaryDialog::getExtentA2() {
  return extent_a[1];
}

long int SBoundaryDialog::getExtentB1() {
  return extent_b[0];
}

long int SBoundaryDialog::getExtentB2() {
  return extent_b[1];
}

QString SBoundaryDialog::getBoundaryText() {
  return boundaryText;
}

void SBoundaryDialog::on_select_currentIndexChanged(QString str) {

  if(str == "fixed velocity") {
    ui.units->setText("m/s");
  }
  else if (str == "mass outflow" ){
    ui.units->setText("m3/s");
  }
  else {
    ui.units->setText("m");
  }
  
  if(str == "virtual weir") {
    ui.labelValue->setText("Weir length");
    ui.labelTurbulence->setText("Weir invert");
    ui.calculate->setEnabled(false);
  } else {
    ui.labelValue->setText("Value");
    ui.labelTurbulence->setText("Turbulence");
    ui.calculate->setEnabled(true);
  }
  

}