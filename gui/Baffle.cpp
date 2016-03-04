/*
 * Baffle.cpp
 *
 * Implements the baffle select dialog
 *
 */

#include <QDebug>
#include <QMessageBox>
#include <QInputDialog>

#include "Baffle.h"

#include <math.h>

BaffleDialog::BaffleDialog(QString type, int wall, long int extentA1, long int extentA2, long int extentB1, long int extentB2, double value, long int pos) {
  ui.setupUi(this);

  ui.extentA1->setPlainText(QString::number(extentA1));
  ui.extentA2->setPlainText(QString::number(extentA2));
  ui.extentB1->setPlainText(QString::number(extentB1));
  ui.extentB2->setPlainText(QString::number(extentB2));
 
  ui.value->setPlainText(QString::number(value));
  ui.pos->setPlainText(QString::number(pos));

	ui.select->setCurrentIndex(ui.select->findText(type));

  setLabels(wall);
}

BaffleDialog::BaffleDialog(int wall) {
  ui.setupUi(this);

  setLabels(wall);
}

void BaffleDialog::setLabels(int wall) {

  switch(wall) {
  case 0: // x
    ui.A1->setText("min Y");
    ui.B1->setText("max Y");
    ui.A2->setText("min Z");
    ui.B2->setText("max Z");
    break;
  case 1: // y
    ui.A1->setText("min X");
    ui.B1->setText("max X");
    ui.A2->setText("min Z");
    ui.B2->setText("max Z");
    break;
  case 2: // z
    ui.A1->setText("min X");
    ui.B1->setText("max X");
    ui.A2->setText("min Y");
    ui.B2->setText("max Y");
    break;
  }
}

void BaffleDialog::accept() {
  bool ok_value, ok_pos, ok_extent_a[2], ok_extent_b[2];

  baffleText = ui.select->currentText();
  value = ui.value->toPlainText().toDouble(&ok_value);
  pos = ui.turbulence->toPlainText().toLong(&ok_pos);
  extent_a[0] = ui.extentA1->toPlainText().toLong(&ok_extent_a[0]);
  extent_a[1] = ui.extentA2->toPlainText().toLong(&ok_extent_a[1]);
  extent_b[0] = ui.extentB1->toPlainText().toLong(&ok_extent_b[0]);
  extent_b[1] = ui.extentB2->toPlainText().toLong(&ok_extent_b[1]);

  if(!ok_value || !ok_pos || !ok_extent_a[0] ||
     !ok_extent_a[1] || !ok_extent_b[0] || !ok_extent_b[1]) {
    QMessageBox::warning(this, "Civil CFD", "Type checking error: baffle extents must be integers.  All fields must be numeric.", QMessageBox::Ok, QMessageBox::Ok);
    return;
  }

  QDialog::accept();
}


double BaffleDialog::getValue() {
  return value;
}

double BaffleDialog::getPos() {
  return pos;
}

long int BaffleDialog::getExtentA1() {
  return extent_a[0];
}

long int BaffleDialog::getExtentA2() {
  return extent_a[1];
}

long int BaffleDialog::getExtentB1() {
  return extent_b[0];
}

long int BaffleDialog::getExtentB2() {
  return extent_b[1];
}

QString BaffleDialog::getBaffleText() {
  return baffleText;
}
