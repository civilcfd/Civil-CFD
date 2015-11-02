/*
 * ResultList.cpp: maintain the result list
 */

#include <QDebug>
#include <QFileDialog>
#include <QMessageBox>

#include "MainWindow.h"

void MainWindow::buildResultList() {
  sim.trackRewind();

  while(sim.getTrackNext() >= 0) {
    ui.ResultList->addItem(sim.getTrackT());
  }

}

