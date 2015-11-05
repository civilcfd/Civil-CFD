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

void MainWindow::on_SelectAll_clicked() {
	ui.ResultList->selectAll();
}

void MainWindow::on_Clear_clicked() {
	ui.ResultList->clearSelection();
}


void MainWindow::on_Delete_clicked() {
	
	foreach(QListWidgetItem *item, ui.ResultList->selectedItems()) {
		if(sim.deleteTrack(item->text()))
			delete item;
	}
	
  buildTimesteps();
}
