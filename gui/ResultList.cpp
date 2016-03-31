/*
 * ResultList.cpp: maintain the result list
 */

#include <QDebug>
#include <QFileDialog>
#include <QMessageBox>

#include "MainWindow.h"

void MainWindow::buildResultList() {
  sim.trackRewind();

  ui.t->clear();
  ui.t->addItem("0");
  
  bool ok;
  QString t_last;
  
  while(sim.getTrackNext() >= 0) {
    ui.ResultList->addItem(sim.getTrackT());
    if(sim.getTrackT().toDouble(&ok) > 0.00000001)
      if(ok) t_last = sim.getTrackT();
  }

  if(t_last.toDouble(&ok) > 0.0000001)
    if(ok) ui.t->addItem(t_last);

}

void MainWindow::on_SelectAll_clicked() {
	ui.ResultList->selectAll();
}

void MainWindow::on_Clear_clicked() {
	ui.ResultList->clearSelection();
}


void MainWindow::on_Delete_clicked() {
	
	int index;
	int flg=0;
	
	foreach(QListWidgetItem *item, ui.ResultList->selectedItems()) {
		if(sim.deleteTrack(item->text())) {
			index = ui.t->findText(item->text());
			if(index != -1) {
				flg = 1;
			}
			
			delete item;
		}
	}
	
  buildTimesteps();
  build3dTimesteps();
  
  if(flg == 1) {
    sim.trackRewind();

    ui.t->clear();
    ui.t->addItem("0");
  
    bool ok;
    QString t_last;
  
    while(sim.getTrackNext() >= 0) {
      if(sim.getTrackT().toDouble(&ok) > 0.00000001)
        if(ok) t_last = sim.getTrackT();
    }

    if(t_last.toDouble(&ok) > 0.0000001)
      if(ok) ui.t->addItem(t_last);    
  }
}
