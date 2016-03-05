/*
 * MainWindow_Boundaries.cpp
 *
 * Implements the Boundaries tab
 */

#include <QTreeWidget>
#include <QMessageBox>
#include <QDebug>

#include "MainWindow.h"

void MainWindow::buildBoundaryList() {
  QTreeWidgetItem *item = ui.BoundaryTree->topLevelItem(0);
  item->child(0)->setText(1,sim.getBoundaryText(0));
  item->child(1)->setText(1,sim.getBoundaryText(1));
  item->child(2)->setText(1,sim.getBoundaryText(2));
  item->child(3)->setText(1,sim.getBoundaryText(3));
  item->child(4)->setText(1,sim.getBoundaryText(4));
  item->child(5)->setText(1,sim.getBoundaryText(5));

  QTreeWidgetItem *new_item;
  long int extent_a[2];
  long int extent_b[2];
  QString type;
  QString extents;
  double value;
  double turbulence;
  int n;
  for(int i=0; i<6; i++) {
    item = ui.BoundaryTree->topLevelItem(1)->child(i);

    qDeleteAll(item->takeChildren());

    sim.resetSpecialBoundary(i);
    while(sim.getNextSpecialBoundary(type, extent_a, extent_b, value, turbulence)) {
      new_item = new QTreeWidgetItem(QTreeWidgetItem::Type);
      
      extents = QString::number(extent_a[0]) + 
                "," + QString::number(extent_a[1]) +
                " to " + QString::number(extent_b[0]) +
                "," + QString::number(extent_b[1]);
      n = item->childCount();

      item->insertChild(n, new_item);
      new_item->setText(0, QString::number(n+1));
      new_item->setText(1, type);
      new_item->setText(2, extents);
      new_item->setText(3, QString::number(value));
      new_item->setText(4, QString::number(turbulence));
    }
  }  


}

void MainWindow::boundariesUpdate() {

  buildBoundaryList();

  boundaryDisplay->update(sim.getDelx().toDouble(), sim.getDely().toDouble(),
                      sim.getDelz().toDouble(),
                      sim.getImax().toDouble(), sim.getJmax().toDouble(),
                      sim.getKmax().toDouble(),
                      sim.getOrigin(0).toDouble(), sim.getOrigin(1).toDouble(), sim.getOrigin(2).toDouble());
  
  boundaryDisplay->connectVTK("vtk/fv_0.vtk");
  ui.VTKBoundaries->update();
  
  boundaryDisplay->clearRectangle();
 
}

void MainWindow::editBoundary(QTreeWidgetItem *item) {
  int wb;

  if(item == NULL) return;

  if(item->childCount() > 0)
    return;
  
  // This code is for a wall boundary 
  if(ui.BoundaryTree->indexOfTopLevelItem(item->parent()) == 0) {

    boundaryDialog = new BoundaryDialog(item->text(1));
    boundaryDialog->exec();

    if(boundaryDialog->result() == QDialog::Accepted) {
      wb = item->parent()->indexOfChild(item);
      sim.setBoundaryText(wb, boundaryDialog->getBoundaryText());
    }

    delete boundaryDialog;
  }
  // This code is for a special boundary
  else {
    QTreeWidgetItem *item = ui.BoundaryTree->currentItem();
    QString text;
    QString type;
    long int extent_a[2];
    long int extent_b[2];
    double value, turbulence;

    if(item->parent() == NULL) text = "";
    else  text = item->parent()->text(0);
  
    if((text != "West"   && text !="East"   &&
      text != "South"  && text != "North" &&
      text != "Bottom" && text != "Top"  ) ||
      ui.BoundaryTree->indexOfTopLevelItem(item->parent()) == 0) {
        QMessageBox::information(this,"Civil CFD",
                             "Please select from the list of special boundaries", QMessageBox::Ok, QMessageBox::Ok);
                                 
        return;
    }

    int wall = item->parent()->parent()->indexOfChild(item->parent());
    sim.resetSpecialBoundary(wall);

    for(int i = 0; i < item->text(0).toInt(); i++) {
       sim.getNextSpecialBoundary(type, extent_a, extent_b, value, turbulence);
    }

    sBoundaryDialog = new SBoundaryDialog(type, wall, extent_a[0], extent_a[1], extent_b[0], extent_b[1], value, turbulence);
    sBoundaryDialog->exec();

    if(sBoundaryDialog->result() == QDialog::Accepted) {
      sim.resetSpecialBoundary(wall);
      for(int i = 1; i < item->text(0).toInt(); i++) {
        sim.nextSpecialBoundary();
      }     
      sim.editSpecialBoundary( sBoundaryDialog->getBoundaryText(),
      sBoundaryDialog->getExtentA1(), sBoundaryDialog->getExtentA2(),
      sBoundaryDialog->getExtentB1(), sBoundaryDialog->getExtentB2(),
      sBoundaryDialog->getValue(),    sBoundaryDialog->getTurbulence());
    }

    delete sBoundaryDialog;
  }

  //update();
  buildBoundaryList();
  boundaryDisplay->clearRectangle();
}

void MainWindow::on_EditBoundary_clicked() {
  QTreeWidgetItem *item = ui.BoundaryTree->currentItem();
  editBoundary(item);
}

void MainWindow::on_BoundaryTree_itemDoubleClicked(QTreeWidgetItem *item, int column) {
  editBoundary(item);
}

void MainWindow::on_AddSpecialBoundary_clicked() {
  int wall;
  QString text;
  QTreeWidgetItem *item = ui.BoundaryTree->currentItem();

  if(item == NULL) {
    QMessageBox::information(this,"Civil CFD",
                             "Please select a wall from the list of special boundaries", QMessageBox::Ok, QMessageBox::Ok);
    return;
  }

  text = item->text(0);
  if((text != "West"   && text !="East"   &&
     text != "South"  && text != "North" &&
     text != "Bottom" && text != "Top"  ) ||
     ui.BoundaryTree->indexOfTopLevelItem(item->parent()) == 0) {
    QMessageBox::information(this,"Civil CFD",
                             "Please select a wall from the list of special boundaries", QMessageBox::Ok, QMessageBox::Ok);
                                 
    return;
  }

  wall = item->parent()->indexOfChild(item);

  sBoundaryDialog = new SBoundaryDialog(wall);
  sBoundaryDialog->exec();

  if(sBoundaryDialog->result() == QDialog::Accepted) {
    wall = item->parent()->indexOfChild(item);
    sim.addSpecialBoundary(wall, 
      sBoundaryDialog->getBoundaryText(),
      sBoundaryDialog->getExtentA1(), sBoundaryDialog->getExtentA2(),
      sBoundaryDialog->getExtentB1(), sBoundaryDialog->getExtentB2(),
      sBoundaryDialog->getValue(),    sBoundaryDialog->getTurbulence());
  }

  delete sBoundaryDialog;

  
  //update();
  buildBoundaryList();
  boundaryDisplay->clearRectangle();
}

void MainWindow::on_RemoveSpecialBoundary_clicked() {
  QTreeWidgetItem *item = ui.BoundaryTree->currentItem();
  QString text;

  if(item->parent() == NULL) text = "";
  else  text = item->parent()->text(0);
  
  if((text != "West"   && text !="East"   &&
     text != "South"  && text != "North" &&
     text != "Bottom" && text != "Top"  ) ||
     ui.BoundaryTree->indexOfTopLevelItem(item->parent()) == 0) {
    QMessageBox::information(this,"Civil CFD",
                             "Please select from the list of special boundaries", QMessageBox::Ok, QMessageBox::Ok);
                                 
    return;
  }

  int wall = item->parent()->parent()->indexOfChild(item->parent());
  sim.resetSpecialBoundary(wall);

  for(int i = 1; i < item->text(0).toInt(); i++) {
    sim.nextSpecialBoundary();
  }

  sim.removeSpecialBoundary(wall);
 
  
  //update();
  buildBoundaryList();
  boundaryDisplay->clearRectangle();
}

void MainWindow::on_BoundaryTree_itemClicked(QTreeWidgetItem *item, int column) {
  on_BoundaryTree_itemActivated(item, column);
}

void MainWindow::on_BoundaryTree_itemActivated(QTreeWidgetItem *item, int column) {
  QString text;
  double a_1, a_2, a_3, b_1, b_2, b_3;
  int normal;

  if(item==NULL) {
    boundaryDisplay->clearRectangle();
    return;
  }

 // This code is for a wall boundary 
  if(ui.BoundaryTree->indexOfTopLevelItem(item->parent()) == 0) {
    text = item->text(0);

    if(text == "West") {
      a_1 = a_2 = a_3 = 0;

      b_1 = 0;
      b_2 = (sim.getJmax().toLong()-1) * sim.getDely().toDouble();
      b_3 = (sim.getKmax().toLong()-1) * sim.getDelz().toDouble();

      normal = 0;
    }
    else if(text == "East") {
      a_2 = a_3 = 0;

      a_1 = b_1 = (sim.getImax().toLong()-1) * sim.getDelx().toDouble();
      
      b_2 = (sim.getJmax().toLong()-1) * sim.getDely().toDouble();
      b_3 = (sim.getKmax().toLong()-1) * sim.getDelz().toDouble();

      normal = 0;
    }
    else if(text == "South") {
      a_1 = a_2 = a_3 = 0;

      b_2 = 0;
      b_1 = (sim.getImax().toLong()-1) * sim.getDelx().toDouble();
      b_3 = (sim.getKmax().toLong()-1) * sim.getDelz().toDouble();

      normal = 1;
    }
    else if(text == "North") {
      a_1 = a_3 = 0;

      a_2 = b_2 = (sim.getJmax().toLong()-1) * sim.getDely().toDouble();
      
      b_1 = (sim.getImax().toLong()-1) * sim.getDelx().toDouble();
      b_3 = (sim.getKmax().toLong()-1) * sim.getDelz().toDouble();

      normal = 1;
    }
    else if(text == "Bottom") {
      a_1 = a_2 = a_3 = 0;

      b_3 = 0;
      b_1 = (sim.getImax().toLong()-1) * sim.getDelx().toDouble();
      b_2 = (sim.getJmax().toLong()-1) * sim.getDely().toDouble();

      normal = 2;
    }
    else if(text == "Top") {
      a_1 = a_2 = 0;

      a_3 = b_3 = (sim.getKmax().toLong()-1) * sim.getDelz().toDouble();
      
      b_1 = (sim.getImax().toLong()-1) * sim.getDelx().toDouble();
      b_2 = (sim.getJmax().toLong()-1) * sim.getDely().toDouble();

      normal = 2;
    }
    else {
      boundaryDisplay->clearRectangle();
      ui.VTKBoundaries->update();
      return;
    }

		a_1 +=  + sim.getOrigin(0).toDouble();
		b_1 +=  + sim.getOrigin(0).toDouble();
		a_2 +=  + sim.getOrigin(1).toDouble();
		b_2 +=  + sim.getOrigin(1).toDouble();
		a_3 +=  + sim.getOrigin(2).toDouble();
		b_3 +=  + sim.getOrigin(2).toDouble();

    boundaryDisplay->drawRectangle(a_1, a_2, a_3, b_1, b_2, b_3, normal);
  }
  else {
    QString type;
    long int extent_a[2];
    long int extent_b[2];
    double value, turbulence;
    bool numeric;

    item->text(0).toInt(&numeric);
    if(!numeric) {
      boundaryDisplay->clearRectangle();
    }
    else {
      int wall = item->parent()->parent()->indexOfChild(item->parent());
      sim.resetSpecialBoundary(wall);

      for(int i = 0; i < item->text(0).toInt(); i++) {
        sim.getNextSpecialBoundary(type, extent_a, extent_b, value, turbulence);
      }

      text = item->parent()->text(0);

      if(text == "West") {
        a_1 = 0;
        a_2 = extent_a[0] * sim.getDely().toDouble();
        a_3 = extent_a[1] * sim.getDelz().toDouble();

        b_1 = 0;
        b_2 = (extent_b[0]-1) * sim.getDely().toDouble();
        b_3 = (extent_b[1]-1) * sim.getDelz().toDouble();

        normal = 0;
      }
      else if(text == "East") {
        a_1 = (sim.getImax().toLong() - 1) * sim.getDelx().toDouble();
        a_2 = extent_a[0] * sim.getDely().toDouble();
        a_3 = extent_a[1] * sim.getDelz().toDouble();

        b_1 = a_1;
        b_2 = (extent_b[0]-1) * sim.getDely().toDouble();
        b_3 = (extent_b[1]-1) * sim.getDelz().toDouble();

        normal = 0;
      }
      else if(text == "South") {
        a_2 = 0;
        a_1 = extent_a[0] * sim.getDelx().toDouble();
        a_3 = extent_a[1] * sim.getDelz().toDouble();

        b_2 = 0;
        b_1 = (extent_b[0]-1) * sim.getDelx().toDouble();
        b_3 = (extent_b[1]-1) * sim.getDelz().toDouble();

        normal = 1;
      }
      else if(text == "North") {
        a_2 = (sim.getJmax().toLong() - 1) * sim.getDely().toDouble();
        a_1 = extent_a[0] * sim.getDelx().toDouble();
        a_3 = extent_a[1] * sim.getDelz().toDouble();

        b_2 = a_2;
        b_1 = (extent_b[0]-1) * sim.getDelx().toDouble();
        b_3 = (extent_b[1]-1) * sim.getDelz().toDouble();

        normal = 1;
      }
      else if(text == "Bottom") {
        a_3 = 0;
        a_1 = extent_a[0] * sim.getDelx().toDouble();
        a_2 = extent_a[1] * sim.getDely().toDouble();

        b_3 = 0;
        b_1 = (extent_b[0]-1) * sim.getDelx().toDouble();
        b_2 = (extent_b[1]-1) * sim.getDely().toDouble();
        
        normal = 2;

      }
      else if(text == "Top") {
        a_3 = (sim.getKmax().toLong() - 1) * sim.getDelz().toDouble();
        a_1 = extent_a[0] * sim.getDelx().toDouble();
        a_2 = extent_a[1] * sim.getDely().toDouble();

        b_3 = a_3;
        b_1 = (extent_b[0]-1) * sim.getDelx().toDouble();
        b_2 = (extent_b[1]-1) * sim.getDely().toDouble();

        normal = 2;
      }
      else return;


			a_1 +=  + sim.getOrigin(0).toDouble();
			b_1 +=  + sim.getOrigin(0).toDouble();
			a_2 +=  + sim.getOrigin(1).toDouble();
			b_2 +=  + sim.getOrigin(1).toDouble();
			a_3 +=  + sim.getOrigin(2).toDouble();
			b_3 +=  + sim.getOrigin(2).toDouble();   

      boundaryDisplay->drawRectangle(a_1, a_2, a_3, b_1, b_2, b_3, normal);
    }
  }
  
  ui.VTKBoundaries->update();
}
