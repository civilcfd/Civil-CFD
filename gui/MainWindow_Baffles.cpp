/*
 * MainWindow_Baffles.cpp
 *
 * Implements the Baffles tab
 */

#include <QTreeWidget>
#include <QMessageBox>
#include <QDebug>

#include "MainWindow.h"

void MainWindow::buildBaffleList() {
  QTreeWidgetItem *item, *new_item;
  long int extent_a[2];
  long int extent_b[2];
  QString type;
  QString extents;
  double value;
  long int pos;
  int n;
  
  for(int i=0; i<3; i++) {
    item = ui.BaffleTree->topLevelItem(0)->child(i);

    qDeleteAll(item->takeChildren());

    sim.resetBaffle(i);
    while(sim.getNextBaffle(type, extent_a, extent_b, value, pos)) {
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
      new_item->setText(3, QString::number(pos));
      new_item->setText(4, QString::number(value));
    }
  }  
}

void MainWindow::on_hideBaffleMesh_toggled() {
  if(ui.hideBaffleMesh->isChecked()) baffleDisplay->HideMesh();
  else baffleDisplay->ShowMesh();
	
	ui.VTKBaffles->update();
}

void MainWindow::on_hideBaffleDomain_toggled() {
  if(!ui.hideBaffleDomain->isChecked())  baffleDisplay->connectVTK("vtk/fv_0.vtk");
	else baffleDisplay->connectVTK("");
	
	ui.VTKBaffles->update();
}

void MainWindow::bafflesUpdate() { 
  
  buildBaffleList();

  baffleDisplay->update(sim.getDelx().toDouble(), sim.getDely().toDouble(),
                      sim.getDelz().toDouble(),
                      sim.getImax().toDouble(), sim.getJmax().toDouble(),
                      sim.getKmax().toDouble(),
                      sim.getOrigin(0).toDouble(), sim.getOrigin(1).toDouble(), sim.getOrigin(2).toDouble());
  
  ui.VTKBaffles->update();
  
  if(ui.hideBaffleMesh->isChecked()) baffleDisplay->HideMesh();
  else baffleDisplay->ShowMesh();
    
  if(!ui.hideBaffleDomain->isChecked())  baffleDisplay->connectVTK("vtk/fv_0.vtk");
	else baffleDisplay->connectVTK("");

  ui.VTKBaffles->update();
  baffleDisplay->clearRectangle();
 
}

void MainWindow::editBaffle(QTreeWidgetItem *item) {

  if(item == NULL) return;

  if(item->childCount() > 0)
    return;
  
  QString text;
  QString type;
  long int extent_a[2];
  long int extent_b[2];
  double value;
  long int pos;
  QString maxA, maxB, maxPos;

  if(item->parent() == NULL) text = "";
  else  text = item->parent()->text(0);

  if(text != "x"   && text !="y" && text != "z" ) {
      QMessageBox::information(this,"Civil CFD",
                           "Please select from the list of baffles", QMessageBox::Ok, QMessageBox::Ok);
                               
      return;
    }
    

  int wall = item->parent()->parent()->indexOfChild(item->parent());
  sim.resetBaffle(wall);

  for(int i = 0; i < item->text(0).toInt(); i++) {
     sim.getNextBaffle(type, extent_a, extent_b, value, pos);
  }
  
  switch(wall) {
  	case 0:
  		maxA = sim.getJmax();
  		maxB = sim.getKmax();
  		maxPos = sim.getImax();
  		break;
  	case 1:
  		maxA = sim.getImax();
  		maxB = sim.getKmax();
  		maxPos = sim.getJmax();
  		break;
  	case 2:
  		maxA = sim.getImax();
  		maxB = sim.getJmax();
  		maxPos = sim.getKmax();
  		break;
  }

  baffleDialog = new BaffleDialog(type, wall, extent_a[0], extent_a[1], extent_b[0], extent_b[1], value, pos, maxA, maxB, maxPos);
  baffleDialog->exec();

  if(baffleDialog->result() == QDialog::Accepted) {
    sim.resetBaffle(wall);
    for(int i = 1; i < item->text(0).toInt(); i++) {
      sim.nextBaffle();
    }     
    sim.editBaffle( baffleDialog->getBaffleText(),
    baffleDialog->getExtentA1(), baffleDialog->getExtentA2(),
    baffleDialog->getExtentB1(), baffleDialog->getExtentB2(),
    baffleDialog->getValue(),    baffleDialog->getPos());
  }

  delete baffleDialog;
  

  //update();
  buildBaffleList();
  baffleDisplay->clearRectangle();
}

void MainWindow::on_EditBaffle_clicked() {
  QTreeWidgetItem *item = ui.BaffleTree->currentItem();
  editBaffle(item);
}

void MainWindow::on_BaffleTree_itemDoubleClicked(QTreeWidgetItem *item, int column) {
	if(item == NULL) return;
	
	if(item->parent() == NULL) return;
	
	if(item->parent()->parent() == NULL) { // Double clicked on axis - shortcut to add item
		on_AddBaffle_clicked();
		return;
	}
	

  editBaffle(item);
}

void MainWindow::on_AddBaffle_clicked() {
  int wall;
  QString text;
  QTreeWidgetItem *item = ui.BaffleTree->currentItem();
  long int pos;
  QString maxA, maxB, maxPos;

  if(item == NULL) {
    QMessageBox::information(this,"Civil CFD",
                             "Please select an axis", QMessageBox::Ok, QMessageBox::Ok);
    return;
  }
  if(item->parent() == NULL) {
    QMessageBox::information(this,"Civil CFD",
                             "Please select an axis", QMessageBox::Ok, QMessageBox::Ok);
    return;
  }

	if(item->parent()->parent() != NULL) item = item->parent(); //Allows user to select from an existing baffle when adding

  text = item->text(0);
  if((text != "x"   && text !="y"   &&
      text != "z"  )) {
    QMessageBox::information(this,"Civil CFD",
                             "Please select an axis", QMessageBox::Ok, QMessageBox::Ok);
                                 
    return;
  }

  wall = item->parent()->indexOfChild(item);

  switch(wall) {
  	case 0:
  		maxA = sim.getJmax();
  		maxB = sim.getKmax();
  		maxPos = sim.getImax();
  		break;
  	case 1:
  		maxA = sim.getImax();
  		maxB = sim.getKmax();
  		maxPos = sim.getJmax();
  		break;
  	case 2:
  		maxA = sim.getImax();
  		maxB = sim.getJmax();
  		maxPos = sim.getKmax();
  		break;
  }
  
  baffleDialog = new BaffleDialog(wall, maxA, maxB, maxPos);
  baffleDialog->exec();

  if(baffleDialog->result() == QDialog::Accepted) {
    wall = item->parent()->indexOfChild(item);
    sim.addBaffle(wall, 
      baffleDialog->getBaffleText(),
      baffleDialog->getExtentA1(), baffleDialog->getExtentA2(),
      baffleDialog->getExtentB1(), baffleDialog->getExtentB2(),
      baffleDialog->getValue(),    baffleDialog->getPos());
  }

  delete baffleDialog;

  //update();
  buildBaffleList();
  baffleDisplay->clearRectangle();
}

void MainWindow::on_RemoveBaffle_clicked() {
  QTreeWidgetItem *item = ui.BaffleTree->currentItem();
  QString text;

  if(item->parent() == NULL) text = "";
  else  text = item->parent()->text(0);
  
  if(text != "x"   && text !="y" && text != "z" ) {
    QMessageBox::information(this,"Civil CFD",
                             "Please select from the list of special baffles", QMessageBox::Ok, QMessageBox::Ok);
                                 
    return;
  }

  int wall = item->parent()->parent()->indexOfChild(item->parent());
  sim.resetBaffle(wall);

  for(int i = 1; i < item->text(0).toInt(); i++) {
    sim.nextBaffle();
  }

  sim.removeBaffle(wall);
 
  //update();
  buildBaffleList();
  baffleDisplay->clearRectangle();
}

void MainWindow::on_BaffleTree_itemClicked(QTreeWidgetItem *item, int column) {
  on_BaffleTree_itemActivated(item, column);
}

void MainWindow::on_BaffleTree_itemActivated(QTreeWidgetItem *item, int column) {
  QString text;
  double a_1, a_2, a_3, b_1, b_2, b_3;
  int normal;

    QString type;
    long int extent_a[2];
    long int extent_b[2];
    double value;
    long int pos;
    bool numeric;

  if(item==NULL) {
    baffleDisplay->clearRectangle();
    return;
  }

    item->text(0).toInt(&numeric);
    if(!numeric) {
      baffleDisplay->clearRectangle();
    }
    else {
      int wall = item->parent()->parent()->indexOfChild(item->parent());
      sim.resetBaffle(wall);

      for(int i = 0; i < item->text(0).toInt(); i++) {
        sim.getNextBaffle(type, extent_a, extent_b, value, pos);
      }

    text = item->parent()->text(0);

    if(text == "x") {
      a_1 = (pos ) * sim.getDelx().toDouble();;
      a_2 = extent_a[0] * sim.getDely().toDouble();
      a_3 = extent_a[1] * sim.getDelz().toDouble();

      b_1 = (pos ) * sim.getDelx().toDouble();;
      b_2 = (extent_b[0]-1) * sim.getDely().toDouble();
      b_3 = (extent_b[1]-1) * sim.getDelz().toDouble();

      normal = 0;
    }
    else if(text == "y") {
      a_2 = (pos )  * sim.getDely().toDouble();
      a_1 = extent_a[0] * sim.getDelx().toDouble();
      a_3 = extent_a[1] * sim.getDelz().toDouble();

      b_2 = (pos )  * sim.getDely().toDouble();
      b_1 = (extent_b[0]-1) * sim.getDelx().toDouble();
      b_3 = (extent_b[1]-1) * sim.getDelz().toDouble();

      normal = 1;
    }
    else if(text == "z") {
      a_3 = (pos )  * sim.getDelz().toDouble();
      a_1 = extent_a[0] * sim.getDelx().toDouble();
      a_2 = extent_a[1] * sim.getDely().toDouble();

      b_3 = (pos )  * sim.getDelz().toDouble();
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

    baffleDisplay->drawRectangle(a_1, a_2, a_3, b_1, b_2, b_3, normal);
  }

  
  ui.VTKBaffles->update();
}
