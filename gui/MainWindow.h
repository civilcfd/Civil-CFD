/*
 * MainWindow.h : implements the ui from designer
 */

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QDialog>
#include "ui_civlcfd.h"

#include "Simulation.h"
#include "MeshDisplay.h"
#include "Boundary.h"
#include "GeometryDisplay.h"
#include "Render.h"
#include "BoundaryDisplay.h"
#include "VisualizeDisplay.h"
#include "SolverDialog.h"

class MainWindow : public QMainWindow {
  Q_OBJECT

public:
  MainWindow();

private slots:
  void on_RunSimulation_clicked();
  void on_STLRender_clicked();
  void on_STLOpen_clicked();
  void on_MeshUndo_clicked();
  void on_MeshUpdate_clicked();
  void on_New_clicked();
  void on_Save_clicked();
  void on_actionOpen_triggered();
  void on_actionSave_triggered();
  void on_kEpsilon_toggled();
  void on_Laminar_toggled();
  void on_MeshParameters_itemDoubleClicked(QTreeWidgetItem *item, int column);
  void on_EditBoundary_clicked();
  void on_BoundaryTree_itemDoubleClicked(QTreeWidgetItem *item, int column);
  void on_AddSpecialBoundary_clicked();
  void on_BoundaryTree_itemActivated(QTreeWidgetItem *item, int column);
  void on_BoundaryTree_itemClicked(QTreeWidgetItem *item, int column);
  void on_RemoveSpecialBoundary_clicked();
  void on_inside_x_textChanged();
  void on_inside_y_textChanged();
  void on_inside_z_textChanged();
  void on_origin_valueChanged();
  void on_timesteps_currentItemChanged();
  void on_xNormal_toggled();
  void on_yNormal_toggled();
  void on_zNormal_toggled();
  void on_contourVOF_toggled();
  void on_contourP_toggled();
  void on_contourK_toggled();
  void on_showVectors_toggled();
  void on_blockObstacles_toggled();
  void on_showMesh_toggled();
  void on_earthGravity_clicked();
  void on_water20C_clicked();
  void on_defaultLength_clicked();
  void on_calcRough_clicked();
  void on_SelectAll_clicked();
  void on_Clear_clicked();
  void on_Delete_clicked();
  
private:
  bool saveNotify();
  void editBoundary(QTreeWidgetItem *item);
  void disableAll();
  void enableAll();
  void toggle();
  void update();
  void boundariesUpdate();
  void meshUpdate();
  void visualizeUpdate();
  void buildTimesteps();
  void visualizeRender();
  void buildResultList();
  void write();
  void movePoint();
  Ui::MainWindow ui;
  Simulation sim;
  MeshDisplay *meshDisplay;
  GeometryDisplay *geometryDisplay;
  VisualizeDisplay *visualizeDisplay;
  RenderDialog *renderDialog;
  BoundaryDialog *boundaryDialog;
  SBoundaryDialog *sBoundaryDialog;
  BoundaryDisplay *boundaryDisplay;
  SolverDialog *solverDialog;
  QString appPath;
};

#endif

