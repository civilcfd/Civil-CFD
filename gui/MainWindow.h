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
#include "Baffle.h"
#include "GeometryDisplay.h"
#include "Render.h"
#include "InspectDialog.h"
#include "BoundaryDisplay.h"
#include "VisualizeDisplay.h"
#include "Visualize3DDisplay.h"
#include "SolverDialog.h"
#include "About.h"

bool decompressFile(QString zip_filename , QString filename);

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
  void on_actionNew_triggered();
  void on_actionQuit_triggered();
  void on_actionAbout_triggered();
  void on_kEpsilon_toggled();
  void on_Laminar_toggled();
  void on_MeshParameters_itemDoubleClicked(QTreeWidgetItem *item, int column);
  void on_EditBoundary_clicked();
  void on_BoundaryTree_itemDoubleClicked(QTreeWidgetItem *item, int column);
  void on_AddSpecialBoundary_clicked();
  void on_BoundaryTree_itemActivated(QTreeWidgetItem *item, int column);
  void on_BoundaryTree_itemClicked(QTreeWidgetItem *item, int column);
  void on_RemoveSpecialBoundary_clicked();
  void on_InspectCell_clicked();
  
  void on_EditBaffle_clicked();
  void on_BaffleTree_itemDoubleClicked(QTreeWidgetItem *item, int column);
  void on_AddBaffle_clicked();
  void on_BaffleTree_itemActivated(QTreeWidgetItem *item, int column);
  void on_BaffleTree_itemClicked(QTreeWidgetItem *item, int column);
  void on_RemoveBaffle_clicked();
  
  void on_updateRange_toggled();
  void on_saveJPEG_clicked();
  void on_saveJPEG3d_clicked();
  
  void on_inside_x_textChanged();
  void on_inside_y_textChanged();
  void on_inside_z_textChanged();
  void on_from_textChanged();
  void on_to_textChanged();
  void on_origin_valueChanged();
  void on_origin3d_valueChanged();
  void on_extent3d_valueChanged();
  void on_timesteps_currentItemChanged();
  void on_timesteps3d_currentItemChanged();
  void on_xNormal_toggled();
  void on_yNormal_toggled();
  void on_zNormal_toggled();
  void on_xNormal3d_toggled();
  void on_yNormal3d_toggled();
  void on_zNormal3d_toggled();
  void on_hideBoundaryMesh_toggled();
  void on_hideBoundaryDomain_toggled();
  void on_hideBaffleMesh_toggled();
  void on_hideBaffleDomain_toggled();
  void on_contourVOF_toggled();
  void on_contourP_toggled();
  void on_contourK_toggled();
  void on_showLegend_toggled();
  void on_showAxis_toggled();
  void on_showAxis3d_toggled();
  void on_contourVorticity_toggled();
  void on_showVectors_toggled();
  void on_blockObstacles_toggled();
  void on_showMesh_toggled();
  void on_showMesh3d_toggled();
  void on_blockObstacles3d_toggled();
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
  void editBaffle(QTreeWidgetItem *item);
  void disableAll();
  void enableAll();
  void toggle();
  void update();
  void boundariesUpdate();
  void buildBoundaryList();
  void bafflesUpdate();
  void buildBaffleList();
  void meshUpdate();
  void visualizeUpdate();
  void visualize3dUpdate();
  int  buildTimesteps();
  int  build3dTimesteps();
  void visualizeRender();
  void visualize3dRender();
  void buildResultList();
  void write();
  void movePoint();
  void updateSlider();
  void updateSlider3d();
  Ui::MainWindow ui;
  Simulation sim;
  MeshDisplay *meshDisplay;
  GeometryDisplay *geometryDisplay;
  VisualizeDisplay *visualizeDisplay;
  Visualize3DDisplay *visualize3dDisplay;
  RenderDialog *renderDialog;
  BoundaryDialog *boundaryDialog;
  SBoundaryDialog *sBoundaryDialog;
  BoundaryDisplay *boundaryDisplay;
  BaffleDialog *baffleDialog;
  BoundaryDisplay *baffleDisplay;
  SolverDialog *solverDialog;
  AboutDialog *aboutDialog;
  InspectDialog *inspectDialog;
  QString appPath;
};

#endif

