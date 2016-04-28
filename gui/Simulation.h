/*
 * Simulation.h: connects interface to solver 
 */

#ifndef SIMULATION_H
#define SIMULATION_H

#include <QString>
#include <QDir>
#include <QFileInfo>

extern "C" {
#include "solver.h"
#include "readfile.h"
#include "readsolver.h"
#include "writesolver.h"
#include "kE.h"
#include "laminar.h"
#include "vof.h"
#include "track.h"
#include "vtk.h"
}

bool removeDir(const QString & dirName);

class Simulation {

public:
  Simulation();
  ~Simulation();

  void setName(QString new_name);
  void setPath(QString new_path);
  QString getName();
  QString getPath();

  int newSimulation(QString new_name, QString new_path);
  int loadSimulation(QString new_path);


  int save();
  int load();

  QString getT();
  QString getEndt();
  QString getWritet();
  bool getAutot();
  QString getDelt();

  bool setT(QString str);
  bool setEndt(QString str);
  bool setWritet(QString str);
  bool setDelt(QString str);
  void setAutot(bool autot);

  QString getGx();
  QString getGy();
  QString getGz();
  bool setGx(QString str);
  bool setGy(QString str);
  bool setGz(QString str);

  QString getNu();
  QString getRho();
  bool setNu(QString str);
  bool setRho(QString str);

  bool laminar();
  bool kEpsilon();
  bool setTurbulence(QString model);

  bool GMRES();
  bool parallelGMRES();
  bool SOR();
  bool setImplicit(QString model);
  
  QString getRough();
  QString getLength();
  QString domainLength();
  QString getLength_scale();
  bool setRough(QString str);
  bool setLength_scale(QString str);
  bool setLength(QString str);

  QString getDelx();
  QString getDely();
  QString getDelz();
  QString getImax();
  QString getJmax();
  QString getKmax();
  QString getOrigin(int n);
  
  bool setDelx(QString str);
  bool setDely(QString str);
  bool setDelz(QString str);
  bool setImax(QString str);
  bool setJmax(QString str);
  bool setKmax(QString str);
  bool setOrigin(int n, QString str);

  QString getInside(int n);
  bool setInside(int n, QString str);

  bool loadSTL(QString file);
  QString getStlFilename();

  QString getBoundaryText(int index);
  bool setBoundaryText(int index,
                          QString value);
  bool addSpecialBoundary(int wall,
    QString type,
    long int extent_a_1, long int extent_a_2,
    long int extent_b_1, long int extent_b_2,
    double value, double turbulence);
 
  bool getNextSpecialBoundary(
    QString &type,
    long int (&extent_a)[2], 
    long int (&extent_b)[2],
    double &value, double &turbulence);
  void nextSpecialBoundary();

  void resetSpecialBoundary(int wall);

  void removeSpecialBoundary(int wall);
	void clearSpecialBoundaries();

  void editSpecialBoundary(QString type, long int extent_a_1, long int extent_a_2, long int extent_b_1, long int extent_b_2, double value, double turbulence); 

  bool addBaffle(int wall,
    QString type,
    long int extent_a_1, long int extent_a_2,
    long int extent_b_1, long int extent_b_2,
    double value, long int pos);
 
  bool getNextBaffle(
    QString &type,
    long int (&extent_a)[2], 
    long int (&extent_b)[2],
    double &value, long int &pos);
  void nextBaffle();

  void resetBaffle(int wall);

  void removeBaffle(int wall);
	void clearBaffles();

  void editBaffle(QString type, long int extent_a_1, long int extent_a_2, long int extent_b_1, long int extent_b_2, double value, long int pos); 
  
  
  QString getInitialScalar(QString param);
  void setInitial(QString param, QString s1="", QString s2="", QString s3="");
  int getInitialVector(QString param, double (&vector)[3]);
  
  int getTrackNext();
  QString getTrackT();
  void trackRewind();
  QString getTrackN(QString str);
  bool deleteTrack(QString t);
  

private:
  struct solver_data *solver; 
  struct sb_data *sb_item;
  struct baffle_data *baffle;

  QFileInfo stlFile;
  QString name;
  QDir path;
  int ready;
};
#endif
