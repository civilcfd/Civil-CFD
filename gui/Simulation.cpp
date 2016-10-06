/*
 * Simulation.cpp: implements the connection to the solver
 */

#include <QDebug>
#include <QFile>
#include <cstdlib>
#include <stdlib.h>
#include <QMessageBox>
#include <QByteArray>

#include "Simulation.h"

Simulation::Simulation() {
  solver = solver_init_empty();

  if(solver==NULL) {
    qDebug() << "solver failed to initialize in Simulation.cpp";
    throw;
  }

  ready = 0;
}

Simulation::~Simulation() {
  std::free(solver);

}

bool Simulation::loadSTL(QString file) {
  stlFile.setFile(file);
  return(stlFile.isFile());

}

int Simulation::loadSimulation(QString new_path) {
  if(!path.exists(new_path)) {
    qDebug() << "path does not exist " << new_path;
    return false;
  }
  QDir::setCurrent(new_path);
  
  setPath(new_path);

  ready=1;

  load();

  return true;
}

int Simulation::newSimulation(QString new_name, QString new_path) {
  if(!path.mkpath(new_path)) {
    qDebug() << "could not create directory " << new_path;
    return false;
  }
  if(!path.mkpath(new_path + "/0.000")) {
    qDebug() << "could not create directory " << new_path << "/0.000";
    return false;
  }
  if(!path.mkpath(new_path + "/vtk")) {
    qDebug() << "could not create directory " << new_path << "/vtk";
    return false;
  }
  QDir::setCurrent(new_path);
  
  setPath(new_path);
  setName(new_name);
  setAutot(true);

  track_empty();
  track_write();

  ready=1;

  save();

  return true;
}

int Simulation::save() {
  QFile file("casefile");
  if(!file.open(QIODevice::WriteOnly | QIODevice::Truncate | QIODevice::Text)) {
    qDebug() << "could not save to casefile";
    ready=0;
    return false;
  }
  QTextStream out(&file);
  
  out << name << "\n";
  out << stlFile.filePath() << "\n";
  file.close();

  if(write_solver(solver, "solverfile")) {
    qDebug() << "could not write solverfile";
    ready=0;
    return false;
  }
  // TESTING
  write_solver_xml(solver, "solver.xml");

  
  track_write();

  return true;
}

int Simulation::load() {


  clearSpecialBoundaries();
  clearBaffles();
  track_empty();
  std::free(solver->mesh);
  std::free(solver);
  solver = solver_init_empty();  
  
  QFile file("casefile");
  if(!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
    qDebug() << "could not open casefile";
    ready=0;
    return false;
  }
  QTextStream in(&file);
  QString stlFilename;

  name = in.readLine();
  stlFilename = in.readLine();
  file.close();

	vof_mpi_setup_solver(solver);

  read_mesh_xml(solver->mesh, "solver.xml");
  read_solver_xml(solver, "solver.xml");

  if(solver->turbulence_read("turbulencefile")) {
    qDebug() << "could not read turbulencefile";
    ready=0;
    return false;
  }
  
  
  strcpy(solver->ic[0].param,"end");
  /* if(read_initial(solver, "initials")) {
    qDebug() << "could not read initials";
    ready=0;
    return false;
  } */
  read_initial_xml(solver, "solver.xml");

	track_empty();
  track_read();

  loadSTL(stlFilename);

  return true;
}

void Simulation::setName(QString new_name) {
  name = new_name;
}

void Simulation::setPath(QString new_path) {
  path.setPath(new_path);
  
}

QString Simulation::getStlFilename() {
  return stlFile.filePath();
}

QString Simulation::getName() {
  return name;
}

QString Simulation::getPath() {
  return path.absolutePath();
}

QString Simulation::getGx() {
  return QString::number(solver->gx); 
}

QString Simulation::getGy() {
  return QString::number(solver->gy); 
}

QString Simulation::getGz() {
  return QString::number(solver->gz); 
}

QString Simulation::getNu() {
  return QString::number(solver->nu); 
}

QString Simulation::getRho() {
  return QString::number(solver->rho); 
}

bool Simulation::laminar() {
  if(solver->mesh->turbulence_model == NULL)
    return true;
  return false;
}

bool Simulation::kEpsilon() {
  if(!laminar()) {
    if(kE_check(solver))
      return true;
  }
 return false;
}

bool Simulation::parallelGMRES() {
	if(solver->pressure == vof_pressure_gmres_mpi) return true;
	else return false;
}

bool Simulation::GMRES() {
	return false;
}

bool Simulation::SOR() {
	return false;
}

bool Simulation::setImplicit(QString str) {
  if(str == "GMRES") 
    solver->pressure = vof_pressure_gmres_mpi;
  else if(str == "parallelGMRES") 
    solver->pressure = vof_pressure_gmres_mpi;
  else solver->pressure = vof_pressure_gmres_mpi;
  return true; 
}

QString Simulation::getRough() {
  double rough;
  struct kE_data *kE;
  
  if(kEpsilon()) {
    kE = (struct kE_data *) solver->mesh->turbulence_model;
    rough = kE->rough;
    return QString::number(rough); 
  }
  else return 0;
}

QString Simulation::getLength_scale() {
  double length_scale;
  struct kE_data *kE;
  
  if(kEpsilon()) {
    kE = (struct kE_data *) solver->mesh->turbulence_model;
    length_scale = kE->length_scale;
    return QString::number(length_scale); 
  }

  else return 0;
}

QString Simulation::getLength() {
  double length;
  struct kE_data *kE;
  
  if(kEpsilon()) {
    kE = (struct kE_data *) solver->mesh->turbulence_model;
    length = kE->raw_length;
    return QString::number(length); 
  }

  else return 0;
}

QString Simulation::domainLength() {
	double length;

	if(solver->mesh->delx * solver->mesh->imax < solver->mesh->dely * solver->mesh->jmax) 
		length = solver->mesh->delx * solver->mesh->imax;
	else length = solver->mesh->dely * solver->mesh->jmax;
	if(solver->mesh->delz * solver->mesh->kmax < length) 
		length = solver->mesh->delz * solver->mesh->kmax;
	
	return QString::number(length);
}

bool Simulation::setGx(QString str) {
  bool ok;
  double conv = str.toDouble(&ok);
  if(ok) solver->gx = conv;
  return ok;
}

bool Simulation::setGy(QString str) {
  bool ok;
  double conv = str.toDouble(&ok);
  if(ok) solver->gy = conv;
  return ok;
}

bool Simulation::setGz(QString str) {
  bool ok;
  double conv = str.toDouble(&ok);
  if(ok) solver->gz = conv;
  return ok;
}

bool Simulation::setNu(QString str) {
  bool ok;
  double conv = str.toDouble(&ok);
  if(ok) solver->nu = conv;
  return ok;
}

bool Simulation::setRho(QString str) {
  bool ok;
  double conv = str.toDouble(&ok);
  if(ok) solver->rho = conv;
  return ok;
}

bool Simulation::setTurbulence(QString str) {
  if(str == "kEpsilon") 
    kE_setup(solver);
  else laminar_setup(solver);
  return true; 
}

bool Simulation::setRough(QString str) {
  bool ok;
  double conv = str.toDouble(&ok);
  if(ok) {
    struct kE_data *kE;
  
    if(kEpsilon()) {
      kE = (struct kE_data *) solver->mesh->turbulence_model;
      kE->rough = conv;
    }
    else return false;
  }
  return ok;
}

bool Simulation::setLength_scale(QString str) {
  bool ok;
  double conv = str.toDouble(&ok);
  if(ok) {
    struct kE_data *kE;
  
    if(kEpsilon()) {
      kE = (struct kE_data *) solver->mesh->turbulence_model;
      kE->length_scale = conv;
    }
    else return false;
  }
  return ok;
}

bool Simulation::setLength(QString str) {
  bool ok;
  double conv = str.toDouble(&ok);
  if(ok) {
    struct kE_data *kE;
  
    if(kEpsilon()) {
      kE = (struct kE_data *) solver->mesh->turbulence_model;
      kE->raw_length = conv;
    }
    else return false;
  }
  return ok;
}

QString Simulation::getDelx() {
  return QString::number(solver->mesh->delx);
}

QString Simulation::getDely() {
  return QString::number(solver->mesh->dely);
}

QString Simulation::getDelz() {
  return QString::number(solver->mesh->delz);
}

QString Simulation::getImax() {
  return QString::number(solver->mesh->imax);
}

QString Simulation::getJmax() {
  return QString::number(solver->mesh->jmax);
}

QString Simulation::getKmax() {
  return QString::number(solver->mesh->kmax);
}

QString Simulation::getOrigin(int n) {
  return QString::number(solver->mesh->origin[n]);
}

QString Simulation::getInside(int n) {
  return QString::number(solver->mesh->inside[n]);
}

bool Simulation::setDelx(QString str) {
  bool ok;
  double conv = str.toDouble(&ok);
  if(ok) solver->mesh->delx = conv;
  return ok;
}

bool Simulation::setDely(QString str) {
  bool ok;
  double conv = str.toDouble(&ok);
  if(ok) solver->mesh->dely = conv;
  return ok;
}

bool Simulation::setDelz(QString str) {
  bool ok;
  double conv = str.toDouble(&ok);
  if(ok) solver->mesh->delz = conv;
  return ok;
}

bool Simulation::setImax(QString str) {
  bool ok;
  double conv = str.toDouble(&ok);
  if(ok) solver->mesh->imax = conv;
  return ok;
}

bool Simulation::setJmax(QString str) {
  bool ok;
  double conv = str.toDouble(&ok);
  if(ok) solver->mesh->jmax = conv;
  return ok;
}

bool Simulation::setKmax(QString str) {
  bool ok;
  double conv = str.toDouble(&ok);
  if(ok) solver->mesh->kmax = conv;
  return ok;
}

bool Simulation::setOrigin(int n, QString str) {
  bool ok;
  double conv = str.toDouble(&ok);
  if(ok) solver->mesh->origin[n] = conv;
  return ok;
}

bool Simulation::setInside(int n, QString str) {
  bool ok;
  double conv = str.toDouble(&ok);
  if(ok) solver->mesh->inside[n] = conv;
  return ok;
}

QString Simulation::getBoundaryText(int index) {

  switch(solver->mesh->wb[index]) {
  case 0:
    return "slip";
  case 1:
    return "no slip";
  case 2:
    return "continuous";
  }
  
  return "";
}

bool Simulation::setBoundaryText(int index, QString value) {
  enum wall_boundaries n;

  if(value == "slip") n = slip;
  else if(value == "no slip") n = no_slip;
  else if(value == "continuous") n = zero_gradient;
  else return false;

  solver->mesh->wb[index] = n;

  return true;
}

bool Simulation::addSpecialBoundary(int wall, 
  QString type,
  long int extent_a_1,
  long int extent_a_2, long int extent_b_1,
  long int extent_b_2, double value, double turbulence) {

  int n=0;
  if(type == "fixed velocity") n = 0;
  else if(type == "mass outflow") n = 1;
  else if(type == "hydraulic grade") n = 2;
  else if(type == "virtual weir") n = 3;

  if(mesh_sb_create(solver->mesh, wall, n, value, turbulence)==1)
    return false;

  mesh_sb_extent_a(solver->mesh, wall, extent_a_1, extent_a_2);
  mesh_sb_extent_b(solver->mesh, wall, extent_b_1, extent_b_2);

  return true;

}

void Simulation::editSpecialBoundary(QString type, long int extent_a_1, long int extent_a_2, long int extent_b_1, long int extent_b_2, double value, double turbulence) {
  sb_item->extent_a[0] = extent_a_1;
  sb_item->extent_a[1] = extent_a_2;
  sb_item->extent_b[0] = extent_b_1;
  sb_item->extent_b[1] = extent_b_2;
  sb_item->value = value;
  sb_item->turbulence = turbulence;
  

  if(type == "fixed velocity") sb_item->type = fixed_velocity;
  else if(type == "mass outflow") sb_item->type = mass_outflow;
  else if(type == "hydraulic grade") sb_item->type = hgl;
  else if(type == "virtual weir") sb_item->type = weir;

}

void Simulation::resetSpecialBoundary(int wall) {
  sb_item = solver->mesh->sb[wall];
}

bool Simulation::getNextSpecialBoundary(
      QString &type,
      long int (&extent_a)[2], long int (&extent_b)[2],
      double &value, double &turbulence) {

  if(sb_item == NULL) return false;

  extent_a[0] = sb_item->extent_a[0];
  extent_a[1] = sb_item->extent_a[1];
  extent_b[0] = sb_item->extent_b[0];
  extent_b[1] = sb_item->extent_b[1];
  
  value = sb_item->value;
  turbulence = sb_item->turbulence;

  type = "";
  switch(sb_item->type) {
  case 0:
    type = "fixed velocity";
    break;
  case 1:
    type = "mass outflow";
    break;
  case 2:
    type = "hydraulic grade";
    break;
  case 3:
    type = "virtual weir";
    break;
  }

  sb_item = sb_item->next;
  
  return true;
}

void Simulation::nextSpecialBoundary() {

  if(sb_item == NULL) return;

  sb_item = sb_item->next;
}

void Simulation::removeSpecialBoundary(int wall) {
  struct sb_data *sb_nitem;

  sb_nitem = solver->mesh->sb[wall];

  if(sb_nitem == sb_item) {
    solver->mesh->sb[wall] = sb_item->next;
    free(sb_item);
    return;
  }
  
  while(sb_nitem != NULL) {
    if(sb_nitem->next == sb_item) {
      sb_nitem->next = sb_item->next;
      free(sb_item);
      return;
    }

    sb_nitem = sb_nitem->next;
  }

}

void Simulation::clearSpecialBoundaries() {
	int i;
	
	for(i = 0; i<6; i++) {
		resetSpecialBoundary(i);
		while(sb_item != NULL) {
			removeSpecialBoundary(i);
			resetSpecialBoundary(i);
		}
	}

}

/* baffles */

bool Simulation::addBaffle(int wall, 
  QString type,
  long int extent_a_1,
  long int extent_a_2, long int extent_b_1,
  long int extent_b_2, double value, long int pos) {

  int n=0;
  if(type == "flow measurement") n = 0;
  else if(type == "barrier") n = 1;
  else if(type == "headloss") n = 2;
  else if(type == "swirl angle") n = 3;
  else if(type == "velocity deviation") n = 4;

  if(mesh_baffle_create(solver->mesh, wall, n, value, pos)==1)
    return false;

  mesh_baffle_extent_a(solver->mesh, wall, extent_a_1, extent_a_2);
  mesh_baffle_extent_b(solver->mesh, wall, extent_b_1, extent_b_2);
  
  return true;

}

void Simulation::editBaffle(QString type, long int extent_a_1, long int extent_a_2, long int extent_b_1, long int extent_b_2, double value, long int pos) {
  baffle->extent_a[0] = extent_a_1;
  baffle->extent_a[1] = extent_a_2;
  baffle->extent_b[0] = extent_b_1;
  baffle->extent_b[1] = extent_b_2;
  baffle->value = value;
  baffle->pos = pos;
  

  if(type == "flow measurement") baffle->type = flow;
  else if(type == "barrier") baffle->type = barrier;
  else if(type == "headloss") baffle->type = k;
  else if(type == "swirl angle") baffle->type = swirl_angle;
  else if(type == "velocity deviation") baffle->type = v_deviation;

}

void Simulation::resetBaffle(int wall) {
  baffle = solver->mesh->baffles[wall];
}

bool Simulation::getNextBaffle(
      QString &type,
      long int (&extent_a)[2], long int (&extent_b)[2],
      double &value, long int &pos) {

  if(baffle == NULL) return false;

  extent_a[0] = baffle->extent_a[0];
  extent_a[1] = baffle->extent_a[1];
  extent_b[0] = baffle->extent_b[0];
  extent_b[1] = baffle->extent_b[1];
  
  value = baffle->value;
  pos = baffle->pos;

  type = "";
  switch(baffle->type) {
  case 0:
    type = "flow measurement";
    break;
  case 1:
    type = "barrier";
    break;
  case 2:
    type = "headloss";
    break;
  case 3:
    type = "swirl angle";
    break;
  case 4:
    type = "velocity deviation";
    break;
  }

  baffle = baffle->next;
  
  return true;
}

void Simulation::nextBaffle() {

  if(baffle == NULL) return;

  baffle = baffle->next;
}

void Simulation::removeBaffle(int wall) {
  struct baffle_data *nbaffle;

  nbaffle = solver->mesh->baffles[wall];

  if(nbaffle == baffle) {
    solver->mesh->baffles[wall] = baffle->next;
    free(baffle);
    return;
  }
  
  while(nbaffle != NULL) {
    if(nbaffle->next == baffle) {
      nbaffle->next = baffle->next;
      free(baffle);
      return;
    }

    nbaffle = nbaffle->next;
  }

}

void Simulation::clearBaffles() {
	int i;
	
	for(i = 0; i<3; i++) {
		resetBaffle(i);
		while(baffle != NULL) {
			removeBaffle(i);
			resetBaffle(i);
		}
	}

}

QString Simulation::getInitialScalar(QString param) {
  QByteArray byteArray = param.toUtf8();
  const char* cString = byteArray.constData();
 
  return( QString::number(solver_get_initial_scalar(solver, const_cast<char*>(cString))) );
}

int Simulation::getInitialVector(QString param, double (&vector)[3]) {
  QByteArray byteArray = param.toUtf8();
  const char* cString = byteArray.constData();
  double _vector[3];
   
  if(solver_get_initial_vector(solver, const_cast<char*>(cString), _vector) == 1) {
    vector[0] = 0;
    vector[1] = 0;
    vector[2] = 0;   
    return 1;
  }
  
  vector[0] = _vector[0];
  vector[1] = _vector[1];
  vector[2] = _vector[2];
  
  return 0;
}

void Simulation::setInitial(QString param, QString s1, QString s2, QString s3) {
  double vector[3];
  bool ok;

  vector[0] = s1.toDouble(&ok);
  if(!ok) vector[0] = 0;
  vector[1] = s2.toDouble(&ok);
  if(!ok) vector[1] = 0;
  vector[2] = s3.toDouble(&ok);
  if(!ok) vector[2] = 0;

  QByteArray byteArray = param.toUtf8();
  const char* cString = byteArray.constData();
  
  solver_store_initial(solver, const_cast<char *>(cString), 3, vector);

}

QString Simulation::getT() {
  return QString::number(solver->t);
}

QString Simulation::getEndt() {
  return QString::number(solver->endt);
}

QString Simulation::getWritet() {
  return QString::number(solver->writet);
}

QString Simulation::getDelt() {
  return QString::number(solver->delt);
}

bool Simulation::getAutot() {
  if(solver->deltcal == NULL) return false;
  else return true;
}

bool Simulation::setT(QString str) {
  bool ok;
  double conv = str.toDouble(&ok);
  if(ok) solver->t = conv;
  return ok;
}

bool Simulation::setEndt(QString str) {
  bool ok;
  double conv = str.toDouble(&ok);
  if(ok) solver->endt = conv;
  return ok;
}

bool Simulation::setWritet(QString str) {
  bool ok;
  double conv = str.toDouble(&ok);
  if(ok) solver->writet = conv;
  return ok;
}

bool Simulation::setDelt(QString str) {
  bool ok;
  double conv = str.toDouble(&ok);
  if(ok) solver->delt = conv;
  return ok;
}

void Simulation::setAutot(bool autot) {

  if(autot) solver->deltcal = vof_mpi_deltcal;
  else solver->deltcal = NULL;
}

int Simulation::getTrackNext() {
  return track_next();
}

QString Simulation::getTrackT() {
  return QString::number(track_t(),'f',3);
}

void Simulation::trackRewind() {
  track_rewind();
}

QString Simulation::getTrackN(QString value) {
	int n;
	QString ret;
  track_rewind();

  while(value != getTrackT()) {
    if(getTrackNext() < 0) return QString::number(-1);
  }

	n = track_get_n();
  track_rewind();
  ret = QString::number(n);
  return ret;
}

bool removeDir(const QString & dirName)
{
    bool result = true;
    QDir dir(dirName);

    if (dir.exists(dirName)) {
        Q_FOREACH(QFileInfo info, dir.entryInfoList(QDir::NoDotAndDotDot | QDir::System | QDir::Hidden  | QDir::AllDirs | QDir::Files, QDir::DirsFirst)) {
            if (info.isDir()) {
                result = removeDir(info.absoluteFilePath());    
                if(!result) qDebug() << info.absoluteFilePath();

            }
            else {
                result = QFile::remove(info.absoluteFilePath());
                if(!result) qDebug() << info.absoluteFilePath();
            }

            if (!result) {
                return result;
            }
        }
        result = dir.rmdir(dirName);
                if(!result) qDebug() << dirName;
    }
          
    return result;
}
 
bool Simulation::deleteTrack(QString t) {
	int n = getTrackN(t).toInt();
	if(n == -1) return false;
	
	if(track_delete(n)) {
		track_write();
		removeDir(getPath() + "/" + t);
				
		QDir dir(getPath() + "/vtk");
		dir.setNameFilters(QStringList() << "*_" + QString::number(n) + ".vtk");
		dir.setFilter(QDir::Files);
		foreach(QString dirFile, dir.entryList())
		{
    	dir.remove(dirFile);
		}
		return true;
	}
	else return false;
}

