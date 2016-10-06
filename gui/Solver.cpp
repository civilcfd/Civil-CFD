/*
 * Solver.cpp
 *
 * Implementation of solver dialog
 */

#include "SolverDialog.h"
#include "qcustomplot.h"


SolverDialog::SolverDialog(Simulation &sim, QString appPath, QString t, QString prefix) {
  QString cmd;
  bool ok;

  stopT = sim.getEndt().toDouble();
  stopped = false;

  ui.setupUi(this);

  ui.status->setText("Simulation Running");
  ui.Return->setEnabled(false);

  process = new QProcess(this);

  connect(process,
            SIGNAL(error(QProcess::ProcessError)),
            SLOT(error(QProcess::ProcessError)));
  connect(process,
            SIGNAL(finished(int, QProcess::ExitStatus)),
            SLOT(finished(int, QProcess::ExitStatus)));
  connect(process,
            SIGNAL(stateChanged(QProcess::ProcessState)),
            SLOT(stateChanged(QProcess::ProcessState)));
  connect(process,
            SIGNAL(readyReadStandardOutput()),
            SLOT(readyReadStandardOutput()));

#ifdef _WIN32
  cmd = appPath + "/solver3d.exe";
#else
  cmd = appPath + "/solver3d";
#endif
  
  mint=0;
  
  if(t.toDouble(&ok) < 0.00001) t = "";
  if(!ok) t="";
  
  
  if(!QFile::exists(cmd)) {
    ui.status->setText("Could not find executable: " + cmd);
    ui.Return->setEnabled(true);
  }
  else {
    if(t != "") {
      cmd = cmd + " " + t;
      mint = t.toDouble();
    }
    cmd = prefix + " " + cmd;
    process->start(cmd);
  }
  
  progressVal = mint;
   
  ui.plot->addGraph();
  ui.plot->graph(0)->setPen(QPen(Qt::blue));
  ui.plot->graph(0)->setBrush(QBrush(QColor(240, 255, 200)));
  ui.plot->xAxis->setRange(mint,mint+1);
  ui.plot->xAxis->setLabel("Timestep (s)");
  ui.plot->yAxis->setRange(0,0.5);
  ui.plot->yAxis->setLabel("Timestep size (s)");
  connect(ui.plot->xAxis, SIGNAL(rangeChanged(QCPRange)), ui.plot->xAxis2, SLOT(setRange(QCPRange)));
  connect(ui.plot->yAxis, SIGNAL(rangeChanged(QCPRange)), ui.plot->yAxis2, SLOT(setRange(QCPRange)));
  
  ui.flow->addGraph();
  ui.flow->graph(0)->setPen(QPen(Qt::blue));
  ui.flow->xAxis->setRange(mint,mint+1);
  ui.flow->xAxis->setLabel("Timestep (s)");
  ui.flow->yAxis->setRange(0,0.5);
  ui.flow->yAxis->setLabel("Flow (L/s)");
  connect(ui.flow->xAxis, SIGNAL(rangeChanged(QCPRange)), ui.flow->xAxis2, SLOT(setRange(QCPRange)));
  connect(ui.flow->yAxis, SIGNAL(rangeChanged(QCPRange)), ui.flow->yAxis2, SLOT(setRange(QCPRange)));
/*
  ui.plot->plotLayout()->insertRow(0);
  ui.plot->plotLayout()->addElement(0, 0, new QCPPlotTitle(ui.plot, "Timestep Size")); */
  
  maxDelt = 0;
  maxFlow = 1;
  minFlow = 0;
  nFlow   = 0;
}

void SolverDialog::on_Return_clicked() {
  close(); 
}

void SolverDialog::on_Stop_clicked() {

  stopped = true;
  
#ifdef _WIN32
  process->kill();
#else
  process->terminate();
#endif

}

void SolverDialog::readyReadStandardOutput() {
  QString str = process->readAllStandardOutput(); 
  ui.output->appendPlainText(str);
  int n;
  bool ok;
  double val;

  QFont legendFont;
  
  if(str.contains("timestep") || str.contains("delt")) {

    QStringList list = str.split(" ");
    for(int i=0; i < list.count()-2; i++) {
      if(list[i].contains("timestep")) {
      	
        val = list[i+1].toDouble(&ok);
        if(ok) progressVal = val; 
        else continue;
      }
      
      if(list[i].contains("delt")) {
        val = list[i+1].toDouble(&ok);
        if(ok) delt = val; 
        else continue;
        if(list[i+2] != "|") continue;
        if(delt > maxDelt) maxDelt = delt;
        
        ui.plot->graph(0)->addData(progressVal, delt);
      }
      
      if(list[i].contains("Flow")) {
        if(i+4 >= list.count()-3) continue;
        if(list[i+5] != "L/s") continue;
        
        flow = list[i+4].toDouble(&ok);
        if(!ok) continue;
        if(flow > maxFlow) maxFlow = flow;
        if(flow < minFlow) minFlow = flow;
        
        
        n = list[i+2].toInt(&ok);
        if(!ok) continue;
        if(n>8) continue;
        if(n > nFlow) {
          nFlow = n;
          ui.flow->addGraph();
          
          
          switch(n) {
          case 1:
          ui.flow->graph(n)->setPen(QPen(Qt::red));
          
          ui.flow->legend->setVisible(true);
          legendFont = font();
          legendFont.setPointSize(9);
          ui.flow->legend->setFont(legendFont);
          ui.flow->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignLeft|Qt::AlignTop);
          
          ui.flow->graph(0)->setName("Baffle 0");
          ui.flow->graph(1)->setName("Baffle 1");
          break;
          case 2:
          ui.flow->graph(n)->setPen(QPen(Qt::green));
          ui.flow->graph(n)->setName("Baffle 2");
          break;
          case 3:
          ui.flow->graph(n)->setPen(QPen(Qt::black));
          ui.flow->graph(n)->setName("Baffle 3");
          break;
          case 4:
          ui.flow->graph(n)->setPen(QPen(Qt::cyan));
          ui.flow->graph(n)->setName("Baffle 4");
          break;
          case 5:
          ui.flow->graph(n)->setPen(QPen(Qt::darkGray));
          ui.flow->graph(n)->setName("Baffle 5");
          break;
          case 6:
          ui.flow->graph(n)->setPen(QPen(Qt::darkYellow));
          ui.flow->graph(n)->setName("Baffle 6");
          break;
          case 7:
          ui.flow->graph(n)->setPen(QPen(Qt::darkRed));
          ui.flow->graph(n)->setName("Baffle 7");
          break;
          case 8:
          ui.flow->graph(n)->setPen(QPen(Qt::darkGreen));
          ui.flow->graph(n)->setName("Baffle 8");
          break;
          }
        }
        
        ui.flow->graph(n)->addData(progressVal, flow);
        
      }
    }

    if(progressVal > 0) {
      ui.progressBar->setValue(progressVal * 100 / stopT);
    }
  }
  
  ui.plot->xAxis->setRange(mint, progressVal);
  ui.plot->yAxis->setRange(0, maxDelt);
  ui.plot->replot();
  ui.flow->xAxis->setRange(mint, progressVal);
  ui.flow->yAxis->setRange(minFlow, maxFlow);
  ui.flow->replot();

}

void SolverDialog::error(QProcess::ProcessError error) {

  if(stopped == true) {
    ui.status->setText("Stopped by user");  
  } else if(error == QProcess::FailedToStart) {
    ui.status->setText("Error: solver3d could not start");
  } else if(error == QProcess::Crashed) {
    ui.status->setText("Error: solver3d has crashed");
  }

  ui.Return->setEnabled(true);
}

void SolverDialog::finished(int exitCode, QProcess::ExitStatus status) {
  ui.progressBar->setValue(100);

  if(exitCode == 0) {
    if(stopped == false) ui.status->setText("Finished successfully");
    else ui.status->setText("Stopped by user");
  }
  
  if(exitCode == 1) {
    if(stopped == false) ui.status->setText("Process failed");
    else ui.status->setText("Stopped by user");
  }

  ui.Return->setEnabled(true);
}

void SolverDialog::stateChanged(QProcess::ProcessState state) {
}
