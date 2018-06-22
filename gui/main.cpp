/*
 * main.cpp: main loop for interface
 */

#include <vtkOutputWindow.h>
#include <vtkFileOutputWindow.h>
 
#include "MainWindow.h"

int main(int argc, char **argv) {

    vtkOpenGLRenderWindow::SetGlobalMaximumNumberOfMultiSamples( 0 );
	QSurfaceFormat::setDefaultFormat(QVTKOpenGLWidget::defaultFormat());

	QApplication app(argc, argv);
    
#if defined(Q_OS_MAC)
    QDir dir(QApplication::applicationDirPath());
    dir.cdUp();
    dir.cd("PlugIns");
    QApplication::setLibraryPaths(QStringList(dir.absolutePath()));
#endif

    vtkFileOutputWindow *w = vtkFileOutputWindow::New();
    w->SetFileName("vtk_log.txt");
    vtkOutputWindow::SetInstance(w);
    w->Delete(); // now SetInstance owns the reference

    MainWindow window;

    window.show();

    return app.exec();
}
