/*
 * main.cpp: main loop for interface
 */

#include "MainWindow.h"

int main(int argc, char **argv) {
    QApplication app(argc, argv);
    
#if defined(Q_OS_MAC)
QDir dir(QApplication::applicationDirPath());
dir.cdUp();
dir.cd("PlugIns");
QApplication::setLibraryPaths(QStringList(dir.absolutePath()));
#endif

    MainWindow window;

    window.show();

    return app.exec();
}
