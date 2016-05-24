######################################################################
# Automatically generated by qmake (2.01a) Mon Jul 21 20:31:48 2014
######################################################################

TEMPLATE = app
TARGET =  ../bin/civil-cfd
DEPENDPATH += .

QMAKE_CXXFLAGS += -pthread -fopenmp
macx {
ICON = icon.icns
QMAKE_CXX = clang-omp++
QMAKE_LINK = clang-omp++
QMAKE_LINK_SHLIB = clang-omp++
}
unix:!macx {
QMAKE_CXX = g++
QMAKE_LINK = g++
QMAKE_LINK_SHLIB = g++
}
win32 {
RC_FILE = civilcfd.rc
QMAKE_CXX = cl.exe 
QMAKE_LINK = link.exe
QMAKE_LINK_SHLIB = link.exe
TARGET = ../../bin/gui
}
INCLUDEPATH += ../src/mesh3d ../src/solver3d . /usr/include/qt4 /usr/include/vtk-5.8 /opt/local/include/vtk-5.10 /usr/local/opt/vtk5/include/vtk-5.10 /usr/include/vtk-5.10
!win32 {
LIBS += -L/usr/local/lib -L/opt/local/lib -L/opt/local/lib/vtk-5.10 -fopenmp -L/usr/local/opt/vtk5/lib/vtk-5.10 -L../lib -lsolver3d -lmesh3d -lm -lqhull -lvtkCommon -lvtksys -lQVTK -lvtkViews -lvtkWidgets -lvtkInfovis -lvtkRendering -lvtkGraphics -lvtkImaging -lvtkIO -lvtkFiltering -lvtkDICOMParser -lvtkalglib -lvtkverdict -lvtkmetaio -lvtkexoIIc -lvtkftgl -lvtkHybrid -lvtkVolumeRendering -lz -lquazip -lpetsc -lmpi
CONFIG += debug x86_64
}
win32 {
INCLUDEPATH += "C:/qt-everywhere-opensource-src-4.8.7/lib" "C:/Program Files/VTK/include/vtk-5.10" "C:/SDKs/zlib-1.2.8/"
DEPENDPATH += ../bin "C:/Program Files/VTK/lib/vtk-5.10" "C:/qhull-2015.2/lib" "C:/qt-everywhere-opensource-src-4.8.7/lib"
LIBS += kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib opengl32.lib  ../bin/libsolver3d.lib ../bin/libmesh3d.lib "C:/qhull-2015.2/lib/qhullstatic.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkCommon.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtksys.lib" "C:/Program Files/VTK/lib/vtk-5.10/QVTK.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkViews.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkWidgets.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkInfovis.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkRendering.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkGraphics.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkImaging.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkIO.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkFiltering.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkDICOMParser.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkalglib.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkverdict.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkmetaio.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkexoIIc.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkftgl.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkHybrid.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkVolumeRendering.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkfreetype.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkjpeg.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkpng.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkftgl.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtktiff.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkexpat.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkzlib.lib" "C:/SDKs/zlib-1.2.8/zlib.lib"
CONFIG += release
# kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib 
}
# Input
HEADERS += MainWindow.h MeshDisplay.h GeometryDisplay.h Render.h RenderDisplay.h Boundary.h BoundaryDisplay.h SolverDialog.h VisualizeDisplay.h Baffle.h qcustomplot.h Visualize3DDisplay.h About.h InspectDialog.h Simulation.h
FORMS += civlcfd.ui render.ui boundary.ui sboundary.ui solver.ui inspect_cell.ui baffle.ui about.ui 
SOURCES += main.cpp MainWindow.cpp MeshDisplay.cpp GeometryDisplay.cpp Render.cpp RenderDisplay.cpp MainWindow_Boundaries.cpp Boundary.cpp BoundaryDisplay.cpp Solver.cpp ResultList.cpp Visualize.cpp VisualizeDisplay.cpp Baffle.cpp MainWindow_Baffles.cpp qcustomplot.cpp Visualize3D.cpp Visualize3DDisplay.cpp About.cpp InspectDialog.cpp Simulation.cpp

RESOURCES += \
    gui.qrc
