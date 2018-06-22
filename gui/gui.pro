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
QMAKE_CXX = mpic++
QMAKE_LINK = mpic++
QMAKE_LINK_SHLIB = mpic++
}
win32 {
CONFIG += console
RC_FILE = civilcfd.rc
QMAKE_CXX = cl.exe 
QMAKE_LINK = link.exe
QMAKE_LINK_SHLIB = link.exe
TARGET = ../../bin/gui
}
INCLUDEPATH += ../src/mesh3d ../src/solver3d . /usr/include/qt4 /usr/include/vtk-5.8  /usr/include/vtk-6.3 /opt/local/include/vtk-5.10 /usr/local/opt/vtk5/include/vtk-5.10 /usr/include/vtk-5.10 /usr/include/mpi /usr/include/petsc /usr/local/include /usr/include/libxml2 /usr/local/opt/libxml2/include/libxml2
!win32 {
LIBS += -L/usr/local/opt/libxml2/lib -L/usr/local/lib -L/opt/local/lib -L/opt/local/lib/vtk-5.10 -L/usr/lib/x86_64-linux-gnu -fopenmp -L/usr/local/opt/vtk5/lib/vtk-5.10 -L../lib -lsolver3d -lmesh3d -lm -lqhull -lvtkCommonCore-6.3 -lvtksys-6.3 -lvtkViewsCore-6.3 -lvtkRenderingCore-6.3 -lvtkImagingCore-6.3 -lvtkIOCore-6.3 -lvtkFiltersCore-6.3 -lvtkDICOMParser-6.3 -lvtkalglib-6.3 -lvtkverdict-6.3 -lvtkmetaio-6.3 -lvtkexoIIc-6.3 -lvtkftgl-6.3 -lvtkRenderingVolume-6.3 -lz -lquazip -lpetsc -lmpi -lxml2
CONFIG += debug x86_64
}
win32 {
INCLUDEPATH += "C:/qt/4.8.6/lib" "C:/Program Files/VTK/include/vtk-5.10" "C:/SDKs/zlib-1.2.8/" "C:/SDKs/libxml2-2.9.4/include" "C:/MSMPI/Include" "C:/SDKs/libiconv/include" "C:/quazip-0.7.2"
DEPENDPATH += ../bin "C:/Program Files/VTK/lib/vtk-5.10" "C:/qhull-2015.2/lib" "C:/qt/4.8.6/lib" "C:/MSMPI/Lib/x64" "C:/SDKs/libxml2-2.9.4/lib" "C:/SDKs/libiconv" "C:/quazip-0.7.2/quazip/release"
LIBS += kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib opengl32.lib  ../bin/libsolver3d.lib ../bin/libmesh3d.lib "C:/qhull-2015.2/lib/qhullstatic.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkCommon.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtksys.lib" "C:/Program Files/VTK/lib/vtk-5.10/QVTK.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkViews.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkWidgets.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkInfovis.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkRendering.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkGraphics.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkImaging.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkIO.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkFiltering.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkDICOMParser.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkalglib.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkverdict.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkmetaio.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkexoIIc.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkftgl.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkHybrid.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkVolumeRendering.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkfreetype.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkjpeg.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkpng.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkftgl.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtktiff.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkexpat.lib" "C:/Program Files/VTK/lib/vtk-5.10/vtkzlib.lib" "C:/SDKs/zlib-1.2.8/zlib.lib" "C:/MSMPI/Lib/x64/msmpi.lib" "C:/SDKs/libxml2-2.9.4/lib/libxml2.lib" "C:/petsc-3.7.3/arch-mswin-c-opt/lib/libpetsc.lib" "C:/petsc-3.7.3/arch-mswin-c-opt/lib/libf2cblas.lib" "C:/petsc-3.7.3/arch-mswin-c-opt/lib/libf2clapack.lib"  "C:/SDKs/libiconv/libiconv.lib" "C:/quazip-0.5.1/quazip/release/quazip.lib" ws2_32.lib
CONFIG += release
QMAKE_LFLAGS += /NODEFAULTLIB:LIBCMT 
# kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib 
}
# Input
HEADERS += MainWindow.h MeshDisplay.h GeometryDisplay.h Render.h RenderDisplay.h Boundary.h BoundaryDisplay.h SolverDialog.h VisualizeDisplay.h Baffle.h qcustomplot.h Visualize3DDisplay.h About.h InspectDialog.h Simulation.h
FORMS += civlcfd.ui render.ui boundary.ui sboundary.ui solver.ui inspect_cell.ui baffle.ui about.ui 
SOURCES += main.cpp MainWindow.cpp MeshDisplay.cpp GeometryDisplay.cpp Render.cpp RenderDisplay.cpp MainWindow_Boundaries.cpp Boundary.cpp BoundaryDisplay.cpp Solver.cpp ResultList.cpp Visualize.cpp VisualizeDisplay.cpp Baffle.cpp MainWindow_Baffles.cpp qcustomplot.cpp Visualize3D.cpp Visualize3DDisplay.cpp About.cpp InspectDialog.cpp Simulation.cpp

RESOURCES += \
    gui.qrc
