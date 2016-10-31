/********************************************************************************
** Form generated from reading UI file 'civlcfd.ui'
**
** Created by: Qt User Interface Compiler version 4.8.7
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_CIVLCFD_H
#define UI_CIVLCFD_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QComboBox>
#include <QtGui/QFrame>
#include <QtGui/QGridLayout>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QListWidget>
#include <QtGui/QMainWindow>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QPlainTextEdit>
#include <QtGui/QPushButton>
#include <QtGui/QRadioButton>
#include <QtGui/QSlider>
#include <QtGui/QSpacerItem>
#include <QtGui/QStatusBar>
#include <QtGui/QTabWidget>
#include <QtGui/QTableWidget>
#include <QtGui/QTreeWidget>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>
#include "QVTKWidget.h"

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QAction *actionNew;
    QAction *actionOpen;
    QAction *actionSave;
    QAction *actionQuit;
    QAction *actionAbout;
    QWidget *centralwidget;
    QGridLayout *gridLayout_3;
    QTabWidget *tabWidget;
    QWidget *Information;
    QVBoxLayout *verticalLayout;
    QGridLayout *gridLayout;
    QPlainTextEdit *ElementNumber;
    QPlainTextEdit *Dims;
    QLabel *label_2;
    QListWidget *ResultList;
    QFrame *line;
    QLabel *label_5;
    QPushButton *New;
    QPlainTextEdit *ElementSize;
    QSpacerItem *verticalSpacer;
    QSpacerItem *verticalSpacer_2;
    QPushButton *Save;
    QPushButton *Clear;
    QPlainTextEdit *Name;
    QLabel *label_9;
    QPushButton *InspectCell;
    QPushButton *Delete;
    QPushButton *SelectAll;
    QLabel *label_10;
    QSpacerItem *verticalSpacer_4;
    QLabel *label;
    QPlainTextEdit *STLFile;
    QLabel *label_4;
    QSpacerItem *horizontalSpacer_3;
    QLabel *label_6;
    QPlainTextEdit *Path;
    QSpacerItem *horizontalSpacer_2;
    QSpacerItem *verticalSpacer_3;
    QLabel *label_3;
    QLabel *label_7;
    QWidget *Models;
    QGridLayout *gridLayout_2;
    QGridLayout *gridLayout_4;
    QPlainTextEdit *gx;
    QPushButton *calcRough;
    QPlainTextEdit *rho;
    QSpacerItem *horizontalSpacer_7;
    QRadioButton *kEpsilon;
    QPlainTextEdit *length_scale;
    QSpacerItem *horizontalSpacer_8;
    QLabel *label_16;
    QLabel *label_20;
    QLabel *label_14;
    QSpacerItem *verticalSpacer_5;
    QLabel *label_18;
    QRadioButton *Laminar;
    QPlainTextEdit *gz;
    QSpacerItem *verticalSpacer_6;
    QLabel *label_11;
    QLabel *label_19;
    QSpacerItem *horizontalSpacer_4;
    QFrame *line_2;
    QPushButton *water20C;
    QPushButton *defaultLength;
    QLabel *label_15;
    QLabel *label_12;
    QLabel *label_13;
    QSpacerItem *verticalSpacer_7;
    QPlainTextEdit *gy;
    QPlainTextEdit *nu;
    QSpacerItem *horizontalSpacer_5;
    QLabel *label_17;
    QPlainTextEdit *rough;
    QSpacerItem *horizontalSpacer;
    QPushButton *earthGravity;
    QPlainTextEdit *length;
    QLabel *label_8;
    QPushButton *domainLength;
    QWidget *Mesh;
    QGridLayout *gridLayout_6;
    QGridLayout *gridLayout_5;
    QWidget *widget;
    QGridLayout *gridLayout_7;
    QPushButton *MeshUndo;
    QVTKWidget *VTKMesh;
    QTreeWidget *MeshParameters;
    QPushButton *MeshUpdate;
    QWidget *Geometry;
    QGridLayout *gridLayout_10;
    QGridLayout *gridLayout_9;
    QLabel *label_22;
    QPlainTextEdit *STLFile_2;
    QVTKWidget *VTKGeometry;
    QLabel *label_25;
    QLabel *label_24;
    QLabel *label_23;
    QPlainTextEdit *inside_x;
    QLabel *label_21;
    QPlainTextEdit *inside_z;
    QPlainTextEdit *inside_y;
    QSpacerItem *horizontalSpacer_6;
    QSpacerItem *verticalSpacer_8;
    QPushButton *STLRender;
    QSpacerItem *horizontalSpacer_9;
    QPushButton *STLOpen;
    QWidget *Boundaries;
    QGridLayout *gridLayout_8;
    QGridLayout *gridLayout_12;
    QTreeWidget *BoundaryTree;
    QHBoxLayout *horizontalLayout;
    QPushButton *EditBoundary;
    QPushButton *AddSpecialBoundary;
    QPushButton *RemoveSpecialBoundary;
    QVTKWidget *VTKBoundaries;
    QHBoxLayout *horizontalLayout_6;
    QSpacerItem *horizontalSpacer_24;
    QCheckBox *hideBoundaryMesh;
    QSpacerItem *horizontalSpacer_26;
    QCheckBox *hideBoundaryDomain;
    QSpacerItem *horizontalSpacer_25;
    QWidget *Baffles;
    QGridLayout *gridLayout_20;
    QGridLayout *gridLayout_19;
    QHBoxLayout *horizontalLayout_2;
    QPushButton *EditBaffle;
    QPushButton *AddBaffle;
    QPushButton *RemoveBaffle;
    QVTKWidget *VTKBaffles;
    QTreeWidget *BaffleTree;
    QHBoxLayout *horizontalLayout_5;
    QSpacerItem *horizontalSpacer_21;
    QCheckBox *hideBaffleMesh;
    QSpacerItem *horizontalSpacer_23;
    QCheckBox *hideBaffleDomain;
    QSpacerItem *horizontalSpacer_22;
    QWidget *InitialConditions;
    QGridLayout *gridLayout_14;
    QGridLayout *gridLayout_13;
    QSpacerItem *verticalSpacer_12;
    QLabel *label_26;
    QSpacerItem *horizontalSpacer_10;
    QFrame *line_5;
    QLabel *label_30;
    QSpacerItem *horizontalSpacer_12;
    QLabel *label_29;
    QFrame *line_4;
    QSpacerItem *verticalSpacer_11;
    QPlainTextEdit *initialU;
    QSpacerItem *horizontalSpacer_11;
    QPlainTextEdit *initialW;
    QLabel *label_28;
    QPlainTextEdit *height;
    QSpacerItem *verticalSpacer_10;
    QLabel *label_32;
    QFrame *line_3;
    QLabel *label_27;
    QPlainTextEdit *initialV;
    QLabel *label_31;
    QSpacerItem *horizontalSpacer_13;
    QCheckBox *hydrostatic;
    QSpacerItem *verticalSpacer_9;
    QLabel *label_33;
    QSpacerItem *verticalSpacer_13;
    QLabel *label_34;
    QPlainTextEdit *initialk;
    QLabel *label_43;
    QSpacerItem *verticalSpacer_24;
    QTableWidget *fillPoints;
    QWidget *Simulate;
    QGridLayout *gridLayout_16;
    QGridLayout *gridLayout_15;
    QLabel *label_37;
    QSpacerItem *horizontalSpacer_15;
    QSpacerItem *verticalSpacer_16;
    QComboBox *t;
    QSpacerItem *verticalSpacer_15;
    QSpacerItem *horizontalSpacer_16;
    QLabel *label_39;
    QSpacerItem *verticalSpacer_14;
    QSpacerItem *verticalSpacer_29;
    QFrame *line_6;
    QCheckBox *autot;
    QSpacerItem *verticalSpacer_17;
    QPlainTextEdit *writet;
    QPushButton *RunSimulation;
    QPlainTextEdit *delt;
    QLabel *label_38;
    QSpacerItem *horizontalSpacer_14;
    QLabel *label_47;
    QLabel *label_36;
    QPlainTextEdit *endt;
    QLabel *label_35;
    QFrame *line_7;
    QLabel *label_49;
    QComboBox *processes;
    QWidget *Visualize;
    QGridLayout *gridLayout_18;
    QGridLayout *gridLayout_17;
    QHBoxLayout *horizontalLayout_3;
    QSpacerItem *horizontalSpacer_17;
    QCheckBox *updateRange;
    QPlainTextEdit *from;
    QLabel *label_48;
    QPlainTextEdit *to;
    QSpacerItem *horizontalSpacer_18;
    QSpacerItem *verticalSpacer_19;
    QSpacerItem *verticalSpacer_21;
    QLabel *label_41;
    QSlider *origin;
    QRadioButton *xNormal;
    QSpacerItem *verticalSpacer_22;
    QLabel *label_44;
    QRadioButton *contourVOF;
    QSpacerItem *verticalSpacer_20;
    QLabel *label_45;
    QFrame *line_8;
    QLabel *label_42;
    QRadioButton *zNormal;
    QRadioButton *yNormal;
    QGridLayout *gridLayout_11;
    QRadioButton *contourVorticity;
    QRadioButton *contourVelocity;
    QLabel *originText;
    QRadioButton *contourP;
    QVTKWidget *vis;
    QRadioButton *contourK;
    QPushButton *saveJPEG;
    QListWidget *timesteps;
    QCheckBox *showMesh;
    QSpacerItem *verticalSpacer_18;
    QCheckBox *showVectors;
    QCheckBox *blockObstacles;
    QLabel *label_46;
    QLabel *label_40;
    QCheckBox *showAxis;
    QCheckBox *showLegend;
    QWidget *Visualize3d;
    QGridLayout *gridLayout_23;
    QGridLayout *gridLayout_21;
    QCheckBox *showAxis3d;
    QLabel *extent3dText;
    QVTKWidget *vis3d;
    QLabel *origin3dText;
    QLabel *label_54;
    QCheckBox *blockObstacles3d;
    QCheckBox *showMesh3d;
    QPushButton *saveJPEG3d;
    QRadioButton *yNormal3d;
    QSpacerItem *verticalSpacer_27;
    QSpacerItem *verticalSpacer_23;
    QSpacerItem *verticalSpacer_25;
    QRadioButton *xNormal3d;
    QSpacerItem *verticalSpacer_26;
    QLabel *label_51;
    QSlider *origin3d;
    QLabel *label_53;
    QFrame *line_9;
    QHBoxLayout *horizontalLayout_4;
    QSpacerItem *horizontalSpacer_19;
    QSpacerItem *horizontalSpacer_20;
    QRadioButton *zNormal3d;
    QListWidget *timesteps3d;
    QSpacerItem *verticalSpacer_28;
    QLabel *label_55;
    QLabel *label_56;
    QLabel *label_52;
    QSlider *extent3d;
    QMenuBar *menubar;
    QMenu *menuFile;
    QMenu *menuHelp;
    QStatusBar *statusbar;
    QButtonGroup *buttonGroup_2;
    QButtonGroup *buttonGroup;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->setEnabled(true);
        MainWindow->resize(950, 751);
        QSizePolicy sizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(MainWindow->sizePolicy().hasHeightForWidth());
        MainWindow->setSizePolicy(sizePolicy);
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/icons/icon.svg"), QSize(), QIcon::Normal, QIcon::Off);
        MainWindow->setWindowIcon(icon);
        actionNew = new QAction(MainWindow);
        actionNew->setObjectName(QString::fromUtf8("actionNew"));
        actionOpen = new QAction(MainWindow);
        actionOpen->setObjectName(QString::fromUtf8("actionOpen"));
        actionSave = new QAction(MainWindow);
        actionSave->setObjectName(QString::fromUtf8("actionSave"));
        actionQuit = new QAction(MainWindow);
        actionQuit->setObjectName(QString::fromUtf8("actionQuit"));
        actionAbout = new QAction(MainWindow);
        actionAbout->setObjectName(QString::fromUtf8("actionAbout"));
        centralwidget = new QWidget(MainWindow);
        centralwidget->setObjectName(QString::fromUtf8("centralwidget"));
        centralwidget->setEnabled(true);
        QSizePolicy sizePolicy1(QSizePolicy::Expanding, QSizePolicy::Expanding);
        sizePolicy1.setHorizontalStretch(100);
        sizePolicy1.setVerticalStretch(100);
        sizePolicy1.setHeightForWidth(centralwidget->sizePolicy().hasHeightForWidth());
        centralwidget->setSizePolicy(sizePolicy1);
        gridLayout_3 = new QGridLayout(centralwidget);
        gridLayout_3->setObjectName(QString::fromUtf8("gridLayout_3"));
        tabWidget = new QTabWidget(centralwidget);
        tabWidget->setObjectName(QString::fromUtf8("tabWidget"));
        tabWidget->setEnabled(true);
        sizePolicy1.setHeightForWidth(tabWidget->sizePolicy().hasHeightForWidth());
        tabWidget->setSizePolicy(sizePolicy1);
        tabWidget->setMinimumSize(QSize(850, 0));
        tabWidget->setMaximumSize(QSize(2000, 16777215));
        Information = new QWidget();
        Information->setObjectName(QString::fromUtf8("Information"));
        verticalLayout = new QVBoxLayout(Information);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        gridLayout = new QGridLayout();
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        gridLayout->setSizeConstraint(QLayout::SetDefaultConstraint);
        gridLayout->setContentsMargins(-1, 0, -1, -1);
        ElementNumber = new QPlainTextEdit(Information);
        ElementNumber->setObjectName(QString::fromUtf8("ElementNumber"));
        ElementNumber->setEnabled(true);
        ElementNumber->setMaximumSize(QSize(16777215, 28));
        ElementNumber->setReadOnly(true);

        gridLayout->addWidget(ElementNumber, 6, 2, 1, 1);

        Dims = new QPlainTextEdit(Information);
        Dims->setObjectName(QString::fromUtf8("Dims"));
        Dims->setEnabled(true);
        Dims->setMaximumSize(QSize(16777215, 28));
        Dims->setReadOnly(true);

        gridLayout->addWidget(Dims, 4, 2, 1, 1);

        label_2 = new QLabel(Information);
        label_2->setObjectName(QString::fromUtf8("label_2"));

        gridLayout->addWidget(label_2, 1, 1, 1, 1, Qt::AlignRight);

        ResultList = new QListWidget(Information);
        ResultList->setObjectName(QString::fromUtf8("ResultList"));
        ResultList->setSelectionMode(QAbstractItemView::ExtendedSelection);

        gridLayout->addWidget(ResultList, 10, 2, 6, 1);

        line = new QFrame(Information);
        line->setObjectName(QString::fromUtf8("line"));
        line->setFrameShape(QFrame::HLine);
        line->setFrameShadow(QFrame::Sunken);

        gridLayout->addWidget(line, 2, 0, 1, 4);

        label_5 = new QLabel(Information);
        label_5->setObjectName(QString::fromUtf8("label_5"));

        gridLayout->addWidget(label_5, 5, 1, 1, 1, Qt::AlignRight);

        New = new QPushButton(Information);
        New->setObjectName(QString::fromUtf8("New"));
        QIcon icon1;
        icon1.addFile(QString::fromUtf8(":/icons/resources/document-new.svg"), QSize(), QIcon::Normal, QIcon::Off);
        New->setIcon(icon1);
        New->setIconSize(QSize(22, 22));

        gridLayout->addWidget(New, 0, 3, 1, 1);

        ElementSize = new QPlainTextEdit(Information);
        ElementSize->setObjectName(QString::fromUtf8("ElementSize"));
        ElementSize->setEnabled(true);
        ElementSize->setMaximumSize(QSize(16777215, 28));
        ElementSize->setReadOnly(true);

        gridLayout->addWidget(ElementSize, 5, 2, 1, 1);

        verticalSpacer = new QSpacerItem(20, 20, QSizePolicy::Minimum, QSizePolicy::Fixed);

        gridLayout->addItem(verticalSpacer, 13, 3, 1, 1);

        verticalSpacer_2 = new QSpacerItem(20, 20, QSizePolicy::Minimum, QSizePolicy::Fixed);

        gridLayout->addItem(verticalSpacer_2, 8, 0, 1, 1);

        Save = new QPushButton(Information);
        Save->setObjectName(QString::fromUtf8("Save"));
        Save->setEnabled(false);
        QIcon icon2;
        icon2.addFile(QString::fromUtf8(":/icons/resources/document-save.svg"), QSize(), QIcon::Normal, QIcon::Off);
        Save->setIcon(icon2);
        Save->setIconSize(QSize(22, 22));

        gridLayout->addWidget(Save, 1, 3, 1, 1);

        Clear = new QPushButton(Information);
        Clear->setObjectName(QString::fromUtf8("Clear"));
        QIcon icon3;
        icon3.addFile(QString::fromUtf8(":/icons/resources/edit-clear.svg"), QSize(), QIcon::Normal, QIcon::Off);
        Clear->setIcon(icon3);
        Clear->setIconSize(QSize(22, 22));

        gridLayout->addWidget(Clear, 11, 3, 1, 1);

        Name = new QPlainTextEdit(Information);
        Name->setObjectName(QString::fromUtf8("Name"));
        Name->setMaximumSize(QSize(16777215, 28));

        gridLayout->addWidget(Name, 0, 2, 1, 1);

        label_9 = new QLabel(Information);
        label_9->setObjectName(QString::fromUtf8("label_9"));
        label_9->setMaximumSize(QSize(16777215, 28));

        gridLayout->addWidget(label_9, 9, 0, 1, 2);

        InspectCell = new QPushButton(Information);
        InspectCell->setObjectName(QString::fromUtf8("InspectCell"));
        QIcon icon4;
        icon4.addFile(QString::fromUtf8(":/icons/resources/system-search.svg"), QSize(), QIcon::Normal, QIcon::Off);
        InspectCell->setIcon(icon4);
        InspectCell->setIconSize(QSize(22, 22));

        gridLayout->addWidget(InspectCell, 12, 3, 1, 1);

        Delete = new QPushButton(Information);
        Delete->setObjectName(QString::fromUtf8("Delete"));
        QIcon icon5;
        icon5.addFile(QString::fromUtf8(":/icons/resources/edit-delete.svg"), QSize(), QIcon::Normal, QIcon::Off);
        Delete->setIcon(icon5);
        Delete->setIconSize(QSize(22, 22));

        gridLayout->addWidget(Delete, 14, 3, 1, 1);

        SelectAll = new QPushButton(Information);
        SelectAll->setObjectName(QString::fromUtf8("SelectAll"));
        QIcon icon6;
        icon6.addFile(QString::fromUtf8(":/icons/resources/edit-select-all.svg"), QSize(), QIcon::Normal, QIcon::Off);
        SelectAll->setIcon(icon6);
        SelectAll->setIconSize(QSize(22, 22));

        gridLayout->addWidget(SelectAll, 10, 3, 1, 1);

        label_10 = new QLabel(Information);
        label_10->setObjectName(QString::fromUtf8("label_10"));

        gridLayout->addWidget(label_10, 10, 1, 1, 1);

        verticalSpacer_4 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Preferred);

        gridLayout->addItem(verticalSpacer_4, 16, 2, 1, 1);

        label = new QLabel(Information);
        label->setObjectName(QString::fromUtf8("label"));

        gridLayout->addWidget(label, 0, 1, 1, 1, Qt::AlignRight);

        STLFile = new QPlainTextEdit(Information);
        STLFile->setObjectName(QString::fromUtf8("STLFile"));
        STLFile->setEnabled(true);
        STLFile->setMaximumSize(QSize(16777215, 28));
        STLFile->setReadOnly(true);

        gridLayout->addWidget(STLFile, 7, 2, 1, 1);

        label_4 = new QLabel(Information);
        label_4->setObjectName(QString::fromUtf8("label_4"));

        gridLayout->addWidget(label_4, 4, 1, 1, 1, Qt::AlignRight);

        horizontalSpacer_3 = new QSpacerItem(40, 20, QSizePolicy::Fixed, QSizePolicy::Minimum);

        gridLayout->addItem(horizontalSpacer_3, 4, 0, 1, 1);

        label_6 = new QLabel(Information);
        label_6->setObjectName(QString::fromUtf8("label_6"));

        gridLayout->addWidget(label_6, 6, 1, 1, 1, Qt::AlignRight);

        Path = new QPlainTextEdit(Information);
        Path->setObjectName(QString::fromUtf8("Path"));
        Path->setEnabled(true);
        Path->setMaximumSize(QSize(16777215, 28));
        Path->setReadOnly(true);

        gridLayout->addWidget(Path, 1, 2, 1, 1);

        horizontalSpacer_2 = new QSpacerItem(40, 20, QSizePolicy::Fixed, QSizePolicy::Minimum);

        gridLayout->addItem(horizontalSpacer_2, 0, 0, 1, 1);

        verticalSpacer_3 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        gridLayout->addItem(verticalSpacer_3, 15, 3, 1, 1);

        label_3 = new QLabel(Information);
        label_3->setObjectName(QString::fromUtf8("label_3"));
        label_3->setMaximumSize(QSize(16777215, 28));

        gridLayout->addWidget(label_3, 3, 0, 1, 3);

        label_7 = new QLabel(Information);
        label_7->setObjectName(QString::fromUtf8("label_7"));

        gridLayout->addWidget(label_7, 7, 1, 1, 1, Qt::AlignRight);


        verticalLayout->addLayout(gridLayout);

        tabWidget->addTab(Information, QString());
        Models = new QWidget();
        Models->setObjectName(QString::fromUtf8("Models"));
        Models->setEnabled(false);
        gridLayout_2 = new QGridLayout(Models);
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        gridLayout_4 = new QGridLayout();
        gridLayout_4->setObjectName(QString::fromUtf8("gridLayout_4"));
        gx = new QPlainTextEdit(Models);
        gx->setObjectName(QString::fromUtf8("gx"));
        gx->setMaximumSize(QSize(16777215, 28));

        gridLayout_4->addWidget(gx, 1, 3, 1, 1);

        calcRough = new QPushButton(Models);
        calcRough->setObjectName(QString::fromUtf8("calcRough"));
        QIcon icon7;
        icon7.addFile(QString::fromUtf8(":/icons/resources/accessories-calculator.svg"), QSize(), QIcon::Normal, QIcon::Off);
        calcRough->setIcon(icon7);

        gridLayout_4->addWidget(calcRough, 11, 12, 1, 1);

        rho = new QPlainTextEdit(Models);
        rho->setObjectName(QString::fromUtf8("rho"));
        rho->setMaximumSize(QSize(16777215, 28));

        gridLayout_4->addWidget(rho, 3, 10, 1, 1);

        horizontalSpacer_7 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_7, 1, 8, 1, 1);

        kEpsilon = new QRadioButton(Models);
        kEpsilon->setObjectName(QString::fromUtf8("kEpsilon"));
        kEpsilon->setChecked(true);

        gridLayout_4->addWidget(kEpsilon, 8, 1, 1, 1);

        length_scale = new QPlainTextEdit(Models);
        length_scale->setObjectName(QString::fromUtf8("length_scale"));
        length_scale->setMaximumSize(QSize(16777215, 28));

        gridLayout_4->addWidget(length_scale, 10, 7, 1, 1);

        horizontalSpacer_8 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_8, 1, 13, 1, 1);

        label_16 = new QLabel(Models);
        label_16->setObjectName(QString::fromUtf8("label_16"));

        gridLayout_4->addWidget(label_16, 1, 5, 1, 2);

        label_20 = new QLabel(Models);
        label_20->setObjectName(QString::fromUtf8("label_20"));
        label_20->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_4->addWidget(label_20, 11, 3, 1, 3);

        label_14 = new QLabel(Models);
        label_14->setObjectName(QString::fromUtf8("label_14"));
        label_14->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_4->addWidget(label_14, 3, 5, 1, 5);

        verticalSpacer_5 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        gridLayout_4->addItem(verticalSpacer_5, 12, 3, 1, 1);

        label_18 = new QLabel(Models);
        label_18->setObjectName(QString::fromUtf8("label_18"));
        QFont font;
        font.setBold(true);
        font.setWeight(75);
        label_18->setFont(font);

        gridLayout_4->addWidget(label_18, 6, 0, 1, 1);

        Laminar = new QRadioButton(Models);
        Laminar->setObjectName(QString::fromUtf8("Laminar"));

        gridLayout_4->addWidget(Laminar, 7, 1, 1, 1);

        gz = new QPlainTextEdit(Models);
        gz->setObjectName(QString::fromUtf8("gz"));
        gz->setMaximumSize(QSize(16777215, 28));

        gridLayout_4->addWidget(gz, 1, 10, 1, 1);

        verticalSpacer_6 = new QSpacerItem(20, 10, QSizePolicy::Minimum, QSizePolicy::Fixed);

        gridLayout_4->addItem(verticalSpacer_6, 2, 1, 1, 1);

        label_11 = new QLabel(Models);
        label_11->setObjectName(QString::fromUtf8("label_11"));
        label_11->setMaximumSize(QSize(16777215, 28));
        label_11->setFont(font);

        gridLayout_4->addWidget(label_11, 0, 0, 1, 1);

        label_19 = new QLabel(Models);
        label_19->setObjectName(QString::fromUtf8("label_19"));
        label_19->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_4->addWidget(label_19, 10, 3, 1, 3);

        horizontalSpacer_4 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_4, 1, 4, 1, 1);

        line_2 = new QFrame(Models);
        line_2->setObjectName(QString::fromUtf8("line_2"));
        line_2->setFrameShape(QFrame::HLine);
        line_2->setFrameShadow(QFrame::Sunken);

        gridLayout_4->addWidget(line_2, 5, 0, 1, 13);

        water20C = new QPushButton(Models);
        water20C->setObjectName(QString::fromUtf8("water20C"));
        QIcon icon8;
        icon8.addFile(QString::fromUtf8(":/icons/resources/water.png"), QSize(), QIcon::Normal, QIcon::Off);
        water20C->setIcon(icon8);

        gridLayout_4->addWidget(water20C, 3, 12, 1, 1);

        defaultLength = new QPushButton(Models);
        defaultLength->setObjectName(QString::fromUtf8("defaultLength"));
        defaultLength->setMinimumSize(QSize(0, 0));
        QIcon icon9;
        icon9.addFile(QString::fromUtf8(":/icons/resources/view-refresh.svg"), QSize(), QIcon::Normal, QIcon::Off);
        defaultLength->setIcon(icon9);

        gridLayout_4->addWidget(defaultLength, 10, 12, 1, 1);

        label_15 = new QLabel(Models);
        label_15->setObjectName(QString::fromUtf8("label_15"));

        gridLayout_4->addWidget(label_15, 1, 2, 1, 1);

        label_12 = new QLabel(Models);
        label_12->setObjectName(QString::fromUtf8("label_12"));

        gridLayout_4->addWidget(label_12, 1, 1, 1, 1);

        label_13 = new QLabel(Models);
        label_13->setObjectName(QString::fromUtf8("label_13"));
        label_13->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_4->addWidget(label_13, 3, 0, 1, 3);

        verticalSpacer_7 = new QSpacerItem(20, 10, QSizePolicy::Minimum, QSizePolicy::Fixed);

        gridLayout_4->addItem(verticalSpacer_7, 4, 1, 1, 1);

        gy = new QPlainTextEdit(Models);
        gy->setObjectName(QString::fromUtf8("gy"));
        gy->setMaximumSize(QSize(16777215, 28));

        gridLayout_4->addWidget(gy, 1, 7, 1, 1);

        nu = new QPlainTextEdit(Models);
        nu->setObjectName(QString::fromUtf8("nu"));
        nu->setMaximumSize(QSize(16777215, 28));

        gridLayout_4->addWidget(nu, 3, 3, 1, 1);

        horizontalSpacer_5 = new QSpacerItem(40, 20, QSizePolicy::Minimum, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_5, 1, 11, 1, 1);

        label_17 = new QLabel(Models);
        label_17->setObjectName(QString::fromUtf8("label_17"));

        gridLayout_4->addWidget(label_17, 1, 9, 1, 1);

        rough = new QPlainTextEdit(Models);
        rough->setObjectName(QString::fromUtf8("rough"));
        rough->setMaximumSize(QSize(16777215, 28));

        gridLayout_4->addWidget(rough, 11, 7, 1, 1);

        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Fixed, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer, 1, 0, 1, 1);

        earthGravity = new QPushButton(Models);
        earthGravity->setObjectName(QString::fromUtf8("earthGravity"));
        QIcon icon10;
        icon10.addFile(QString::fromUtf8(":/icons/resources/internet-web-browser.svg"), QSize(), QIcon::Normal, QIcon::Off);
        earthGravity->setIcon(icon10);

        gridLayout_4->addWidget(earthGravity, 1, 12, 1, 1);

        length = new QPlainTextEdit(Models);
        length->setObjectName(QString::fromUtf8("length"));
        length->setMaximumSize(QSize(16777215, 28));

        gridLayout_4->addWidget(length, 9, 7, 1, 1);

        label_8 = new QLabel(Models);
        label_8->setObjectName(QString::fromUtf8("label_8"));

        gridLayout_4->addWidget(label_8, 9, 3, 1, 3, Qt::AlignRight);

        domainLength = new QPushButton(Models);
        domainLength->setObjectName(QString::fromUtf8("domainLength"));
        QIcon icon11;
        icon11.addFile(QString::fromUtf8(":/icons/resources/view-fullscreen.svg"), QSize(), QIcon::Normal, QIcon::Off);
        domainLength->setIcon(icon11);

        gridLayout_4->addWidget(domainLength, 9, 12, 1, 1);


        gridLayout_2->addLayout(gridLayout_4, 1, 2, 1, 1);

        tabWidget->addTab(Models, QString());
        Mesh = new QWidget();
        Mesh->setObjectName(QString::fromUtf8("Mesh"));
        Mesh->setEnabled(false);
        gridLayout_6 = new QGridLayout(Mesh);
        gridLayout_6->setObjectName(QString::fromUtf8("gridLayout_6"));
        gridLayout_5 = new QGridLayout();
        gridLayout_5->setObjectName(QString::fromUtf8("gridLayout_5"));
        widget = new QWidget(Mesh);
        widget->setObjectName(QString::fromUtf8("widget"));
        widget->setMinimumSize(QSize(500, 0));
        gridLayout_7 = new QGridLayout(widget);
        gridLayout_7->setContentsMargins(0, 0, 0, 0);
        gridLayout_7->setObjectName(QString::fromUtf8("gridLayout_7"));
        MeshUndo = new QPushButton(widget);
        MeshUndo->setObjectName(QString::fromUtf8("MeshUndo"));
        MeshUndo->setMaximumSize(QSize(160, 16777215));
        QIcon icon12;
        icon12.addFile(QString::fromUtf8(":/icons/resources/edit-undo.svg"), QSize(), QIcon::Normal, QIcon::Off);
        MeshUndo->setIcon(icon12);
        MeshUndo->setIconSize(QSize(22, 22));

        gridLayout_7->addWidget(MeshUndo, 0, 1, 1, 1, Qt::AlignLeft);

        VTKMesh = new QVTKWidget(widget);
        VTKMesh->setObjectName(QString::fromUtf8("VTKMesh"));
        sizePolicy.setHeightForWidth(VTKMesh->sizePolicy().hasHeightForWidth());
        VTKMesh->setSizePolicy(sizePolicy);
        VTKMesh->setMinimumSize(QSize(500, 0));

        gridLayout_7->addWidget(VTKMesh, 1, 1, 1, 4);

        MeshParameters = new QTreeWidget(widget);
        QTreeWidgetItem *__qtreewidgetitem = new QTreeWidgetItem(MeshParameters);
        __qtreewidgetitem->setFont(0, font);
        QTreeWidgetItem *__qtreewidgetitem1 = new QTreeWidgetItem(__qtreewidgetitem);
        __qtreewidgetitem1->setFlags(Qt::ItemIsSelectable|Qt::ItemIsEditable|Qt::ItemIsDragEnabled|Qt::ItemIsUserCheckable|Qt::ItemIsEnabled);
        QTreeWidgetItem *__qtreewidgetitem2 = new QTreeWidgetItem(__qtreewidgetitem);
        __qtreewidgetitem2->setFlags(Qt::ItemIsSelectable|Qt::ItemIsEditable|Qt::ItemIsDragEnabled|Qt::ItemIsUserCheckable|Qt::ItemIsEnabled);
        QTreeWidgetItem *__qtreewidgetitem3 = new QTreeWidgetItem(__qtreewidgetitem);
        __qtreewidgetitem3->setFlags(Qt::ItemIsSelectable|Qt::ItemIsEditable|Qt::ItemIsDragEnabled|Qt::ItemIsUserCheckable|Qt::ItemIsEnabled);
        QTreeWidgetItem *__qtreewidgetitem4 = new QTreeWidgetItem(MeshParameters);
        __qtreewidgetitem4->setFont(0, font);
        QTreeWidgetItem *__qtreewidgetitem5 = new QTreeWidgetItem(__qtreewidgetitem4);
        __qtreewidgetitem5->setFlags(Qt::ItemIsSelectable|Qt::ItemIsEditable|Qt::ItemIsDragEnabled|Qt::ItemIsUserCheckable|Qt::ItemIsEnabled);
        QTreeWidgetItem *__qtreewidgetitem6 = new QTreeWidgetItem(__qtreewidgetitem4);
        __qtreewidgetitem6->setFlags(Qt::ItemIsSelectable|Qt::ItemIsEditable|Qt::ItemIsDragEnabled|Qt::ItemIsUserCheckable|Qt::ItemIsEnabled);
        QTreeWidgetItem *__qtreewidgetitem7 = new QTreeWidgetItem(__qtreewidgetitem4);
        __qtreewidgetitem7->setFlags(Qt::ItemIsSelectable|Qt::ItemIsEditable|Qt::ItemIsDragEnabled|Qt::ItemIsUserCheckable|Qt::ItemIsEnabled);
        QTreeWidgetItem *__qtreewidgetitem8 = new QTreeWidgetItem(MeshParameters);
        __qtreewidgetitem8->setFont(0, font);
        QTreeWidgetItem *__qtreewidgetitem9 = new QTreeWidgetItem(__qtreewidgetitem8);
        __qtreewidgetitem9->setFlags(Qt::ItemIsSelectable|Qt::ItemIsEditable|Qt::ItemIsDragEnabled|Qt::ItemIsUserCheckable|Qt::ItemIsEnabled);
        QTreeWidgetItem *__qtreewidgetitem10 = new QTreeWidgetItem(__qtreewidgetitem8);
        __qtreewidgetitem10->setFlags(Qt::ItemIsSelectable|Qt::ItemIsEditable|Qt::ItemIsDragEnabled|Qt::ItemIsUserCheckable|Qt::ItemIsEnabled);
        QTreeWidgetItem *__qtreewidgetitem11 = new QTreeWidgetItem(__qtreewidgetitem8);
        __qtreewidgetitem11->setFlags(Qt::ItemIsSelectable|Qt::ItemIsEditable|Qt::ItemIsDragEnabled|Qt::ItemIsUserCheckable|Qt::ItemIsEnabled);
        MeshParameters->setObjectName(QString::fromUtf8("MeshParameters"));
        MeshParameters->setMinimumSize(QSize(200, 0));
        MeshParameters->setMaximumSize(QSize(400, 16777215));
        MeshParameters->setEditTriggers(QAbstractItemView::NoEditTriggers);
        MeshParameters->setAlternatingRowColors(true);
        MeshParameters->setSelectionBehavior(QAbstractItemView::SelectItems);
        MeshParameters->setRootIsDecorated(true);
        MeshParameters->setItemsExpandable(true);
        MeshParameters->setAllColumnsShowFocus(false);
        MeshParameters->setHeaderHidden(false);
        MeshParameters->header()->setCascadingSectionResizes(false);
        MeshParameters->header()->setDefaultSectionSize(200);
        MeshParameters->header()->setHighlightSections(false);
        MeshParameters->header()->setProperty("showSortIndicator", QVariant(false));

        gridLayout_7->addWidget(MeshParameters, 0, 0, 2, 1);

        MeshUpdate = new QPushButton(widget);
        MeshUpdate->setObjectName(QString::fromUtf8("MeshUpdate"));
        MeshUpdate->setMaximumSize(QSize(160, 16777215));
        MeshUpdate->setIcon(icon9);
        MeshUpdate->setIconSize(QSize(22, 22));
        MeshUpdate->setCheckable(false);
        MeshUpdate->setFlat(false);

        gridLayout_7->addWidget(MeshUpdate, 0, 2, 1, 1);


        gridLayout_5->addWidget(widget, 0, 0, 1, 1);


        gridLayout_6->addLayout(gridLayout_5, 1, 0, 1, 1);

        tabWidget->addTab(Mesh, QString());
        Geometry = new QWidget();
        Geometry->setObjectName(QString::fromUtf8("Geometry"));
        Geometry->setEnabled(false);
        gridLayout_10 = new QGridLayout(Geometry);
        gridLayout_10->setObjectName(QString::fromUtf8("gridLayout_10"));
        gridLayout_9 = new QGridLayout();
        gridLayout_9->setObjectName(QString::fromUtf8("gridLayout_9"));
        label_22 = new QLabel(Geometry);
        label_22->setObjectName(QString::fromUtf8("label_22"));

        gridLayout_9->addWidget(label_22, 1, 1, 1, 1, Qt::AlignRight);

        STLFile_2 = new QPlainTextEdit(Geometry);
        STLFile_2->setObjectName(QString::fromUtf8("STLFile_2"));
        STLFile_2->setMaximumSize(QSize(16777215, 28));
        STLFile_2->setReadOnly(true);

        gridLayout_9->addWidget(STLFile_2, 0, 3, 1, 5);

        VTKGeometry = new QVTKWidget(Geometry);
        VTKGeometry->setObjectName(QString::fromUtf8("VTKGeometry"));

        gridLayout_9->addWidget(VTKGeometry, 3, 1, 1, 9);

        label_25 = new QLabel(Geometry);
        label_25->setObjectName(QString::fromUtf8("label_25"));
        label_25->setMinimumSize(QSize(20, 0));
        label_25->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_9->addWidget(label_25, 1, 6, 1, 1, Qt::AlignRight);

        label_24 = new QLabel(Geometry);
        label_24->setObjectName(QString::fromUtf8("label_24"));
        label_24->setMinimumSize(QSize(20, 0));
        label_24->setLayoutDirection(Qt::LeftToRight);
        label_24->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_9->addWidget(label_24, 1, 4, 1, 1, Qt::AlignRight);

        label_23 = new QLabel(Geometry);
        label_23->setObjectName(QString::fromUtf8("label_23"));
        label_23->setMinimumSize(QSize(20, 0));
        label_23->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_9->addWidget(label_23, 1, 2, 1, 1, Qt::AlignRight);

        inside_x = new QPlainTextEdit(Geometry);
        inside_x->setObjectName(QString::fromUtf8("inside_x"));
        inside_x->setMaximumSize(QSize(100, 28));

        gridLayout_9->addWidget(inside_x, 1, 3, 1, 1);

        label_21 = new QLabel(Geometry);
        label_21->setObjectName(QString::fromUtf8("label_21"));

        gridLayout_9->addWidget(label_21, 0, 1, 1, 1, Qt::AlignRight);

        inside_z = new QPlainTextEdit(Geometry);
        inside_z->setObjectName(QString::fromUtf8("inside_z"));
        inside_z->setMaximumSize(QSize(100, 28));

        gridLayout_9->addWidget(inside_z, 1, 7, 1, 1);

        inside_y = new QPlainTextEdit(Geometry);
        inside_y->setObjectName(QString::fromUtf8("inside_y"));
        inside_y->setMaximumSize(QSize(100, 28));

        gridLayout_9->addWidget(inside_y, 1, 5, 1, 1);

        horizontalSpacer_6 = new QSpacerItem(60, 20, QSizePolicy::Fixed, QSizePolicy::Minimum);

        gridLayout_9->addItem(horizontalSpacer_6, 1, 8, 1, 1);

        verticalSpacer_8 = new QSpacerItem(20, 10, QSizePolicy::Minimum, QSizePolicy::Maximum);

        gridLayout_9->addItem(verticalSpacer_8, 2, 5, 1, 1);

        STLRender = new QPushButton(Geometry);
        STLRender->setObjectName(QString::fromUtf8("STLRender"));
        STLRender->setMinimumSize(QSize(100, 0));
        QIcon icon13;
        icon13.addFile(QString::fromUtf8(":/icons/resources/media-playback-start.svg"), QSize(), QIcon::Normal, QIcon::Off);
        STLRender->setIcon(icon13);
        STLRender->setIconSize(QSize(22, 22));

        gridLayout_9->addWidget(STLRender, 1, 9, 1, 1);

        horizontalSpacer_9 = new QSpacerItem(40, 20, QSizePolicy::Fixed, QSizePolicy::Minimum);

        gridLayout_9->addItem(horizontalSpacer_9, 0, 2, 1, 1);

        STLOpen = new QPushButton(Geometry);
        STLOpen->setObjectName(QString::fromUtf8("STLOpen"));
        STLOpen->setMinimumSize(QSize(100, 0));
        QIcon icon14;
        icon14.addFile(QString::fromUtf8(":/icons/resources/document-open.svg"), QSize(), QIcon::Normal, QIcon::Off);
        STLOpen->setIcon(icon14);
        STLOpen->setIconSize(QSize(22, 22));

        gridLayout_9->addWidget(STLOpen, 0, 9, 1, 1);


        gridLayout_10->addLayout(gridLayout_9, 0, 0, 1, 1);

        tabWidget->addTab(Geometry, QString());
        Boundaries = new QWidget();
        Boundaries->setObjectName(QString::fromUtf8("Boundaries"));
        Boundaries->setEnabled(false);
        gridLayout_8 = new QGridLayout(Boundaries);
        gridLayout_8->setObjectName(QString::fromUtf8("gridLayout_8"));
        gridLayout_12 = new QGridLayout();
        gridLayout_12->setObjectName(QString::fromUtf8("gridLayout_12"));
        BoundaryTree = new QTreeWidget(Boundaries);
        QTreeWidgetItem *__qtreewidgetitem12 = new QTreeWidgetItem(BoundaryTree);
        __qtreewidgetitem12->setFlags(Qt::ItemIsDragEnabled|Qt::ItemIsEnabled);
        __qtreewidgetitem12->setFont(0, font);
        new QTreeWidgetItem(__qtreewidgetitem12);
        new QTreeWidgetItem(__qtreewidgetitem12);
        new QTreeWidgetItem(__qtreewidgetitem12);
        new QTreeWidgetItem(__qtreewidgetitem12);
        new QTreeWidgetItem(__qtreewidgetitem12);
        new QTreeWidgetItem(__qtreewidgetitem12);
        QTreeWidgetItem *__qtreewidgetitem13 = new QTreeWidgetItem(BoundaryTree);
        __qtreewidgetitem13->setFlags(Qt::ItemIsDragEnabled|Qt::ItemIsEnabled);
        __qtreewidgetitem13->setFont(0, font);
        new QTreeWidgetItem(__qtreewidgetitem13);
        new QTreeWidgetItem(__qtreewidgetitem13);
        new QTreeWidgetItem(__qtreewidgetitem13);
        new QTreeWidgetItem(__qtreewidgetitem13);
        new QTreeWidgetItem(__qtreewidgetitem13);
        new QTreeWidgetItem(__qtreewidgetitem13);
        BoundaryTree->setObjectName(QString::fromUtf8("BoundaryTree"));
        sizePolicy.setHeightForWidth(BoundaryTree->sizePolicy().hasHeightForWidth());
        BoundaryTree->setSizePolicy(sizePolicy);
        BoundaryTree->setMinimumSize(QSize(320, 0));
        BoundaryTree->setAlternatingRowColors(true);
        BoundaryTree->header()->setDefaultSectionSize(165);

        gridLayout_12->addWidget(BoundaryTree, 1, 1, 1, 1);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        horizontalLayout->setSizeConstraint(QLayout::SetFixedSize);
        EditBoundary = new QPushButton(Boundaries);
        EditBoundary->setObjectName(QString::fromUtf8("EditBoundary"));
        EditBoundary->setMinimumSize(QSize(0, 0));
        EditBoundary->setMaximumSize(QSize(140, 16777215));
        QIcon icon15;
        icon15.addFile(QString::fromUtf8(":/icons/resources/edit-find-replace.svg"), QSize(), QIcon::Normal, QIcon::Off);
        EditBoundary->setIcon(icon15);
        EditBoundary->setIconSize(QSize(22, 22));

        horizontalLayout->addWidget(EditBoundary);

        AddSpecialBoundary = new QPushButton(Boundaries);
        AddSpecialBoundary->setObjectName(QString::fromUtf8("AddSpecialBoundary"));
        AddSpecialBoundary->setMaximumSize(QSize(190, 16777215));
        QIcon icon16;
        icon16.addFile(QString::fromUtf8(":/icons/resources/list-add.svg"), QSize(), QIcon::Normal, QIcon::Off);
        AddSpecialBoundary->setIcon(icon16);
        AddSpecialBoundary->setIconSize(QSize(22, 22));

        horizontalLayout->addWidget(AddSpecialBoundary);

        RemoveSpecialBoundary = new QPushButton(Boundaries);
        RemoveSpecialBoundary->setObjectName(QString::fromUtf8("RemoveSpecialBoundary"));
        RemoveSpecialBoundary->setMaximumSize(QSize(220, 16777215));
        RemoveSpecialBoundary->setIcon(icon5);
        RemoveSpecialBoundary->setIconSize(QSize(22, 22));

        horizontalLayout->addWidget(RemoveSpecialBoundary);


        gridLayout_12->addLayout(horizontalLayout, 0, 1, 1, 1);

        VTKBoundaries = new QVTKWidget(Boundaries);
        VTKBoundaries->setObjectName(QString::fromUtf8("VTKBoundaries"));
        sizePolicy.setHeightForWidth(VTKBoundaries->sizePolicy().hasHeightForWidth());
        VTKBoundaries->setSizePolicy(sizePolicy);
        VTKBoundaries->setMinimumSize(QSize(400, 310));

        gridLayout_12->addWidget(VTKBoundaries, 3, 1, 1, 1);

        horizontalLayout_6 = new QHBoxLayout();
        horizontalLayout_6->setObjectName(QString::fromUtf8("horizontalLayout_6"));
        horizontalSpacer_24 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_6->addItem(horizontalSpacer_24);

        hideBoundaryMesh = new QCheckBox(Boundaries);
        hideBoundaryMesh->setObjectName(QString::fromUtf8("hideBoundaryMesh"));

        horizontalLayout_6->addWidget(hideBoundaryMesh);

        horizontalSpacer_26 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_6->addItem(horizontalSpacer_26);

        hideBoundaryDomain = new QCheckBox(Boundaries);
        hideBoundaryDomain->setObjectName(QString::fromUtf8("hideBoundaryDomain"));

        horizontalLayout_6->addWidget(hideBoundaryDomain);

        horizontalSpacer_25 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_6->addItem(horizontalSpacer_25);


        gridLayout_12->addLayout(horizontalLayout_6, 2, 1, 1, 1);


        gridLayout_8->addLayout(gridLayout_12, 0, 0, 1, 1);

        tabWidget->addTab(Boundaries, QString());
        Baffles = new QWidget();
        Baffles->setObjectName(QString::fromUtf8("Baffles"));
        Baffles->setEnabled(false);
        Baffles->setAutoFillBackground(false);
        gridLayout_20 = new QGridLayout(Baffles);
        gridLayout_20->setObjectName(QString::fromUtf8("gridLayout_20"));
        gridLayout_20->setContentsMargins(12, -1, 12, 12);
        gridLayout_19 = new QGridLayout();
        gridLayout_19->setObjectName(QString::fromUtf8("gridLayout_19"));
        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        horizontalLayout_2->setSizeConstraint(QLayout::SetFixedSize);
        EditBaffle = new QPushButton(Baffles);
        EditBaffle->setObjectName(QString::fromUtf8("EditBaffle"));
        EditBaffle->setMinimumSize(QSize(0, 0));
        EditBaffle->setMaximumSize(QSize(140, 16777215));
        EditBaffle->setIcon(icon15);
        EditBaffle->setIconSize(QSize(22, 22));

        horizontalLayout_2->addWidget(EditBaffle);

        AddBaffle = new QPushButton(Baffles);
        AddBaffle->setObjectName(QString::fromUtf8("AddBaffle"));
        AddBaffle->setMaximumSize(QSize(190, 16777215));
        AddBaffle->setIcon(icon16);
        AddBaffle->setIconSize(QSize(22, 22));

        horizontalLayout_2->addWidget(AddBaffle);

        RemoveBaffle = new QPushButton(Baffles);
        RemoveBaffle->setObjectName(QString::fromUtf8("RemoveBaffle"));
        RemoveBaffle->setMaximumSize(QSize(220, 16777215));
        RemoveBaffle->setIcon(icon5);
        RemoveBaffle->setIconSize(QSize(22, 22));

        horizontalLayout_2->addWidget(RemoveBaffle);


        gridLayout_19->addLayout(horizontalLayout_2, 0, 0, 1, 1);

        VTKBaffles = new QVTKWidget(Baffles);
        VTKBaffles->setObjectName(QString::fromUtf8("VTKBaffles"));
        sizePolicy.setHeightForWidth(VTKBaffles->sizePolicy().hasHeightForWidth());
        VTKBaffles->setSizePolicy(sizePolicy);
        VTKBaffles->setMinimumSize(QSize(400, 350));

        gridLayout_19->addWidget(VTKBaffles, 7, 0, 1, 1);

        BaffleTree = new QTreeWidget(Baffles);
        QTreeWidgetItem *__qtreewidgetitem14 = new QTreeWidgetItem(BaffleTree);
        __qtreewidgetitem14->setFlags(Qt::ItemIsDragEnabled|Qt::ItemIsEnabled);
        __qtreewidgetitem14->setFont(0, font);
        new QTreeWidgetItem(__qtreewidgetitem14);
        new QTreeWidgetItem(__qtreewidgetitem14);
        new QTreeWidgetItem(__qtreewidgetitem14);
        BaffleTree->setObjectName(QString::fromUtf8("BaffleTree"));
        sizePolicy.setHeightForWidth(BaffleTree->sizePolicy().hasHeightForWidth());
        BaffleTree->setSizePolicy(sizePolicy);
        BaffleTree->setMinimumSize(QSize(320, 0));
        BaffleTree->setAlternatingRowColors(true);
        BaffleTree->header()->setDefaultSectionSize(165);

        gridLayout_19->addWidget(BaffleTree, 1, 0, 1, 1);

        horizontalLayout_5 = new QHBoxLayout();
        horizontalLayout_5->setObjectName(QString::fromUtf8("horizontalLayout_5"));
        horizontalSpacer_21 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_5->addItem(horizontalSpacer_21);

        hideBaffleMesh = new QCheckBox(Baffles);
        hideBaffleMesh->setObjectName(QString::fromUtf8("hideBaffleMesh"));

        horizontalLayout_5->addWidget(hideBaffleMesh);

        horizontalSpacer_23 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_5->addItem(horizontalSpacer_23);

        hideBaffleDomain = new QCheckBox(Baffles);
        hideBaffleDomain->setObjectName(QString::fromUtf8("hideBaffleDomain"));

        horizontalLayout_5->addWidget(hideBaffleDomain);

        horizontalSpacer_22 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_5->addItem(horizontalSpacer_22);


        gridLayout_19->addLayout(horizontalLayout_5, 4, 0, 1, 1);


        gridLayout_20->addLayout(gridLayout_19, 0, 0, 1, 1);

        tabWidget->addTab(Baffles, QString());
        InitialConditions = new QWidget();
        InitialConditions->setObjectName(QString::fromUtf8("InitialConditions"));
        InitialConditions->setEnabled(false);
        gridLayout_14 = new QGridLayout(InitialConditions);
        gridLayout_14->setObjectName(QString::fromUtf8("gridLayout_14"));
        gridLayout_13 = new QGridLayout();
        gridLayout_13->setObjectName(QString::fromUtf8("gridLayout_13"));
        verticalSpacer_12 = new QSpacerItem(20, 20, QSizePolicy::Minimum, QSizePolicy::Fixed);

        gridLayout_13->addItem(verticalSpacer_12, 12, 2, 1, 1);

        label_26 = new QLabel(InitialConditions);
        label_26->setObjectName(QString::fromUtf8("label_26"));
        label_26->setFont(font);

        gridLayout_13->addWidget(label_26, 0, 0, 1, 1);

        horizontalSpacer_10 = new QSpacerItem(40, 20, QSizePolicy::Fixed, QSizePolicy::Minimum);

        gridLayout_13->addItem(horizontalSpacer_10, 1, 0, 1, 1);

        line_5 = new QFrame(InitialConditions);
        line_5->setObjectName(QString::fromUtf8("line_5"));
        line_5->setFrameShape(QFrame::HLine);
        line_5->setFrameShadow(QFrame::Sunken);

        gridLayout_13->addWidget(line_5, 13, 0, 1, 10);

        label_30 = new QLabel(InitialConditions);
        label_30->setObjectName(QString::fromUtf8("label_30"));

        gridLayout_13->addWidget(label_30, 1, 4, 1, 1, Qt::AlignRight);

        horizontalSpacer_12 = new QSpacerItem(40, 20, QSizePolicy::Fixed, QSizePolicy::Minimum);

        gridLayout_13->addItem(horizontalSpacer_12, 1, 3, 1, 1);

        label_29 = new QLabel(InitialConditions);
        label_29->setObjectName(QString::fromUtf8("label_29"));

        gridLayout_13->addWidget(label_29, 1, 1, 1, 1, Qt::AlignRight);

        line_4 = new QFrame(InitialConditions);
        line_4->setObjectName(QString::fromUtf8("line_4"));
        line_4->setFrameShape(QFrame::HLine);
        line_4->setFrameShadow(QFrame::Sunken);

        gridLayout_13->addWidget(line_4, 9, 0, 1, 10);

        verticalSpacer_11 = new QSpacerItem(20, 20, QSizePolicy::Minimum, QSizePolicy::Fixed);

        gridLayout_13->addItem(verticalSpacer_11, 6, 2, 1, 1);

        initialU = new QPlainTextEdit(InitialConditions);
        initialU->setObjectName(QString::fromUtf8("initialU"));
        initialU->setMaximumSize(QSize(16777215, 28));

        gridLayout_13->addWidget(initialU, 1, 2, 1, 1);

        horizontalSpacer_11 = new QSpacerItem(40, 20, QSizePolicy::Fixed, QSizePolicy::Minimum);

        gridLayout_13->addItem(horizontalSpacer_11, 1, 9, 1, 1);

        initialW = new QPlainTextEdit(InitialConditions);
        initialW->setObjectName(QString::fromUtf8("initialW"));
        initialW->setMaximumSize(QSize(16777215, 28));

        gridLayout_13->addWidget(initialW, 1, 8, 1, 1);

        label_28 = new QLabel(InitialConditions);
        label_28->setObjectName(QString::fromUtf8("label_28"));
        label_28->setFont(font);

        gridLayout_13->addWidget(label_28, 10, 0, 1, 1);

        height = new QPlainTextEdit(InitialConditions);
        height->setObjectName(QString::fromUtf8("height"));
        height->setMaximumSize(QSize(16777215, 28));

        gridLayout_13->addWidget(height, 5, 2, 1, 1);

        verticalSpacer_10 = new QSpacerItem(20, 20, QSizePolicy::Minimum, QSizePolicy::Fixed);

        gridLayout_13->addItem(verticalSpacer_10, 2, 1, 1, 1);

        label_32 = new QLabel(InitialConditions);
        label_32->setObjectName(QString::fromUtf8("label_32"));

        gridLayout_13->addWidget(label_32, 5, 0, 1, 2, Qt::AlignRight);

        line_3 = new QFrame(InitialConditions);
        line_3->setObjectName(QString::fromUtf8("line_3"));
        line_3->setFrameShape(QFrame::HLine);
        line_3->setFrameShadow(QFrame::Sunken);

        gridLayout_13->addWidget(line_3, 3, 0, 1, 10);

        label_27 = new QLabel(InitialConditions);
        label_27->setObjectName(QString::fromUtf8("label_27"));
        label_27->setFont(font);

        gridLayout_13->addWidget(label_27, 4, 0, 1, 1);

        initialV = new QPlainTextEdit(InitialConditions);
        initialV->setObjectName(QString::fromUtf8("initialV"));
        initialV->setMaximumSize(QSize(16777215, 28));

        gridLayout_13->addWidget(initialV, 1, 5, 1, 1);

        label_31 = new QLabel(InitialConditions);
        label_31->setObjectName(QString::fromUtf8("label_31"));

        gridLayout_13->addWidget(label_31, 1, 7, 1, 1, Qt::AlignRight);

        horizontalSpacer_13 = new QSpacerItem(40, 20, QSizePolicy::Fixed, QSizePolicy::Minimum);

        gridLayout_13->addItem(horizontalSpacer_13, 1, 6, 1, 1);

        hydrostatic = new QCheckBox(InitialConditions);
        hydrostatic->setObjectName(QString::fromUtf8("hydrostatic"));

        gridLayout_13->addWidget(hydrostatic, 11, 2, 1, 1);

        verticalSpacer_9 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        gridLayout_13->addItem(verticalSpacer_9, 17, 5, 1, 1);

        label_33 = new QLabel(InitialConditions);
        label_33->setObjectName(QString::fromUtf8("label_33"));
        label_33->setFont(font);

        gridLayout_13->addWidget(label_33, 14, 0, 1, 1);

        verticalSpacer_13 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        gridLayout_13->addItem(verticalSpacer_13, 16, 5, 1, 1);

        label_34 = new QLabel(InitialConditions);
        label_34->setObjectName(QString::fromUtf8("label_34"));

        gridLayout_13->addWidget(label_34, 15, 0, 1, 2, Qt::AlignRight);

        initialk = new QPlainTextEdit(InitialConditions);
        initialk->setObjectName(QString::fromUtf8("initialk"));
        initialk->setMaximumSize(QSize(16777215, 28));

        gridLayout_13->addWidget(initialk, 15, 2, 1, 1);

        label_43 = new QLabel(InitialConditions);
        label_43->setObjectName(QString::fromUtf8("label_43"));
        label_43->setLayoutDirection(Qt::LeftToRight);
        label_43->setAlignment(Qt::AlignRight|Qt::AlignTop|Qt::AlignTrailing);

        gridLayout_13->addWidget(label_43, 7, 0, 1, 2);

        verticalSpacer_24 = new QSpacerItem(20, 20, QSizePolicy::Minimum, QSizePolicy::Fixed);

        gridLayout_13->addItem(verticalSpacer_24, 8, 0, 1, 1);

        fillPoints = new QTableWidget(InitialConditions);
        if (fillPoints->columnCount() < 3)
            fillPoints->setColumnCount(3);
        QTableWidgetItem *__qtablewidgetitem = new QTableWidgetItem();
        fillPoints->setHorizontalHeaderItem(0, __qtablewidgetitem);
        QTableWidgetItem *__qtablewidgetitem1 = new QTableWidgetItem();
        fillPoints->setHorizontalHeaderItem(1, __qtablewidgetitem1);
        QTableWidgetItem *__qtablewidgetitem2 = new QTableWidgetItem();
        fillPoints->setHorizontalHeaderItem(2, __qtablewidgetitem2);
        if (fillPoints->rowCount() < 8)
            fillPoints->setRowCount(8);
        QTableWidgetItem *__qtablewidgetitem3 = new QTableWidgetItem();
        fillPoints->setVerticalHeaderItem(0, __qtablewidgetitem3);
        QTableWidgetItem *__qtablewidgetitem4 = new QTableWidgetItem();
        __qtablewidgetitem4->setFlags(Qt::ItemIsSelectable|Qt::ItemIsEditable|Qt::ItemIsEnabled);
        fillPoints->setItem(0, 0, __qtablewidgetitem4);
        fillPoints->setObjectName(QString::fromUtf8("fillPoints"));
        fillPoints->setRowCount(8);
        fillPoints->setColumnCount(3);
        fillPoints->horizontalHeader()->setDefaultSectionSize(75);
        fillPoints->horizontalHeader()->setHighlightSections(false);
        fillPoints->horizontalHeader()->setMinimumSectionSize(30);
        fillPoints->verticalHeader()->setVisible(false);
        fillPoints->verticalHeader()->setHighlightSections(false);

        gridLayout_13->addWidget(fillPoints, 8, 3, 1, 5);


        gridLayout_14->addLayout(gridLayout_13, 0, 0, 1, 1);

        tabWidget->addTab(InitialConditions, QString());
        Simulate = new QWidget();
        Simulate->setObjectName(QString::fromUtf8("Simulate"));
        Simulate->setEnabled(false);
        gridLayout_16 = new QGridLayout(Simulate);
        gridLayout_16->setObjectName(QString::fromUtf8("gridLayout_16"));
        gridLayout_15 = new QGridLayout();
        gridLayout_15->setObjectName(QString::fromUtf8("gridLayout_15"));
        label_37 = new QLabel(Simulate);
        label_37->setObjectName(QString::fromUtf8("label_37"));

        gridLayout_15->addWidget(label_37, 1, 4, 1, 1, Qt::AlignRight);

        horizontalSpacer_15 = new QSpacerItem(40, 20, QSizePolicy::Fixed, QSizePolicy::Minimum);

        gridLayout_15->addItem(horizontalSpacer_15, 1, 3, 1, 1);

        verticalSpacer_16 = new QSpacerItem(20, 10, QSizePolicy::Minimum, QSizePolicy::Fixed);

        gridLayout_15->addItem(verticalSpacer_16, 5, 2, 1, 1);

        t = new QComboBox(Simulate);
        t->setObjectName(QString::fromUtf8("t"));
        QSizePolicy sizePolicy2(QSizePolicy::Expanding, QSizePolicy::Fixed);
        sizePolicy2.setHorizontalStretch(0);
        sizePolicy2.setVerticalStretch(0);
        sizePolicy2.setHeightForWidth(t->sizePolicy().hasHeightForWidth());
        t->setSizePolicy(sizePolicy2);
        t->setMaximumSize(QSize(200, 28));

        gridLayout_15->addWidget(t, 1, 2, 1, 1);

        verticalSpacer_15 = new QSpacerItem(20, 10, QSizePolicy::Minimum, QSizePolicy::Fixed);

        gridLayout_15->addItem(verticalSpacer_15, 2, 2, 1, 1);

        horizontalSpacer_16 = new QSpacerItem(80, 20, QSizePolicy::MinimumExpanding, QSizePolicy::Minimum);

        gridLayout_15->addItem(horizontalSpacer_16, 1, 6, 1, 1);

        label_39 = new QLabel(Simulate);
        label_39->setObjectName(QString::fromUtf8("label_39"));

        gridLayout_15->addWidget(label_39, 3, 4, 1, 1, Qt::AlignRight);

        verticalSpacer_14 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        gridLayout_15->addItem(verticalSpacer_14, 12, 1, 1, 1);

        verticalSpacer_29 = new QSpacerItem(20, 10, QSizePolicy::Minimum, QSizePolicy::Fixed);

        gridLayout_15->addItem(verticalSpacer_29, 9, 2, 1, 1);

        line_6 = new QFrame(Simulate);
        line_6->setObjectName(QString::fromUtf8("line_6"));
        line_6->setFrameShape(QFrame::HLine);
        line_6->setFrameShadow(QFrame::Sunken);

        gridLayout_15->addWidget(line_6, 10, 0, 1, 7);

        autot = new QCheckBox(Simulate);
        autot->setObjectName(QString::fromUtf8("autot"));
        autot->setChecked(true);

        gridLayout_15->addWidget(autot, 4, 2, 1, 1);

        verticalSpacer_17 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Maximum);

        gridLayout_15->addItem(verticalSpacer_17, 11, 1, 1, 1);

        writet = new QPlainTextEdit(Simulate);
        writet->setObjectName(QString::fromUtf8("writet"));
        writet->setMaximumSize(QSize(200, 28));

        gridLayout_15->addWidget(writet, 3, 5, 1, 1);

        RunSimulation = new QPushButton(Simulate);
        RunSimulation->setObjectName(QString::fromUtf8("RunSimulation"));
        RunSimulation->setMinimumSize(QSize(120, 0));
        RunSimulation->setMaximumSize(QSize(120, 16777215));
        RunSimulation->setIcon(icon13);
        RunSimulation->setIconSize(QSize(24, 24));

        gridLayout_15->addWidget(RunSimulation, 12, 5, 1, 2, Qt::AlignRight|Qt::AlignTop);

        delt = new QPlainTextEdit(Simulate);
        delt->setObjectName(QString::fromUtf8("delt"));
        delt->setMaximumSize(QSize(200, 28));

        gridLayout_15->addWidget(delt, 3, 2, 1, 1);

        label_38 = new QLabel(Simulate);
        label_38->setObjectName(QString::fromUtf8("label_38"));

        gridLayout_15->addWidget(label_38, 3, 1, 1, 1, Qt::AlignRight);

        horizontalSpacer_14 = new QSpacerItem(40, 20, QSizePolicy::MinimumExpanding, QSizePolicy::Minimum);

        gridLayout_15->addItem(horizontalSpacer_14, 1, 0, 1, 1);

        label_47 = new QLabel(Simulate);
        label_47->setObjectName(QString::fromUtf8("label_47"));
        label_47->setFont(font);

        gridLayout_15->addWidget(label_47, 7, 0, 1, 1);

        label_36 = new QLabel(Simulate);
        label_36->setObjectName(QString::fromUtf8("label_36"));

        gridLayout_15->addWidget(label_36, 1, 1, 1, 1, Qt::AlignRight);

        endt = new QPlainTextEdit(Simulate);
        endt->setObjectName(QString::fromUtf8("endt"));
        sizePolicy2.setHeightForWidth(endt->sizePolicy().hasHeightForWidth());
        endt->setSizePolicy(sizePolicy2);
        endt->setMaximumSize(QSize(200, 28));

        gridLayout_15->addWidget(endt, 1, 5, 1, 1);

        label_35 = new QLabel(Simulate);
        label_35->setObjectName(QString::fromUtf8("label_35"));
        label_35->setFont(font);

        gridLayout_15->addWidget(label_35, 0, 0, 1, 1);

        line_7 = new QFrame(Simulate);
        line_7->setObjectName(QString::fromUtf8("line_7"));
        line_7->setFrameShape(QFrame::HLine);
        line_7->setFrameShadow(QFrame::Sunken);

        gridLayout_15->addWidget(line_7, 6, 0, 1, 7);

        label_49 = new QLabel(Simulate);
        label_49->setObjectName(QString::fromUtf8("label_49"));

        gridLayout_15->addWidget(label_49, 8, 1, 1, 1);

        processes = new QComboBox(Simulate);
        processes->setObjectName(QString::fromUtf8("processes"));

        gridLayout_15->addWidget(processes, 8, 2, 1, 1);


        gridLayout_16->addLayout(gridLayout_15, 0, 0, 1, 1);

        tabWidget->addTab(Simulate, QString());
        Visualize = new QWidget();
        Visualize->setObjectName(QString::fromUtf8("Visualize"));
        Visualize->setEnabled(false);
        gridLayout_18 = new QGridLayout(Visualize);
        gridLayout_18->setObjectName(QString::fromUtf8("gridLayout_18"));
        gridLayout_17 = new QGridLayout();
        gridLayout_17->setObjectName(QString::fromUtf8("gridLayout_17"));
        horizontalLayout_3 = new QHBoxLayout();
        horizontalLayout_3->setObjectName(QString::fromUtf8("horizontalLayout_3"));
        horizontalSpacer_17 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_3->addItem(horizontalSpacer_17);

        updateRange = new QCheckBox(Visualize);
        updateRange->setObjectName(QString::fromUtf8("updateRange"));

        horizontalLayout_3->addWidget(updateRange);

        from = new QPlainTextEdit(Visualize);
        from->setObjectName(QString::fromUtf8("from"));
        QSizePolicy sizePolicy3(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy3.setHorizontalStretch(0);
        sizePolicy3.setVerticalStretch(0);
        sizePolicy3.setHeightForWidth(from->sizePolicy().hasHeightForWidth());
        from->setSizePolicy(sizePolicy3);
        from->setMaximumSize(QSize(100, 28));

        horizontalLayout_3->addWidget(from);

        label_48 = new QLabel(Visualize);
        label_48->setObjectName(QString::fromUtf8("label_48"));
        label_48->setMaximumSize(QSize(17, 16777215));

        horizontalLayout_3->addWidget(label_48, 0, Qt::AlignHCenter);

        to = new QPlainTextEdit(Visualize);
        to->setObjectName(QString::fromUtf8("to"));
        QSizePolicy sizePolicy4(QSizePolicy::Minimum, QSizePolicy::Fixed);
        sizePolicy4.setHorizontalStretch(0);
        sizePolicy4.setVerticalStretch(0);
        sizePolicy4.setHeightForWidth(to->sizePolicy().hasHeightForWidth());
        to->setSizePolicy(sizePolicy4);
        to->setMaximumSize(QSize(100, 28));

        horizontalLayout_3->addWidget(to, 0, Qt::AlignLeft);

        horizontalSpacer_18 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_3->addItem(horizontalSpacer_18);


        gridLayout_17->addLayout(horizontalLayout_3, 12, 2, 1, 3);

        verticalSpacer_19 = new QSpacerItem(20, 5, QSizePolicy::Minimum, QSizePolicy::Fixed);

        gridLayout_17->addItem(verticalSpacer_19, 2, 2, 1, 1);

        verticalSpacer_21 = new QSpacerItem(20, 10, QSizePolicy::Minimum, QSizePolicy::Fixed);

        gridLayout_17->addItem(verticalSpacer_21, 7, 3, 1, 1);

        label_41 = new QLabel(Visualize);
        label_41->setObjectName(QString::fromUtf8("label_41"));
        QSizePolicy sizePolicy5(QSizePolicy::Preferred, QSizePolicy::Fixed);
        sizePolicy5.setHorizontalStretch(0);
        sizePolicy5.setVerticalStretch(0);
        sizePolicy5.setHeightForWidth(label_41->sizePolicy().hasHeightForWidth());
        label_41->setSizePolicy(sizePolicy5);
        label_41->setMaximumSize(QSize(16777215, 15));
        label_41->setFont(font);

        gridLayout_17->addWidget(label_41, 0, 0, 1, 1);

        origin = new QSlider(Visualize);
        origin->setObjectName(QString::fromUtf8("origin"));
        origin->setSliderPosition(50);
        origin->setTracking(false);
        origin->setOrientation(Qt::Horizontal);
        origin->setInvertedAppearance(false);
        origin->setInvertedControls(false);

        gridLayout_17->addWidget(origin, 3, 1, 1, 3);

        xNormal = new QRadioButton(Visualize);
        buttonGroup = new QButtonGroup(MainWindow);
        buttonGroup->setObjectName(QString::fromUtf8("buttonGroup"));
        buttonGroup->addButton(xNormal);
        xNormal->setObjectName(QString::fromUtf8("xNormal"));

        gridLayout_17->addWidget(xNormal, 1, 1, 1, 1, Qt::AlignHCenter);

        verticalSpacer_22 = new QSpacerItem(20, 13, QSizePolicy::Minimum, QSizePolicy::Fixed);

        gridLayout_17->addItem(verticalSpacer_22, 9, 1, 1, 1);

        label_44 = new QLabel(Visualize);
        label_44->setObjectName(QString::fromUtf8("label_44"));
        sizePolicy5.setHeightForWidth(label_44->sizePolicy().hasHeightForWidth());
        label_44->setSizePolicy(sizePolicy5);
        label_44->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_17->addWidget(label_44, 1, 0, 1, 1);

        contourVOF = new QRadioButton(Visualize);
        buttonGroup_2 = new QButtonGroup(MainWindow);
        buttonGroup_2->setObjectName(QString::fromUtf8("buttonGroup_2"));
        buttonGroup_2->addButton(contourVOF);
        contourVOF->setObjectName(QString::fromUtf8("contourVOF"));
        contourVOF->setChecked(true);

        gridLayout_17->addWidget(contourVOF, 5, 1, 1, 1, Qt::AlignHCenter);

        verticalSpacer_20 = new QSpacerItem(20, 15, QSizePolicy::Minimum, QSizePolicy::Fixed);

        gridLayout_17->addItem(verticalSpacer_20, 4, 3, 1, 1);

        label_45 = new QLabel(Visualize);
        label_45->setObjectName(QString::fromUtf8("label_45"));
        sizePolicy5.setHeightForWidth(label_45->sizePolicy().hasHeightForWidth());
        label_45->setSizePolicy(sizePolicy5);
        label_45->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_17->addWidget(label_45, 5, 0, 1, 1);

        line_8 = new QFrame(Visualize);
        line_8->setObjectName(QString::fromUtf8("line_8"));
        line_8->setFrameShape(QFrame::HLine);
        line_8->setFrameShadow(QFrame::Sunken);

        gridLayout_17->addWidget(line_8, 10, 0, 1, 5);

        label_42 = new QLabel(Visualize);
        label_42->setObjectName(QString::fromUtf8("label_42"));
        label_42->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_17->addWidget(label_42, 3, 0, 1, 1);

        zNormal = new QRadioButton(Visualize);
        buttonGroup->addButton(zNormal);
        zNormal->setObjectName(QString::fromUtf8("zNormal"));

        gridLayout_17->addWidget(zNormal, 1, 3, 1, 1, Qt::AlignHCenter);

        yNormal = new QRadioButton(Visualize);
        buttonGroup->addButton(yNormal);
        yNormal->setObjectName(QString::fromUtf8("yNormal"));
        yNormal->setChecked(true);

        gridLayout_17->addWidget(yNormal, 1, 2, 1, 1, Qt::AlignHCenter);

        gridLayout_11 = new QGridLayout();
        gridLayout_11->setObjectName(QString::fromUtf8("gridLayout_11"));
        contourVorticity = new QRadioButton(Visualize);
        buttonGroup_2->addButton(contourVorticity);
        contourVorticity->setObjectName(QString::fromUtf8("contourVorticity"));

        gridLayout_11->addWidget(contourVorticity, 0, 0, 1, 1, Qt::AlignHCenter);

        contourVelocity = new QRadioButton(Visualize);
        buttonGroup_2->addButton(contourVelocity);
        contourVelocity->setObjectName(QString::fromUtf8("contourVelocity"));

        gridLayout_11->addWidget(contourVelocity, 0, 1, 1, 1, Qt::AlignHCenter);


        gridLayout_17->addLayout(gridLayout_11, 6, 1, 1, 3);

        originText = new QLabel(Visualize);
        originText->setObjectName(QString::fromUtf8("originText"));
        originText->setMaximumSize(QSize(16777215, 15));
        originText->setAlignment(Qt::AlignCenter);

        gridLayout_17->addWidget(originText, 4, 2, 1, 1);

        contourP = new QRadioButton(Visualize);
        buttonGroup_2->addButton(contourP);
        contourP->setObjectName(QString::fromUtf8("contourP"));

        gridLayout_17->addWidget(contourP, 5, 2, 1, 1, Qt::AlignHCenter);

        vis = new QVTKWidget(Visualize);
        vis->setObjectName(QString::fromUtf8("vis"));

        gridLayout_17->addWidget(vis, 13, 1, 2, 4);

        contourK = new QRadioButton(Visualize);
        buttonGroup_2->addButton(contourK);
        contourK->setObjectName(QString::fromUtf8("contourK"));

        gridLayout_17->addWidget(contourK, 5, 3, 1, 1, Qt::AlignHCenter);

        saveJPEG = new QPushButton(Visualize);
        saveJPEG->setObjectName(QString::fromUtf8("saveJPEG"));
        saveJPEG->setIcon(icon2);

        gridLayout_17->addWidget(saveJPEG, 13, 0, 1, 1);

        timesteps = new QListWidget(Visualize);
        timesteps->setObjectName(QString::fromUtf8("timesteps"));
        QSizePolicy sizePolicy6(QSizePolicy::Fixed, QSizePolicy::Expanding);
        sizePolicy6.setHorizontalStretch(0);
        sizePolicy6.setVerticalStretch(0);
        sizePolicy6.setHeightForWidth(timesteps->sizePolicy().hasHeightForWidth());
        timesteps->setSizePolicy(sizePolicy6);
        timesteps->setMaximumSize(QSize(150, 1600000));
        timesteps->setLayoutMode(QListView::SinglePass);

        gridLayout_17->addWidget(timesteps, 14, 0, 1, 1);

        showMesh = new QCheckBox(Visualize);
        showMesh->setObjectName(QString::fromUtf8("showMesh"));
        showMesh->setChecked(false);

        gridLayout_17->addWidget(showMesh, 8, 3, 1, 1, Qt::AlignHCenter);

        verticalSpacer_18 = new QSpacerItem(20, 20, QSizePolicy::Minimum, QSizePolicy::Fixed);

        gridLayout_17->addItem(verticalSpacer_18, 0, 2, 1, 1);

        showVectors = new QCheckBox(Visualize);
        showVectors->setObjectName(QString::fromUtf8("showVectors"));

        gridLayout_17->addWidget(showVectors, 8, 1, 1, 1, Qt::AlignHCenter);

        blockObstacles = new QCheckBox(Visualize);
        blockObstacles->setObjectName(QString::fromUtf8("blockObstacles"));
        blockObstacles->setChecked(true);

        gridLayout_17->addWidget(blockObstacles, 8, 2, 1, 1, Qt::AlignHCenter);

        label_46 = new QLabel(Visualize);
        label_46->setObjectName(QString::fromUtf8("label_46"));
        label_46->setFont(font);
        label_46->setAlignment(Qt::AlignCenter);

        gridLayout_17->addWidget(label_46, 12, 1, 1, 1);

        label_40 = new QLabel(Visualize);
        label_40->setObjectName(QString::fromUtf8("label_40"));
        sizePolicy5.setHeightForWidth(label_40->sizePolicy().hasHeightForWidth());
        label_40->setSizePolicy(sizePolicy5);
        label_40->setMaximumSize(QSize(16777215, 15));
        label_40->setFont(font);

        gridLayout_17->addWidget(label_40, 12, 0, 1, 1);

        showAxis = new QCheckBox(Visualize);
        showAxis->setObjectName(QString::fromUtf8("showAxis"));
        showAxis->setChecked(true);

        gridLayout_17->addWidget(showAxis, 8, 0, 1, 1, Qt::AlignHCenter);

        showLegend = new QCheckBox(Visualize);
        showLegend->setObjectName(QString::fromUtf8("showLegend"));
        showLegend->setChecked(true);

        gridLayout_17->addWidget(showLegend, 8, 4, 1, 1, Qt::AlignHCenter);


        gridLayout_18->addLayout(gridLayout_17, 0, 0, 1, 1);

        tabWidget->addTab(Visualize, QString());
        Visualize3d = new QWidget();
        Visualize3d->setObjectName(QString::fromUtf8("Visualize3d"));
        Visualize3d->setEnabled(false);
        gridLayout_23 = new QGridLayout(Visualize3d);
        gridLayout_23->setObjectName(QString::fromUtf8("gridLayout_23"));
        gridLayout_21 = new QGridLayout();
        gridLayout_21->setObjectName(QString::fromUtf8("gridLayout_21"));
        showAxis3d = new QCheckBox(Visualize3d);
        showAxis3d->setObjectName(QString::fromUtf8("showAxis3d"));
        showAxis3d->setChecked(true);

        gridLayout_21->addWidget(showAxis3d, 8, 1, 1, 1, Qt::AlignHCenter);

        extent3dText = new QLabel(Visualize3d);
        extent3dText->setObjectName(QString::fromUtf8("extent3dText"));
        extent3dText->setAlignment(Qt::AlignCenter);

        gridLayout_21->addWidget(extent3dText, 6, 2, 1, 1);

        vis3d = new QVTKWidget(Visualize3d);
        vis3d->setObjectName(QString::fromUtf8("vis3d"));

        gridLayout_21->addWidget(vis3d, 13, 1, 2, 4);

        origin3dText = new QLabel(Visualize3d);
        origin3dText->setObjectName(QString::fromUtf8("origin3dText"));
        origin3dText->setMaximumSize(QSize(16777215, 15));
        origin3dText->setAlignment(Qt::AlignCenter);

        gridLayout_21->addWidget(origin3dText, 4, 2, 1, 1);

        label_54 = new QLabel(Visualize3d);
        label_54->setObjectName(QString::fromUtf8("label_54"));
        label_54->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_21->addWidget(label_54, 3, 0, 1, 1);

        blockObstacles3d = new QCheckBox(Visualize3d);
        blockObstacles3d->setObjectName(QString::fromUtf8("blockObstacles3d"));
        blockObstacles3d->setChecked(true);

        gridLayout_21->addWidget(blockObstacles3d, 8, 2, 1, 1, Qt::AlignHCenter);

        showMesh3d = new QCheckBox(Visualize3d);
        showMesh3d->setObjectName(QString::fromUtf8("showMesh3d"));
        showMesh3d->setChecked(false);

        gridLayout_21->addWidget(showMesh3d, 8, 3, 1, 1, Qt::AlignHCenter);

        saveJPEG3d = new QPushButton(Visualize3d);
        saveJPEG3d->setObjectName(QString::fromUtf8("saveJPEG3d"));
        saveJPEG3d->setIcon(icon2);

        gridLayout_21->addWidget(saveJPEG3d, 13, 0, 1, 1);

        yNormal3d = new QRadioButton(Visualize3d);
        yNormal3d->setObjectName(QString::fromUtf8("yNormal3d"));
        yNormal3d->setChecked(true);

        gridLayout_21->addWidget(yNormal3d, 1, 2, 1, 1, Qt::AlignHCenter);

        verticalSpacer_27 = new QSpacerItem(20, 15, QSizePolicy::Minimum, QSizePolicy::Fixed);

        gridLayout_21->addItem(verticalSpacer_27, 4, 3, 1, 1);

        verticalSpacer_23 = new QSpacerItem(20, 5, QSizePolicy::Minimum, QSizePolicy::Fixed);

        gridLayout_21->addItem(verticalSpacer_23, 2, 2, 1, 1);

        verticalSpacer_25 = new QSpacerItem(20, 10, QSizePolicy::Minimum, QSizePolicy::Fixed);

        gridLayout_21->addItem(verticalSpacer_25, 7, 3, 1, 1);

        xNormal3d = new QRadioButton(Visualize3d);
        xNormal3d->setObjectName(QString::fromUtf8("xNormal3d"));

        gridLayout_21->addWidget(xNormal3d, 1, 1, 1, 1, Qt::AlignHCenter);

        verticalSpacer_26 = new QSpacerItem(20, 13, QSizePolicy::Minimum, QSizePolicy::Fixed);

        gridLayout_21->addItem(verticalSpacer_26, 9, 1, 1, 1);

        label_51 = new QLabel(Visualize3d);
        label_51->setObjectName(QString::fromUtf8("label_51"));
        sizePolicy5.setHeightForWidth(label_51->sizePolicy().hasHeightForWidth());
        label_51->setSizePolicy(sizePolicy5);
        label_51->setMaximumSize(QSize(16777215, 15));
        label_51->setFont(font);

        gridLayout_21->addWidget(label_51, 0, 0, 1, 1);

        origin3d = new QSlider(Visualize3d);
        origin3d->setObjectName(QString::fromUtf8("origin3d"));
        origin3d->setValue(0);
        origin3d->setSliderPosition(0);
        origin3d->setTracking(false);
        origin3d->setOrientation(Qt::Horizontal);
        origin3d->setInvertedAppearance(false);
        origin3d->setInvertedControls(false);

        gridLayout_21->addWidget(origin3d, 3, 1, 1, 3);

        label_53 = new QLabel(Visualize3d);
        label_53->setObjectName(QString::fromUtf8("label_53"));
        sizePolicy5.setHeightForWidth(label_53->sizePolicy().hasHeightForWidth());
        label_53->setSizePolicy(sizePolicy5);
        label_53->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_21->addWidget(label_53, 5, 0, 1, 1);

        line_9 = new QFrame(Visualize3d);
        line_9->setObjectName(QString::fromUtf8("line_9"));
        line_9->setFrameShape(QFrame::HLine);
        line_9->setFrameShadow(QFrame::Sunken);

        gridLayout_21->addWidget(line_9, 10, 0, 1, 5);

        horizontalLayout_4 = new QHBoxLayout();
        horizontalLayout_4->setObjectName(QString::fromUtf8("horizontalLayout_4"));
        horizontalSpacer_19 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_4->addItem(horizontalSpacer_19);

        horizontalSpacer_20 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_4->addItem(horizontalSpacer_20);


        gridLayout_21->addLayout(horizontalLayout_4, 12, 2, 1, 3);

        zNormal3d = new QRadioButton(Visualize3d);
        zNormal3d->setObjectName(QString::fromUtf8("zNormal3d"));

        gridLayout_21->addWidget(zNormal3d, 1, 3, 1, 1, Qt::AlignHCenter);

        timesteps3d = new QListWidget(Visualize3d);
        timesteps3d->setObjectName(QString::fromUtf8("timesteps3d"));
        sizePolicy6.setHeightForWidth(timesteps3d->sizePolicy().hasHeightForWidth());
        timesteps3d->setSizePolicy(sizePolicy6);
        timesteps3d->setMaximumSize(QSize(150, 1600000));
        timesteps3d->setLayoutMode(QListView::SinglePass);

        gridLayout_21->addWidget(timesteps3d, 14, 0, 1, 1);

        verticalSpacer_28 = new QSpacerItem(20, 20, QSizePolicy::Minimum, QSizePolicy::Fixed);

        gridLayout_21->addItem(verticalSpacer_28, 0, 2, 1, 1);

        label_55 = new QLabel(Visualize3d);
        label_55->setObjectName(QString::fromUtf8("label_55"));
        label_55->setFont(font);
        label_55->setAlignment(Qt::AlignCenter);

        gridLayout_21->addWidget(label_55, 12, 1, 1, 1);

        label_56 = new QLabel(Visualize3d);
        label_56->setObjectName(QString::fromUtf8("label_56"));
        sizePolicy5.setHeightForWidth(label_56->sizePolicy().hasHeightForWidth());
        label_56->setSizePolicy(sizePolicy5);
        label_56->setMaximumSize(QSize(16777215, 15));
        label_56->setFont(font);

        gridLayout_21->addWidget(label_56, 12, 0, 1, 1);

        label_52 = new QLabel(Visualize3d);
        label_52->setObjectName(QString::fromUtf8("label_52"));
        sizePolicy5.setHeightForWidth(label_52->sizePolicy().hasHeightForWidth());
        label_52->setSizePolicy(sizePolicy5);
        label_52->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_21->addWidget(label_52, 1, 0, 1, 1);

        extent3d = new QSlider(Visualize3d);
        extent3d->setObjectName(QString::fromUtf8("extent3d"));
        extent3d->setValue(99);
        extent3d->setTracking(false);
        extent3d->setOrientation(Qt::Horizontal);

        gridLayout_21->addWidget(extent3d, 5, 1, 1, 3);


        gridLayout_23->addLayout(gridLayout_21, 0, 0, 1, 1);

        tabWidget->addTab(Visualize3d, QString());

        gridLayout_3->addWidget(tabWidget, 0, 0, 1, 1);

        MainWindow->setCentralWidget(centralwidget);
        menubar = new QMenuBar(MainWindow);
        menubar->setObjectName(QString::fromUtf8("menubar"));
        menubar->setGeometry(QRect(0, 0, 950, 25));
        menuFile = new QMenu(menubar);
        menuFile->setObjectName(QString::fromUtf8("menuFile"));
        menuHelp = new QMenu(menubar);
        menuHelp->setObjectName(QString::fromUtf8("menuHelp"));
        MainWindow->setMenuBar(menubar);
        statusbar = new QStatusBar(MainWindow);
        statusbar->setObjectName(QString::fromUtf8("statusbar"));
        MainWindow->setStatusBar(statusbar);

        menubar->addAction(menuFile->menuAction());
        menubar->addAction(menuHelp->menuAction());
        menuFile->addAction(actionNew);
        menuFile->addAction(actionOpen);
        menuFile->addSeparator();
        menuFile->addAction(actionSave);
        menuFile->addSeparator();
        menuFile->addAction(actionQuit);
        menuHelp->addAction(actionAbout);

        retranslateUi(MainWindow);

        tabWidget->setCurrentIndex(0);


        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "Civil CFD", 0, QApplication::UnicodeUTF8));
        actionNew->setText(QApplication::translate("MainWindow", "New", 0, QApplication::UnicodeUTF8));
        actionNew->setShortcut(QApplication::translate("MainWindow", "Ctrl+N", 0, QApplication::UnicodeUTF8));
        actionOpen->setText(QApplication::translate("MainWindow", "Open", 0, QApplication::UnicodeUTF8));
        actionOpen->setShortcut(QApplication::translate("MainWindow", "Ctrl+O", 0, QApplication::UnicodeUTF8));
        actionSave->setText(QApplication::translate("MainWindow", "Save", 0, QApplication::UnicodeUTF8));
        actionSave->setShortcut(QApplication::translate("MainWindow", "Ctrl+S", 0, QApplication::UnicodeUTF8));
        actionQuit->setText(QApplication::translate("MainWindow", "Quit", 0, QApplication::UnicodeUTF8));
        actionAbout->setText(QApplication::translate("MainWindow", "About", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("MainWindow", "Simulation Path", 0, QApplication::UnicodeUTF8));
        label_5->setText(QApplication::translate("MainWindow", "Mesh Element Size", 0, QApplication::UnicodeUTF8));
        New->setText(QApplication::translate("MainWindow", "New", 0, QApplication::UnicodeUTF8));
        Save->setText(QApplication::translate("MainWindow", "Save", 0, QApplication::UnicodeUTF8));
        Clear->setText(QApplication::translate("MainWindow", "Clear Selection", 0, QApplication::UnicodeUTF8));
        label_9->setText(QApplication::translate("MainWindow", "<html><head/><body><p><span style=\" font-weight:600;\">Results</span></p></body></html>", 0, QApplication::UnicodeUTF8));
        InspectCell->setText(QApplication::translate("MainWindow", "Inspect Cell", 0, QApplication::UnicodeUTF8));
        Delete->setText(QApplication::translate("MainWindow", "Delete", 0, QApplication::UnicodeUTF8));
        SelectAll->setText(QApplication::translate("MainWindow", "Select All", 0, QApplication::UnicodeUTF8));
        label_10->setText(QApplication::translate("MainWindow", "Timesteps Saved", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("MainWindow", "Simulation Name", 0, QApplication::UnicodeUTF8));
        label_4->setText(QApplication::translate("MainWindow", "Mesh Dimensions", 0, QApplication::UnicodeUTF8));
        label_6->setText(QApplication::translate("MainWindow", "Number of Elements", 0, QApplication::UnicodeUTF8));
        label_3->setText(QApplication::translate("MainWindow", "<html><head/><body><p><span style=\" font-weight:600;\">Simulation Details</span></p></body></html>", 0, QApplication::UnicodeUTF8));
        label_7->setText(QApplication::translate("MainWindow", "Geometry File", 0, QApplication::UnicodeUTF8));
        tabWidget->setTabText(tabWidget->indexOf(Information), QApplication::translate("MainWindow", "Information", 0, QApplication::UnicodeUTF8));
        calcRough->setText(QApplication::translate("MainWindow", "Calculate", 0, QApplication::UnicodeUTF8));
        kEpsilon->setText(QApplication::translate("MainWindow", "k-Epsilon", 0, QApplication::UnicodeUTF8));
        label_16->setText(QApplication::translate("MainWindow", "y", 0, QApplication::UnicodeUTF8));
        label_20->setText(QApplication::translate("MainWindow", "Roughness height (m)", 0, QApplication::UnicodeUTF8));
        label_14->setText(QApplication::translate("MainWindow", "Density (kg/m^3)", 0, QApplication::UnicodeUTF8));
        label_18->setText(QApplication::translate("MainWindow", "Turbulence", 0, QApplication::UnicodeUTF8));
        Laminar->setText(QApplication::translate("MainWindow", "Laminar Flow Only", 0, QApplication::UnicodeUTF8));
        label_11->setText(QApplication::translate("MainWindow", "Constants", 0, QApplication::UnicodeUTF8));
        label_19->setText(QApplication::translate("MainWindow", "Length scale (%)", 0, QApplication::UnicodeUTF8));
        water20C->setText(QApplication::translate("MainWindow", "Water at 20 C", 0, QApplication::UnicodeUTF8));
        defaultLength->setText(QApplication::translate("MainWindow", "Default", 0, QApplication::UnicodeUTF8));
        label_15->setText(QApplication::translate("MainWindow", "x", 0, QApplication::UnicodeUTF8));
        label_12->setText(QApplication::translate("MainWindow", "Gravity (m/s^2)", 0, QApplication::UnicodeUTF8));
        label_13->setText(QApplication::translate("MainWindow", "Kinematic Viscocity (m^2/s)", 0, QApplication::UnicodeUTF8));
        label_17->setText(QApplication::translate("MainWindow", "z", 0, QApplication::UnicodeUTF8));
        rough->setPlainText(QApplication::translate("MainWindow", "0.00161", 0, QApplication::UnicodeUTF8));
        earthGravity->setText(QApplication::translate("MainWindow", "Earth Gravity", 0, QApplication::UnicodeUTF8));
        label_8->setText(QApplication::translate("MainWindow", "Length (m)", 0, QApplication::UnicodeUTF8));
        domainLength->setText(QApplication::translate("MainWindow", "Entire Domain", 0, QApplication::UnicodeUTF8));
        tabWidget->setTabText(tabWidget->indexOf(Models), QApplication::translate("MainWindow", "Models", 0, QApplication::UnicodeUTF8));
        MeshUndo->setText(QApplication::translate("MainWindow", "Revert Changes", 0, QApplication::UnicodeUTF8));
        QTreeWidgetItem *___qtreewidgetitem = MeshParameters->headerItem();
        ___qtreewidgetitem->setText(1, QApplication::translate("MainWindow", "Value", 0, QApplication::UnicodeUTF8));
        ___qtreewidgetitem->setText(0, QApplication::translate("MainWindow", "Parameter", 0, QApplication::UnicodeUTF8));

        const bool __sortingEnabled = MeshParameters->isSortingEnabled();
        MeshParameters->setSortingEnabled(false);
        QTreeWidgetItem *___qtreewidgetitem1 = MeshParameters->topLevelItem(0);
        ___qtreewidgetitem1->setText(0, QApplication::translate("MainWindow", "Mesh Cell Dimensions (m)", 0, QApplication::UnicodeUTF8));
        QTreeWidgetItem *___qtreewidgetitem2 = ___qtreewidgetitem1->child(0);
        ___qtreewidgetitem2->setText(0, QApplication::translate("MainWindow", "X", 0, QApplication::UnicodeUTF8));
        QTreeWidgetItem *___qtreewidgetitem3 = ___qtreewidgetitem1->child(1);
        ___qtreewidgetitem3->setText(0, QApplication::translate("MainWindow", "Y", 0, QApplication::UnicodeUTF8));
        QTreeWidgetItem *___qtreewidgetitem4 = ___qtreewidgetitem1->child(2);
        ___qtreewidgetitem4->setText(0, QApplication::translate("MainWindow", "Z", 0, QApplication::UnicodeUTF8));
        QTreeWidgetItem *___qtreewidgetitem5 = MeshParameters->topLevelItem(1);
        ___qtreewidgetitem5->setText(0, QApplication::translate("MainWindow", "Number of Mesh Cells", 0, QApplication::UnicodeUTF8));
        QTreeWidgetItem *___qtreewidgetitem6 = ___qtreewidgetitem5->child(0);
        ___qtreewidgetitem6->setText(0, QApplication::translate("MainWindow", "X", 0, QApplication::UnicodeUTF8));
        QTreeWidgetItem *___qtreewidgetitem7 = ___qtreewidgetitem5->child(1);
        ___qtreewidgetitem7->setText(0, QApplication::translate("MainWindow", "Y", 0, QApplication::UnicodeUTF8));
        QTreeWidgetItem *___qtreewidgetitem8 = ___qtreewidgetitem5->child(2);
        ___qtreewidgetitem8->setText(0, QApplication::translate("MainWindow", "Z", 0, QApplication::UnicodeUTF8));
        QTreeWidgetItem *___qtreewidgetitem9 = MeshParameters->topLevelItem(2);
        ___qtreewidgetitem9->setText(0, QApplication::translate("MainWindow", "Mesh Origin (m)", 0, QApplication::UnicodeUTF8));
        QTreeWidgetItem *___qtreewidgetitem10 = ___qtreewidgetitem9->child(0);
        ___qtreewidgetitem10->setText(0, QApplication::translate("MainWindow", "X", 0, QApplication::UnicodeUTF8));
        QTreeWidgetItem *___qtreewidgetitem11 = ___qtreewidgetitem9->child(1);
        ___qtreewidgetitem11->setText(0, QApplication::translate("MainWindow", "Y", 0, QApplication::UnicodeUTF8));
        QTreeWidgetItem *___qtreewidgetitem12 = ___qtreewidgetitem9->child(2);
        ___qtreewidgetitem12->setText(0, QApplication::translate("MainWindow", "Z", 0, QApplication::UnicodeUTF8));
        MeshParameters->setSortingEnabled(__sortingEnabled);

        MeshUpdate->setText(QApplication::translate("MainWindow", "Update Display", 0, QApplication::UnicodeUTF8));
        tabWidget->setTabText(tabWidget->indexOf(Mesh), QApplication::translate("MainWindow", "Mesh", 0, QApplication::UnicodeUTF8));
        label_22->setText(QApplication::translate("MainWindow", "Mesh Interior Point", 0, QApplication::UnicodeUTF8));
        label_25->setText(QApplication::translate("MainWindow", "Z", 0, QApplication::UnicodeUTF8));
        label_24->setText(QApplication::translate("MainWindow", "Y", 0, QApplication::UnicodeUTF8));
        label_23->setText(QApplication::translate("MainWindow", "X", 0, QApplication::UnicodeUTF8));
        label_21->setText(QApplication::translate("MainWindow", "Geometry File", 0, QApplication::UnicodeUTF8));
        STLRender->setText(QApplication::translate("MainWindow", "Render", 0, QApplication::UnicodeUTF8));
        STLOpen->setText(QApplication::translate("MainWindow", "Open", 0, QApplication::UnicodeUTF8));
        tabWidget->setTabText(tabWidget->indexOf(Geometry), QApplication::translate("MainWindow", "Geometry", 0, QApplication::UnicodeUTF8));
        QTreeWidgetItem *___qtreewidgetitem13 = BoundaryTree->headerItem();
        ___qtreewidgetitem13->setText(4, QApplication::translate("MainWindow", "Turbulence", 0, QApplication::UnicodeUTF8));
        ___qtreewidgetitem13->setText(3, QApplication::translate("MainWindow", "Value", 0, QApplication::UnicodeUTF8));
        ___qtreewidgetitem13->setText(2, QApplication::translate("MainWindow", "Extents", 0, QApplication::UnicodeUTF8));
        ___qtreewidgetitem13->setText(1, QApplication::translate("MainWindow", "Type", 0, QApplication::UnicodeUTF8));
        ___qtreewidgetitem13->setText(0, QApplication::translate("MainWindow", "Parameter", 0, QApplication::UnicodeUTF8));

        const bool __sortingEnabled1 = BoundaryTree->isSortingEnabled();
        BoundaryTree->setSortingEnabled(false);
        QTreeWidgetItem *___qtreewidgetitem14 = BoundaryTree->topLevelItem(0);
        ___qtreewidgetitem14->setText(0, QApplication::translate("MainWindow", "Wall Boundaries", 0, QApplication::UnicodeUTF8));
        QTreeWidgetItem *___qtreewidgetitem15 = ___qtreewidgetitem14->child(0);
        ___qtreewidgetitem15->setText(0, QApplication::translate("MainWindow", "West", 0, QApplication::UnicodeUTF8));
        QTreeWidgetItem *___qtreewidgetitem16 = ___qtreewidgetitem14->child(1);
        ___qtreewidgetitem16->setText(0, QApplication::translate("MainWindow", "East", 0, QApplication::UnicodeUTF8));
        QTreeWidgetItem *___qtreewidgetitem17 = ___qtreewidgetitem14->child(2);
        ___qtreewidgetitem17->setText(0, QApplication::translate("MainWindow", "South", 0, QApplication::UnicodeUTF8));
        QTreeWidgetItem *___qtreewidgetitem18 = ___qtreewidgetitem14->child(3);
        ___qtreewidgetitem18->setText(0, QApplication::translate("MainWindow", "North", 0, QApplication::UnicodeUTF8));
        QTreeWidgetItem *___qtreewidgetitem19 = ___qtreewidgetitem14->child(4);
        ___qtreewidgetitem19->setText(0, QApplication::translate("MainWindow", "Bottom", 0, QApplication::UnicodeUTF8));
        QTreeWidgetItem *___qtreewidgetitem20 = ___qtreewidgetitem14->child(5);
        ___qtreewidgetitem20->setText(0, QApplication::translate("MainWindow", "Top", 0, QApplication::UnicodeUTF8));
        QTreeWidgetItem *___qtreewidgetitem21 = BoundaryTree->topLevelItem(1);
        ___qtreewidgetitem21->setText(0, QApplication::translate("MainWindow", "Special Boundaries", 0, QApplication::UnicodeUTF8));
        QTreeWidgetItem *___qtreewidgetitem22 = ___qtreewidgetitem21->child(0);
        ___qtreewidgetitem22->setText(0, QApplication::translate("MainWindow", "West", 0, QApplication::UnicodeUTF8));
        QTreeWidgetItem *___qtreewidgetitem23 = ___qtreewidgetitem21->child(1);
        ___qtreewidgetitem23->setText(0, QApplication::translate("MainWindow", "East", 0, QApplication::UnicodeUTF8));
        QTreeWidgetItem *___qtreewidgetitem24 = ___qtreewidgetitem21->child(2);
        ___qtreewidgetitem24->setText(0, QApplication::translate("MainWindow", "South", 0, QApplication::UnicodeUTF8));
        QTreeWidgetItem *___qtreewidgetitem25 = ___qtreewidgetitem21->child(3);
        ___qtreewidgetitem25->setText(0, QApplication::translate("MainWindow", "North", 0, QApplication::UnicodeUTF8));
        QTreeWidgetItem *___qtreewidgetitem26 = ___qtreewidgetitem21->child(4);
        ___qtreewidgetitem26->setText(0, QApplication::translate("MainWindow", "Bottom", 0, QApplication::UnicodeUTF8));
        QTreeWidgetItem *___qtreewidgetitem27 = ___qtreewidgetitem21->child(5);
        ___qtreewidgetitem27->setText(0, QApplication::translate("MainWindow", "Top", 0, QApplication::UnicodeUTF8));
        BoundaryTree->setSortingEnabled(__sortingEnabled1);

        EditBoundary->setText(QApplication::translate("MainWindow", "Edit Boundary", 0, QApplication::UnicodeUTF8));
        AddSpecialBoundary->setText(QApplication::translate("MainWindow", "Add Special Boundary", 0, QApplication::UnicodeUTF8));
        RemoveSpecialBoundary->setText(QApplication::translate("MainWindow", "Remove Special Boundary", 0, QApplication::UnicodeUTF8));
        hideBoundaryMesh->setText(QApplication::translate("MainWindow", "Hide Mesh", 0, QApplication::UnicodeUTF8));
        hideBoundaryDomain->setText(QApplication::translate("MainWindow", "Hide Rendered Domain", 0, QApplication::UnicodeUTF8));
        tabWidget->setTabText(tabWidget->indexOf(Boundaries), QApplication::translate("MainWindow", "Boundaries", 0, QApplication::UnicodeUTF8));
        EditBaffle->setText(QApplication::translate("MainWindow", "Edit Baffle", 0, QApplication::UnicodeUTF8));
        AddBaffle->setText(QApplication::translate("MainWindow", "Add Baffle", 0, QApplication::UnicodeUTF8));
        RemoveBaffle->setText(QApplication::translate("MainWindow", "Remove Baffle", 0, QApplication::UnicodeUTF8));
        QTreeWidgetItem *___qtreewidgetitem28 = BaffleTree->headerItem();
        ___qtreewidgetitem28->setText(4, QApplication::translate("MainWindow", "Value", 0, QApplication::UnicodeUTF8));
        ___qtreewidgetitem28->setText(3, QApplication::translate("MainWindow", "Position", 0, QApplication::UnicodeUTF8));
        ___qtreewidgetitem28->setText(2, QApplication::translate("MainWindow", "Extents", 0, QApplication::UnicodeUTF8));
        ___qtreewidgetitem28->setText(1, QApplication::translate("MainWindow", "Type", 0, QApplication::UnicodeUTF8));
        ___qtreewidgetitem28->setText(0, QApplication::translate("MainWindow", "Parameter", 0, QApplication::UnicodeUTF8));

        const bool __sortingEnabled2 = BaffleTree->isSortingEnabled();
        BaffleTree->setSortingEnabled(false);
        QTreeWidgetItem *___qtreewidgetitem29 = BaffleTree->topLevelItem(0);
        ___qtreewidgetitem29->setText(0, QApplication::translate("MainWindow", "Axis", 0, QApplication::UnicodeUTF8));
        QTreeWidgetItem *___qtreewidgetitem30 = ___qtreewidgetitem29->child(0);
        ___qtreewidgetitem30->setText(0, QApplication::translate("MainWindow", "x", 0, QApplication::UnicodeUTF8));
        QTreeWidgetItem *___qtreewidgetitem31 = ___qtreewidgetitem29->child(1);
        ___qtreewidgetitem31->setText(0, QApplication::translate("MainWindow", "y", 0, QApplication::UnicodeUTF8));
        QTreeWidgetItem *___qtreewidgetitem32 = ___qtreewidgetitem29->child(2);
        ___qtreewidgetitem32->setText(0, QApplication::translate("MainWindow", "z", 0, QApplication::UnicodeUTF8));
        BaffleTree->setSortingEnabled(__sortingEnabled2);

        hideBaffleMesh->setText(QApplication::translate("MainWindow", "Hide Mesh", 0, QApplication::UnicodeUTF8));
        hideBaffleDomain->setText(QApplication::translate("MainWindow", "Hide Rendered Domain", 0, QApplication::UnicodeUTF8));
        tabWidget->setTabText(tabWidget->indexOf(Baffles), QApplication::translate("MainWindow", "Baffles", 0, QApplication::UnicodeUTF8));
        label_26->setText(QApplication::translate("MainWindow", "Velocity", 0, QApplication::UnicodeUTF8));
        label_30->setText(QApplication::translate("MainWindow", "V (m/s)", 0, QApplication::UnicodeUTF8));
        label_29->setText(QApplication::translate("MainWindow", "U (m/s)", 0, QApplication::UnicodeUTF8));
        initialU->setPlainText(QApplication::translate("MainWindow", "0", 0, QApplication::UnicodeUTF8));
        initialW->setPlainText(QApplication::translate("MainWindow", "0", 0, QApplication::UnicodeUTF8));
        label_28->setText(QApplication::translate("MainWindow", "Pressure", 0, QApplication::UnicodeUTF8));
        label_32->setText(QApplication::translate("MainWindow", "Height (m)", 0, QApplication::UnicodeUTF8));
        label_27->setText(QApplication::translate("MainWindow", "Fluid", 0, QApplication::UnicodeUTF8));
        initialV->setPlainText(QApplication::translate("MainWindow", "0", 0, QApplication::UnicodeUTF8));
        label_31->setText(QApplication::translate("MainWindow", "W (m/s)", 0, QApplication::UnicodeUTF8));
        hydrostatic->setText(QApplication::translate("MainWindow", "Hydrostatic Pressure", 0, QApplication::UnicodeUTF8));
        label_33->setText(QApplication::translate("MainWindow", "Turbulence", 0, QApplication::UnicodeUTF8));
        label_34->setText(QApplication::translate("MainWindow", "Kinetic energy", 0, QApplication::UnicodeUTF8));
        initialk->setPlainText(QApplication::translate("MainWindow", "0.001", 0, QApplication::UnicodeUTF8));
        label_43->setText(QApplication::translate("MainWindow", "Internal Fill Points", 0, QApplication::UnicodeUTF8));
        QTableWidgetItem *___qtablewidgetitem = fillPoints->horizontalHeaderItem(0);
        ___qtablewidgetitem->setText(QApplication::translate("MainWindow", "x", 0, QApplication::UnicodeUTF8));
        QTableWidgetItem *___qtablewidgetitem1 = fillPoints->horizontalHeaderItem(1);
        ___qtablewidgetitem1->setText(QApplication::translate("MainWindow", "y", 0, QApplication::UnicodeUTF8));
        QTableWidgetItem *___qtablewidgetitem2 = fillPoints->horizontalHeaderItem(2);
        ___qtablewidgetitem2->setText(QApplication::translate("MainWindow", "z", 0, QApplication::UnicodeUTF8));

        const bool __sortingEnabled3 = fillPoints->isSortingEnabled();
        fillPoints->setSortingEnabled(false);
        fillPoints->setSortingEnabled(__sortingEnabled3);

        tabWidget->setTabText(tabWidget->indexOf(InitialConditions), QApplication::translate("MainWindow", "Initial Conditions", 0, QApplication::UnicodeUTF8));
        label_37->setText(QApplication::translate("MainWindow", "End Time", 0, QApplication::UnicodeUTF8));
        t->clear();
        t->insertItems(0, QStringList()
         << QApplication::translate("MainWindow", "0", 0, QApplication::UnicodeUTF8)
        );
        label_39->setText(QApplication::translate("MainWindow", "Write Interval", 0, QApplication::UnicodeUTF8));
        autot->setText(QApplication::translate("MainWindow", "Automatic timestep adjustment", 0, QApplication::UnicodeUTF8));
        writet->setPlainText(QApplication::translate("MainWindow", "1.0", 0, QApplication::UnicodeUTF8));
        RunSimulation->setText(QApplication::translate("MainWindow", "Simulate", 0, QApplication::UnicodeUTF8));
        delt->setPlainText(QApplication::translate("MainWindow", "0.001", 0, QApplication::UnicodeUTF8));
        label_38->setText(QApplication::translate("MainWindow", "Initial Timestep", 0, QApplication::UnicodeUTF8));
        label_47->setText(QApplication::translate("MainWindow", "Implicit Solver", 0, QApplication::UnicodeUTF8));
        label_36->setText(QApplication::translate("MainWindow", "Start Time", 0, QApplication::UnicodeUTF8));
        endt->setPlainText(QApplication::translate("MainWindow", "100", 0, QApplication::UnicodeUTF8));
        label_35->setText(QApplication::translate("MainWindow", "Time", 0, QApplication::UnicodeUTF8));
        label_49->setText(QApplication::translate("MainWindow", "Processes", 0, QApplication::UnicodeUTF8));
        processes->clear();
        processes->insertItems(0, QStringList()
         << QApplication::translate("MainWindow", "2", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "3", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "4", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "6", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "8", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "16", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "32", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "64", 0, QApplication::UnicodeUTF8)
        );
        tabWidget->setTabText(tabWidget->indexOf(Simulate), QApplication::translate("MainWindow", "Simulate", 0, QApplication::UnicodeUTF8));
        updateRange->setText(QApplication::translate("MainWindow", "Rescale Range from:", 0, QApplication::UnicodeUTF8));
        from->setPlainText(QApplication::translate("MainWindow", "0", 0, QApplication::UnicodeUTF8));
        label_48->setText(QApplication::translate("MainWindow", "to", 0, QApplication::UnicodeUTF8));
        to->setPlainText(QApplication::translate("MainWindow", "1", 0, QApplication::UnicodeUTF8));
        label_41->setText(QApplication::translate("MainWindow", "Section Plane", 0, QApplication::UnicodeUTF8));
        xNormal->setText(QApplication::translate("MainWindow", "X-Normal", 0, QApplication::UnicodeUTF8));
        label_44->setText(QApplication::translate("MainWindow", "Orientation", 0, QApplication::UnicodeUTF8));
        contourVOF->setText(QApplication::translate("MainWindow", "Fluid fraction", 0, QApplication::UnicodeUTF8));
        label_45->setText(QApplication::translate("MainWindow", "Contour", 0, QApplication::UnicodeUTF8));
        label_42->setText(QApplication::translate("MainWindow", "Origin", 0, QApplication::UnicodeUTF8));
        zNormal->setText(QApplication::translate("MainWindow", "Z-Normal", 0, QApplication::UnicodeUTF8));
        yNormal->setText(QApplication::translate("MainWindow", "Y-Normal", 0, QApplication::UnicodeUTF8));
        contourVorticity->setText(QApplication::translate("MainWindow", "Vorticity", 0, QApplication::UnicodeUTF8));
        contourVelocity->setText(QApplication::translate("MainWindow", "Velocity", 0, QApplication::UnicodeUTF8));
        originText->setText(QApplication::translate("MainWindow", "0", 0, QApplication::UnicodeUTF8));
        contourP->setText(QApplication::translate("MainWindow", "Pressure", 0, QApplication::UnicodeUTF8));
        contourK->setText(QApplication::translate("MainWindow", "Turbulence (k)", 0, QApplication::UnicodeUTF8));
        saveJPEG->setText(QApplication::translate("MainWindow", "Save Image", 0, QApplication::UnicodeUTF8));
        showMesh->setText(QApplication::translate("MainWindow", "Show Mesh", 0, QApplication::UnicodeUTF8));
        showVectors->setText(QApplication::translate("MainWindow", "Show Velocity Vectors", 0, QApplication::UnicodeUTF8));
        blockObstacles->setText(QApplication::translate("MainWindow", "Show Obstacles", 0, QApplication::UnicodeUTF8));
        label_46->setText(QApplication::translate("MainWindow", "Results", 0, QApplication::UnicodeUTF8));
        label_40->setText(QApplication::translate("MainWindow", "Timestep", 0, QApplication::UnicodeUTF8));
        showAxis->setText(QApplication::translate("MainWindow", "Show Axis", 0, QApplication::UnicodeUTF8));
        showLegend->setText(QApplication::translate("MainWindow", "Show Legend", 0, QApplication::UnicodeUTF8));
        tabWidget->setTabText(tabWidget->indexOf(Visualize), QApplication::translate("MainWindow", "Visualize", 0, QApplication::UnicodeUTF8));
        showAxis3d->setText(QApplication::translate("MainWindow", "Show Axis", 0, QApplication::UnicodeUTF8));
        extent3dText->setText(QApplication::translate("MainWindow", "0", 0, QApplication::UnicodeUTF8));
        origin3dText->setText(QApplication::translate("MainWindow", "0", 0, QApplication::UnicodeUTF8));
        label_54->setText(QApplication::translate("MainWindow", "Origin", 0, QApplication::UnicodeUTF8));
        blockObstacles3d->setText(QApplication::translate("MainWindow", "Show Obstacles", 0, QApplication::UnicodeUTF8));
        showMesh3d->setText(QApplication::translate("MainWindow", "Show Mesh", 0, QApplication::UnicodeUTF8));
        saveJPEG3d->setText(QApplication::translate("MainWindow", "Save Image", 0, QApplication::UnicodeUTF8));
        yNormal3d->setText(QApplication::translate("MainWindow", "Y-Normal", 0, QApplication::UnicodeUTF8));
        xNormal3d->setText(QApplication::translate("MainWindow", "X-Normal", 0, QApplication::UnicodeUTF8));
        label_51->setText(QApplication::translate("MainWindow", "Iso-Surface Clipping", 0, QApplication::UnicodeUTF8));
        label_53->setText(QApplication::translate("MainWindow", "Extent", 0, QApplication::UnicodeUTF8));
        zNormal3d->setText(QApplication::translate("MainWindow", "Z-Normal", 0, QApplication::UnicodeUTF8));
        label_55->setText(QApplication::translate("MainWindow", "Results", 0, QApplication::UnicodeUTF8));
        label_56->setText(QApplication::translate("MainWindow", "Timestep", 0, QApplication::UnicodeUTF8));
        label_52->setText(QApplication::translate("MainWindow", "Orientation", 0, QApplication::UnicodeUTF8));
        tabWidget->setTabText(tabWidget->indexOf(Visualize3d), QApplication::translate("MainWindow", "3D Surface", 0, QApplication::UnicodeUTF8));
        menuFile->setTitle(QApplication::translate("MainWindow", "File", 0, QApplication::UnicodeUTF8));
        menuHelp->setTitle(QApplication::translate("MainWindow", "Help", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_CIVLCFD_H
