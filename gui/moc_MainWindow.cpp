/****************************************************************************
** Meta object code from reading C++ file 'MainWindow.h'
**
** Created by: The Qt Meta Object Compiler version 63 (Qt 4.8.7)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "MainWindow.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'MainWindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.7. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_MainWindow[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
      51,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      11,   38,   38,   38, 0x08,
      39,   38,   38,   38, 0x08,
      62,   38,   38,   38, 0x08,
      83,   38,   38,   38, 0x08,
     105,   38,   38,   38, 0x08,
     129,   38,   38,   38, 0x08,
     146,   38,   38,   38, 0x08,
     164,   38,   38,   38, 0x08,
     190,   38,   38,   38, 0x08,
     216,   38,   38,   38, 0x08,
     242,   38,   38,   38, 0x08,
     264,   38,   38,   38, 0x08,
     285,  343,   38,   38, 0x08,
     355,   38,   38,   38, 0x08,
     381,  343,   38,   38, 0x08,
     437,   38,   38,   38, 0x08,
     469,  343,   38,   38, 0x08,
     521,  343,   38,   38, 0x08,
     571,   38,   38,   38, 0x08,
     606,   38,   38,   38, 0x08,
     630,  343,   38,   38, 0x08,
     684,   38,   38,   38, 0x08,
     707,  343,   38,   38, 0x08,
     757,  343,   38,   38, 0x08,
     805,   38,   38,   38, 0x08,
     831,   38,   38,   38, 0x08,
     856,   38,   38,   38, 0x08,
     878,   38,   38,   38, 0x08,
     904,   38,   38,   38, 0x08,
     930,   38,   38,   38, 0x08,
     956,   38,   38,   38, 0x08,
     981,   38,   38,   38, 0x08,
    1015,   38,   38,   38, 0x08,
    1036,   38,   38,   38, 0x08,
    1057,   38,   38,   38, 0x08,
    1078,   38,   38,   38, 0x08,
    1102,   38,   38,   38, 0x08,
    1124,   38,   38,   38, 0x08,
    1146,   38,   38,   38, 0x08,
    1170,   38,   38,   38, 0x08,
    1192,   38,   38,   38, 0x08,
    1222,   38,   38,   38, 0x08,
    1247,   38,   38,   38, 0x08,
    1275,   38,   38,   38, 0x08,
    1297,   38,   38,   38, 0x08,
    1323,   38,   38,   38, 0x08,
    1345,   38,   38,   38, 0x08,
    1372,   38,   38,   38, 0x08,
    1395,   38,   38,   38, 0x08,
    1418,   38,   38,   38, 0x08,
    1437,   38,   38,   38, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_MainWindow[] = {
    "MainWindow\0on_RunSimulation_clicked()\0"
    "\0on_STLRender_clicked()\0on_STLOpen_clicked()\0"
    "on_MeshUndo_clicked()\0on_MeshUpdate_clicked()\0"
    "on_New_clicked()\0on_Save_clicked()\0"
    "on_actionOpen_triggered()\0"
    "on_actionSave_triggered()\0"
    "on_actionQuit_triggered()\0"
    "on_kEpsilon_toggled()\0on_Laminar_toggled()\0"
    "on_MeshParameters_itemDoubleClicked(QTreeWidgetItem*,int)\0"
    "item,column\0on_EditBoundary_clicked()\0"
    "on_BoundaryTree_itemDoubleClicked(QTreeWidgetItem*,int)\0"
    "on_AddSpecialBoundary_clicked()\0"
    "on_BoundaryTree_itemActivated(QTreeWidgetItem*,int)\0"
    "on_BoundaryTree_itemClicked(QTreeWidgetItem*,int)\0"
    "on_RemoveSpecialBoundary_clicked()\0"
    "on_EditBaffle_clicked()\0"
    "on_BaffleTree_itemDoubleClicked(QTreeWidgetItem*,int)\0"
    "on_AddBaffle_clicked()\0"
    "on_BaffleTree_itemActivated(QTreeWidgetItem*,int)\0"
    "on_BaffleTree_itemClicked(QTreeWidgetItem*,int)\0"
    "on_RemoveBaffle_clicked()\0"
    "on_updateRange_clicked()\0on_saveJPEG_clicked()\0"
    "on_inside_x_textChanged()\0"
    "on_inside_y_textChanged()\0"
    "on_inside_z_textChanged()\0"
    "on_origin_valueChanged()\0"
    "on_timesteps_currentItemChanged()\0"
    "on_xNormal_toggled()\0on_yNormal_toggled()\0"
    "on_zNormal_toggled()\0on_contourVOF_toggled()\0"
    "on_contourP_toggled()\0on_contourK_toggled()\0"
    "on_showLegend_toggled()\0on_showAxis_toggled()\0"
    "on_contourVorticity_toggled()\0"
    "on_showVectors_toggled()\0"
    "on_blockObstacles_toggled()\0"
    "on_showMesh_toggled()\0on_earthGravity_clicked()\0"
    "on_water20C_clicked()\0on_defaultLength_clicked()\0"
    "on_calcRough_clicked()\0on_SelectAll_clicked()\0"
    "on_Clear_clicked()\0on_Delete_clicked()\0"
};

void MainWindow::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        MainWindow *_t = static_cast<MainWindow *>(_o);
        switch (_id) {
        case 0: _t->on_RunSimulation_clicked(); break;
        case 1: _t->on_STLRender_clicked(); break;
        case 2: _t->on_STLOpen_clicked(); break;
        case 3: _t->on_MeshUndo_clicked(); break;
        case 4: _t->on_MeshUpdate_clicked(); break;
        case 5: _t->on_New_clicked(); break;
        case 6: _t->on_Save_clicked(); break;
        case 7: _t->on_actionOpen_triggered(); break;
        case 8: _t->on_actionSave_triggered(); break;
        case 9: _t->on_actionQuit_triggered(); break;
        case 10: _t->on_kEpsilon_toggled(); break;
        case 11: _t->on_Laminar_toggled(); break;
        case 12: _t->on_MeshParameters_itemDoubleClicked((*reinterpret_cast< QTreeWidgetItem*(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        case 13: _t->on_EditBoundary_clicked(); break;
        case 14: _t->on_BoundaryTree_itemDoubleClicked((*reinterpret_cast< QTreeWidgetItem*(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        case 15: _t->on_AddSpecialBoundary_clicked(); break;
        case 16: _t->on_BoundaryTree_itemActivated((*reinterpret_cast< QTreeWidgetItem*(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        case 17: _t->on_BoundaryTree_itemClicked((*reinterpret_cast< QTreeWidgetItem*(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        case 18: _t->on_RemoveSpecialBoundary_clicked(); break;
        case 19: _t->on_EditBaffle_clicked(); break;
        case 20: _t->on_BaffleTree_itemDoubleClicked((*reinterpret_cast< QTreeWidgetItem*(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        case 21: _t->on_AddBaffle_clicked(); break;
        case 22: _t->on_BaffleTree_itemActivated((*reinterpret_cast< QTreeWidgetItem*(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        case 23: _t->on_BaffleTree_itemClicked((*reinterpret_cast< QTreeWidgetItem*(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        case 24: _t->on_RemoveBaffle_clicked(); break;
        case 25: _t->on_updateRange_clicked(); break;
        case 26: _t->on_saveJPEG_clicked(); break;
        case 27: _t->on_inside_x_textChanged(); break;
        case 28: _t->on_inside_y_textChanged(); break;
        case 29: _t->on_inside_z_textChanged(); break;
        case 30: _t->on_origin_valueChanged(); break;
        case 31: _t->on_timesteps_currentItemChanged(); break;
        case 32: _t->on_xNormal_toggled(); break;
        case 33: _t->on_yNormal_toggled(); break;
        case 34: _t->on_zNormal_toggled(); break;
        case 35: _t->on_contourVOF_toggled(); break;
        case 36: _t->on_contourP_toggled(); break;
        case 37: _t->on_contourK_toggled(); break;
        case 38: _t->on_showLegend_toggled(); break;
        case 39: _t->on_showAxis_toggled(); break;
        case 40: _t->on_contourVorticity_toggled(); break;
        case 41: _t->on_showVectors_toggled(); break;
        case 42: _t->on_blockObstacles_toggled(); break;
        case 43: _t->on_showMesh_toggled(); break;
        case 44: _t->on_earthGravity_clicked(); break;
        case 45: _t->on_water20C_clicked(); break;
        case 46: _t->on_defaultLength_clicked(); break;
        case 47: _t->on_calcRough_clicked(); break;
        case 48: _t->on_SelectAll_clicked(); break;
        case 49: _t->on_Clear_clicked(); break;
        case 50: _t->on_Delete_clicked(); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData MainWindow::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject MainWindow::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_MainWindow,
      qt_meta_data_MainWindow, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &MainWindow::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *MainWindow::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *MainWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_MainWindow))
        return static_cast<void*>(const_cast< MainWindow*>(this));
    return QMainWindow::qt_metacast(_clname);
}

int MainWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 51)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 51;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
