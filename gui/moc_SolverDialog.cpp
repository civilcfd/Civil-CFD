/****************************************************************************
** Meta object code from reading C++ file 'SolverDialog.h'
**
** Created by: The Qt Meta Object Compiler version 63 (Qt 4.8.7)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "SolverDialog.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'SolverDialog.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.7. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_SolverDialog[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       6,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      20,   14,   13,   13, 0x08,
      66,   50,   13,   13, 0x08,
     107,  101,   13,   13, 0x08,
     144,   13,   13,   13, 0x08,
     170,   13,   13,   13, 0x08,
     190,   13,   13,   13, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_SolverDialog[] = {
    "SolverDialog\0\0error\0error(QProcess::ProcessError)\0"
    "exitCode,status\0finished(int,QProcess::ExitStatus)\0"
    "state\0stateChanged(QProcess::ProcessState)\0"
    "readyReadStandardOutput()\0on_Return_clicked()\0"
    "on_Stop_clicked()\0"
};

void SolverDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        SolverDialog *_t = static_cast<SolverDialog *>(_o);
        switch (_id) {
        case 0: _t->error((*reinterpret_cast< QProcess::ProcessError(*)>(_a[1]))); break;
        case 1: _t->finished((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< QProcess::ExitStatus(*)>(_a[2]))); break;
        case 2: _t->stateChanged((*reinterpret_cast< QProcess::ProcessState(*)>(_a[1]))); break;
        case 3: _t->readyReadStandardOutput(); break;
        case 4: _t->on_Return_clicked(); break;
        case 5: _t->on_Stop_clicked(); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData SolverDialog::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject SolverDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_SolverDialog,
      qt_meta_data_SolverDialog, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &SolverDialog::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *SolverDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *SolverDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_SolverDialog))
        return static_cast<void*>(const_cast< SolverDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int SolverDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 6)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 6;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
