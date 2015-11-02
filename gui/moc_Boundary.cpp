/****************************************************************************
** Meta object code from reading C++ file 'Boundary.h'
**
** Created by: The Qt Meta Object Compiler version 63 (Qt 4.8.7)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "Boundary.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'Boundary.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.7. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_BoundaryDialog[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       1,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      15,   24,   24,   24, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_BoundaryDialog[] = {
    "BoundaryDialog\0accept()\0\0"
};

void BoundaryDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        BoundaryDialog *_t = static_cast<BoundaryDialog *>(_o);
        switch (_id) {
        case 0: _t->accept(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObjectExtraData BoundaryDialog::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject BoundaryDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_BoundaryDialog,
      qt_meta_data_BoundaryDialog, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &BoundaryDialog::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *BoundaryDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *BoundaryDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_BoundaryDialog))
        return static_cast<void*>(const_cast< BoundaryDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int BoundaryDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 1)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 1;
    }
    return _id;
}
static const uint qt_meta_data_SBoundaryDialog[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       2,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      16,   25,   25,   25, 0x08,
      26,   25,   25,   25, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_SBoundaryDialog[] = {
    "SBoundaryDialog\0accept()\0\0"
    "on_calculate_clicked()\0"
};

void SBoundaryDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        SBoundaryDialog *_t = static_cast<SBoundaryDialog *>(_o);
        switch (_id) {
        case 0: _t->accept(); break;
        case 1: _t->on_calculate_clicked(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObjectExtraData SBoundaryDialog::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject SBoundaryDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_SBoundaryDialog,
      qt_meta_data_SBoundaryDialog, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &SBoundaryDialog::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *SBoundaryDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *SBoundaryDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_SBoundaryDialog))
        return static_cast<void*>(const_cast< SBoundaryDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int SBoundaryDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 2)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 2;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
