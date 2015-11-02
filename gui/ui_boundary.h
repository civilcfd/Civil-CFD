/********************************************************************************
** Form generated from reading UI file 'boundary.ui'
**
** Created by: Qt User Interface Compiler version 4.8.7
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_BOUNDARY_H
#define UI_BOUNDARY_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QComboBox>
#include <QtGui/QDialog>
#include <QtGui/QDialogButtonBox>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>

QT_BEGIN_NAMESPACE

class Ui_BoundaryDialog
{
public:
    QDialogButtonBox *buttonBox;
    QComboBox *select;
    QLabel *label;

    void setupUi(QDialog *BoundaryDialog)
    {
        if (BoundaryDialog->objectName().isEmpty())
            BoundaryDialog->setObjectName(QString::fromUtf8("BoundaryDialog"));
        BoundaryDialog->resize(330, 160);
        buttonBox = new QDialogButtonBox(BoundaryDialog);
        buttonBox->setObjectName(QString::fromUtf8("buttonBox"));
        buttonBox->setGeometry(QRect(30, 120, 281, 32));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);
        select = new QComboBox(BoundaryDialog);
        select->setObjectName(QString::fromUtf8("select"));
        select->setGeometry(QRect(20, 60, 291, 31));
        label = new QLabel(BoundaryDialog);
        label->setObjectName(QString::fromUtf8("label"));
        label->setGeometry(QRect(20, 20, 251, 16));

        retranslateUi(BoundaryDialog);
        QObject::connect(buttonBox, SIGNAL(accepted()), BoundaryDialog, SLOT(accept()));
        QObject::connect(buttonBox, SIGNAL(rejected()), BoundaryDialog, SLOT(reject()));

        QMetaObject::connectSlotsByName(BoundaryDialog);
    } // setupUi

    void retranslateUi(QDialog *BoundaryDialog)
    {
        BoundaryDialog->setWindowTitle(QApplication::translate("BoundaryDialog", "Select Boundary", 0, QApplication::UnicodeUTF8));
        select->clear();
        select->insertItems(0, QStringList()
         << QApplication::translate("BoundaryDialog", "slip", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("BoundaryDialog", "no slip", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("BoundaryDialog", "continuous", 0, QApplication::UnicodeUTF8)
        );
        label->setText(QApplication::translate("BoundaryDialog", "Select Boundary Type", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class BoundaryDialog: public Ui_BoundaryDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_BOUNDARY_H
