/********************************************************************************
** Form generated from reading UI file 'sboundary.ui'
**
** Created by: Qt User Interface Compiler version 4.8.7
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_SBOUNDARY_H
#define UI_SBOUNDARY_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QComboBox>
#include <QtGui/QDialog>
#include <QtGui/QDialogButtonBox>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QPlainTextEdit>
#include <QtGui/QPushButton>

QT_BEGIN_NAMESPACE

class Ui_SBoundaryDialog
{
public:
    QDialogButtonBox *buttonBox;
    QLabel *label;
    QComboBox *select;
    QLabel *labelValue;
    QPlainTextEdit *value;
    QPlainTextEdit *turbulence;
    QLabel *labelTurbulence;
    QLabel *label_4;
    QLabel *A1;
    QPlainTextEdit *extentA1;
    QLabel *B1;
    QPlainTextEdit *extentB1;
    QPlainTextEdit *extentB2;
    QLabel *B2;
    QPlainTextEdit *extentA2;
    QLabel *A2;
    QPushButton *calculate;
    QLabel *units;
    QLabel *maxA;
    QLabel *maxB;

    void setupUi(QDialog *SBoundaryDialog)
    {
        if (SBoundaryDialog->objectName().isEmpty())
            SBoundaryDialog->setObjectName(QString::fromUtf8("SBoundaryDialog"));
        SBoundaryDialog->resize(500, 424);
        QSizePolicy sizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(SBoundaryDialog->sizePolicy().hasHeightForWidth());
        SBoundaryDialog->setSizePolicy(sizePolicy);
        SBoundaryDialog->setMinimumSize(QSize(500, 424));
        SBoundaryDialog->setMaximumSize(QSize(500, 424));
        buttonBox = new QDialogButtonBox(SBoundaryDialog);
        buttonBox->setObjectName(QString::fromUtf8("buttonBox"));
        buttonBox->setGeometry(QRect(40, 380, 291, 32));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);
        label = new QLabel(SBoundaryDialog);
        label->setObjectName(QString::fromUtf8("label"));
        label->setGeometry(QRect(20, 170, 251, 16));
        select = new QComboBox(SBoundaryDialog);
        select->setObjectName(QString::fromUtf8("select"));
        select->setGeometry(QRect(70, 200, 301, 31));
        labelValue = new QLabel(SBoundaryDialog);
        labelValue->setObjectName(QString::fromUtf8("labelValue"));
        labelValue->setGeometry(QRect(30, 270, 251, 16));
        value = new QPlainTextEdit(SBoundaryDialog);
        value->setObjectName(QString::fromUtf8("value"));
        value->setGeometry(QRect(160, 260, 161, 31));
        value->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        value->setTabChangesFocus(true);
        value->setLineWrapMode(QPlainTextEdit::NoWrap);
        turbulence = new QPlainTextEdit(SBoundaryDialog);
        turbulence->setObjectName(QString::fromUtf8("turbulence"));
        turbulence->setGeometry(QRect(160, 320, 161, 31));
        turbulence->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        turbulence->setTabChangesFocus(true);
        turbulence->setLineWrapMode(QPlainTextEdit::NoWrap);
        labelTurbulence = new QLabel(SBoundaryDialog);
        labelTurbulence->setObjectName(QString::fromUtf8("labelTurbulence"));
        labelTurbulence->setGeometry(QRect(30, 330, 71, 16));
        label_4 = new QLabel(SBoundaryDialog);
        label_4->setObjectName(QString::fromUtf8("label_4"));
        label_4->setGeometry(QRect(20, 20, 291, 16));
        A1 = new QLabel(SBoundaryDialog);
        A1->setObjectName(QString::fromUtf8("A1"));
        A1->setGeometry(QRect(30, 60, 51, 16));
        extentA1 = new QPlainTextEdit(SBoundaryDialog);
        extentA1->setObjectName(QString::fromUtf8("extentA1"));
        extentA1->setGeometry(QRect(90, 50, 81, 31));
        extentA1->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        extentA1->setTabChangesFocus(true);
        extentA1->setLineWrapMode(QPlainTextEdit::NoWrap);
        B1 = new QLabel(SBoundaryDialog);
        B1->setObjectName(QString::fromUtf8("B1"));
        B1->setGeometry(QRect(260, 60, 41, 16));
        extentB1 = new QPlainTextEdit(SBoundaryDialog);
        extentB1->setObjectName(QString::fromUtf8("extentB1"));
        extentB1->setGeometry(QRect(310, 50, 81, 31));
        extentB1->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        extentB1->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        extentB1->setTabChangesFocus(true);
        extentB1->setLineWrapMode(QPlainTextEdit::NoWrap);
        extentB2 = new QPlainTextEdit(SBoundaryDialog);
        extentB2->setObjectName(QString::fromUtf8("extentB2"));
        extentB2->setGeometry(QRect(310, 100, 81, 31));
        extentB2->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        extentB2->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        extentB2->setTabChangesFocus(true);
        extentB2->setLineWrapMode(QPlainTextEdit::NoWrap);
        B2 = new QLabel(SBoundaryDialog);
        B2->setObjectName(QString::fromUtf8("B2"));
        B2->setGeometry(QRect(260, 110, 41, 16));
        extentA2 = new QPlainTextEdit(SBoundaryDialog);
        extentA2->setObjectName(QString::fromUtf8("extentA2"));
        extentA2->setGeometry(QRect(90, 100, 81, 31));
        extentA2->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        extentA2->setTabChangesFocus(true);
        extentA2->setLineWrapMode(QPlainTextEdit::NoWrap);
        A2 = new QLabel(SBoundaryDialog);
        A2->setObjectName(QString::fromUtf8("A2"));
        A2->setGeometry(QRect(30, 110, 51, 16));
        calculate = new QPushButton(SBoundaryDialog);
        calculate->setObjectName(QString::fromUtf8("calculate"));
        calculate->setGeometry(QRect(330, 320, 110, 32));
        units = new QLabel(SBoundaryDialog);
        units->setObjectName(QString::fromUtf8("units"));
        units->setGeometry(QRect(340, 270, 59, 16));
        maxA = new QLabel(SBoundaryDialog);
        maxA->setObjectName(QString::fromUtf8("maxA"));
        maxA->setGeometry(QRect(400, 60, 56, 13));
        maxB = new QLabel(SBoundaryDialog);
        maxB->setObjectName(QString::fromUtf8("maxB"));
        maxB->setGeometry(QRect(400, 110, 56, 13));

        retranslateUi(SBoundaryDialog);
        QObject::connect(buttonBox, SIGNAL(accepted()), SBoundaryDialog, SLOT(accept()));
        QObject::connect(buttonBox, SIGNAL(rejected()), SBoundaryDialog, SLOT(reject()));

        QMetaObject::connectSlotsByName(SBoundaryDialog);
    } // setupUi

    void retranslateUi(QDialog *SBoundaryDialog)
    {
        SBoundaryDialog->setWindowTitle(QApplication::translate("SBoundaryDialog", "Special Boundary", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("SBoundaryDialog", "Select Boundary Type", 0, QApplication::UnicodeUTF8));
        select->clear();
        select->insertItems(0, QStringList()
         << QApplication::translate("SBoundaryDialog", "fixed velocity", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("SBoundaryDialog", "mass outflow", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("SBoundaryDialog", "hydraulic grade", 0, QApplication::UnicodeUTF8)
        );
        labelValue->setText(QApplication::translate("SBoundaryDialog", "Set Value", 0, QApplication::UnicodeUTF8));
        value->setPlainText(QApplication::translate("SBoundaryDialog", "0", 0, QApplication::UnicodeUTF8));
        turbulence->setPlainText(QApplication::translate("SBoundaryDialog", "0.001", 0, QApplication::UnicodeUTF8));
        labelTurbulence->setText(QApplication::translate("SBoundaryDialog", "Turbulence", 0, QApplication::UnicodeUTF8));
        label_4->setText(QApplication::translate("SBoundaryDialog", "Boundary Extents - mesh cell range", 0, QApplication::UnicodeUTF8));
        A1->setText(QApplication::translate("SBoundaryDialog", "A1", 0, QApplication::UnicodeUTF8));
        extentA1->setPlainText(QApplication::translate("SBoundaryDialog", "0", 0, QApplication::UnicodeUTF8));
        B1->setText(QApplication::translate("SBoundaryDialog", "B1", 0, QApplication::UnicodeUTF8));
        B2->setText(QApplication::translate("SBoundaryDialog", "B2", 0, QApplication::UnicodeUTF8));
        extentA2->setPlainText(QApplication::translate("SBoundaryDialog", "0", 0, QApplication::UnicodeUTF8));
        A2->setText(QApplication::translate("SBoundaryDialog", "A2", 0, QApplication::UnicodeUTF8));
        calculate->setText(QApplication::translate("SBoundaryDialog", "Calculate", 0, QApplication::UnicodeUTF8));
        units->setText(QApplication::translate("SBoundaryDialog", "m/s", 0, QApplication::UnicodeUTF8));
        maxA->setText(QApplication::translate("SBoundaryDialog", "of ", 0, QApplication::UnicodeUTF8));
        maxB->setText(QApplication::translate("SBoundaryDialog", "of ", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class SBoundaryDialog: public Ui_SBoundaryDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_SBOUNDARY_H
