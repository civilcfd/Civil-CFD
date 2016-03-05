/********************************************************************************
** Form generated from reading UI file 'baffle.ui'
**
** Created by: Qt User Interface Compiler version 4.8.7
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_BAFFLE_H
#define UI_BAFFLE_H

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

QT_BEGIN_NAMESPACE

class Ui_BaffleDialog
{
public:
    QDialogButtonBox *buttonBox;
    QLabel *label;
    QComboBox *select;
    QLabel *kfactor;
    QPlainTextEdit *value;
    QLabel *label_4;
    QLabel *A1;
    QPlainTextEdit *extentA1;
    QLabel *B1;
    QPlainTextEdit *extentB1;
    QPlainTextEdit *extentB2;
    QLabel *B2;
    QPlainTextEdit *extentA2;
    QLabel *A2;
    QLabel *label_pos;
    QPlainTextEdit *pos;

    void setupUi(QDialog *BaffleDialog)
    {
        if (BaffleDialog->objectName().isEmpty())
            BaffleDialog->setObjectName(QString::fromUtf8("BaffleDialog"));
        BaffleDialog->resize(450, 424);
        QSizePolicy sizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(BaffleDialog->sizePolicy().hasHeightForWidth());
        BaffleDialog->setSizePolicy(sizePolicy);
        BaffleDialog->setMinimumSize(QSize(450, 424));
        BaffleDialog->setMaximumSize(QSize(450, 424));
        buttonBox = new QDialogButtonBox(BaffleDialog);
        buttonBox->setObjectName(QString::fromUtf8("buttonBox"));
        buttonBox->setGeometry(QRect(40, 380, 291, 32));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);
        label = new QLabel(BaffleDialog);
        label->setObjectName(QString::fromUtf8("label"));
        label->setGeometry(QRect(20, 230, 251, 16));
        select = new QComboBox(BaffleDialog);
        select->setObjectName(QString::fromUtf8("select"));
        select->setGeometry(QRect(70, 260, 301, 31));
        kfactor = new QLabel(BaffleDialog);
        kfactor->setObjectName(QString::fromUtf8("kfactor"));
        kfactor->setGeometry(QRect(80, 330, 201, 16));
        value = new QPlainTextEdit(BaffleDialog);
        value->setObjectName(QString::fromUtf8("value"));
        value->setGeometry(QRect(160, 320, 161, 31));
        value->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        label_4 = new QLabel(BaffleDialog);
        label_4->setObjectName(QString::fromUtf8("label_4"));
        label_4->setGeometry(QRect(20, 20, 161, 16));
        A1 = new QLabel(BaffleDialog);
        A1->setObjectName(QString::fromUtf8("A1"));
        A1->setGeometry(QRect(30, 60, 51, 16));
        extentA1 = new QPlainTextEdit(BaffleDialog);
        extentA1->setObjectName(QString::fromUtf8("extentA1"));
        extentA1->setGeometry(QRect(90, 50, 81, 31));
        extentA1->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        B1 = new QLabel(BaffleDialog);
        B1->setObjectName(QString::fromUtf8("B1"));
        B1->setGeometry(QRect(260, 60, 41, 16));
        extentB1 = new QPlainTextEdit(BaffleDialog);
        extentB1->setObjectName(QString::fromUtf8("extentB1"));
        extentB1->setGeometry(QRect(310, 50, 81, 31));
        extentB1->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        extentB1->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        extentB2 = new QPlainTextEdit(BaffleDialog);
        extentB2->setObjectName(QString::fromUtf8("extentB2"));
        extentB2->setGeometry(QRect(310, 100, 81, 31));
        extentB2->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        extentB2->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        B2 = new QLabel(BaffleDialog);
        B2->setObjectName(QString::fromUtf8("B2"));
        B2->setGeometry(QRect(260, 110, 41, 16));
        extentA2 = new QPlainTextEdit(BaffleDialog);
        extentA2->setObjectName(QString::fromUtf8("extentA2"));
        extentA2->setGeometry(QRect(90, 100, 81, 31));
        extentA2->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        A2 = new QLabel(BaffleDialog);
        A2->setObjectName(QString::fromUtf8("A2"));
        A2->setGeometry(QRect(30, 110, 51, 16));
        label_pos = new QLabel(BaffleDialog);
        label_pos->setObjectName(QString::fromUtf8("label_pos"));
        label_pos->setGeometry(QRect(70, 170, 121, 16));
        pos = new QPlainTextEdit(BaffleDialog);
        pos->setObjectName(QString::fromUtf8("pos"));
        pos->setGeometry(QRect(200, 160, 81, 31));
        pos->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);

        retranslateUi(BaffleDialog);
        QObject::connect(buttonBox, SIGNAL(accepted()), BaffleDialog, SLOT(accept()));
        QObject::connect(buttonBox, SIGNAL(rejected()), BaffleDialog, SLOT(reject()));

        QMetaObject::connectSlotsByName(BaffleDialog);
    } // setupUi

    void retranslateUi(QDialog *BaffleDialog)
    {
        BaffleDialog->setWindowTitle(QApplication::translate("BaffleDialog", "Special Boundary", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("BaffleDialog", "Select Boundary Type", 0, QApplication::UnicodeUTF8));
        select->clear();
        select->insertItems(0, QStringList()
         << QApplication::translate("BaffleDialog", "flow measurement", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("BaffleDialog", "barrier", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("BaffleDialog", "headloss", 0, QApplication::UnicodeUTF8)
        );
        kfactor->setText(QApplication::translate("BaffleDialog", "K-factor", 0, QApplication::UnicodeUTF8));
        value->setPlainText(QApplication::translate("BaffleDialog", "0", 0, QApplication::UnicodeUTF8));
        label_4->setText(QApplication::translate("BaffleDialog", "Baffle Extents", 0, QApplication::UnicodeUTF8));
        A1->setText(QApplication::translate("BaffleDialog", "A1", 0, QApplication::UnicodeUTF8));
        extentA1->setPlainText(QApplication::translate("BaffleDialog", "0", 0, QApplication::UnicodeUTF8));
        B1->setText(QApplication::translate("BaffleDialog", "B1", 0, QApplication::UnicodeUTF8));
        B2->setText(QApplication::translate("BaffleDialog", "B2", 0, QApplication::UnicodeUTF8));
        extentA2->setPlainText(QApplication::translate("BaffleDialog", "0", 0, QApplication::UnicodeUTF8));
        A2->setText(QApplication::translate("BaffleDialog", "A2", 0, QApplication::UnicodeUTF8));
        label_pos->setText(QApplication::translate("BaffleDialog", "Position along axis", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class BaffleDialog: public Ui_BaffleDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_BAFFLE_H
