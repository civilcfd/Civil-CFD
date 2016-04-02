/********************************************************************************
** Form generated from reading UI file 'about.ui'
**
** Created by: Qt User Interface Compiler version 4.8.6
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_ABOUT_H
#define UI_ABOUT_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QDialog>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QPushButton>

QT_BEGIN_NAMESPACE

class Ui_AboutDialog
{
public:
    QLabel *label;
    QLabel *label_2;
    QPushButton *OK;

    void setupUi(QDialog *AboutDialog)
    {
        if (AboutDialog->objectName().isEmpty())
            AboutDialog->setObjectName(QString::fromUtf8("AboutDialog"));
        AboutDialog->resize(740, 343);
        AboutDialog->setModal(true);
        label = new QLabel(AboutDialog);
        label->setObjectName(QString::fromUtf8("label"));
        label->setGeometry(QRect(20, 20, 51, 51));
        label->setPixmap(QPixmap(QString::fromUtf8(":/icons/icon.svg")));
        label_2 = new QLabel(AboutDialog);
        label_2->setObjectName(QString::fromUtf8("label_2"));
        label_2->setGeometry(QRect(30, 10, 681, 271));
        label_2->setTextFormat(Qt::RichText);
        label_2->setWordWrap(true);
        OK = new QPushButton(AboutDialog);
        OK->setObjectName(QString::fromUtf8("OK"));
        OK->setGeometry(QRect(610, 310, 110, 32));
        label_2->raise();
        label->raise();
        OK->raise();

        retranslateUi(AboutDialog);

        QMetaObject::connectSlotsByName(AboutDialog);
    } // setupUi

    void retranslateUi(QDialog *AboutDialog)
    {
        AboutDialog->setWindowTitle(QApplication::translate("AboutDialog", "About Civil CFD", 0, QApplication::UnicodeUTF8));
        label->setText(QString());
        label_2->setText(QApplication::translate("AboutDialog", "<html><head/><body><p align=\"center\"><span style=\" font-weight:600;\">Civil CFD v. 0.1<br/>Free Surface CFD Solver for Civil Engineering Problems</span></p><p align=\"center\">Copyright (c) 2016 by Michael Celli<br/>Released under the GNU LGPL (<a href=\"www.gnu.org\"><span style=\" text-decoration: underline; color:#0000ff;\">www.gnu.org</span></a>)<br/><br/>VTK is Copyright (c) 1993-2016 Ken Martin, Will Schroeder, Bill Lorensen<br/>The Qt GUI Toolkit is Copyright (c) 2016 The Qt Company Ltd.<br/>Qhull is Copyright (c) 1993-2016 C.B. Barber and the The Geometry Centre<br/>QCustomPlot is Copyright (c) 2016 Emanuel Eichhammer</p><p align=\"center\">The author(s) of Civil CFD provide the program '<span style=\" font-weight:600;\">as is'</span> without warranty of any kind, expressed or implied, including, but not limited to, the <span style=\" font-weight:600;\">implied warranties of merchantability and fitness for a particular purpose</span>.</p><p align=\"center\">Civil CFD is an open source framework for "
                        "numerical simulation of differential equations using finite difference methods. The methods, equations, and theory are outlined in the source files that are made available to the user for study or modification. Practicioners using this software must apply their own professional judgement when interpreting the results, and take sole responsibility for professional decisions made based on use of this software. The author(s) of this software shall not be liable for professional decisions made based on results obtained from this software, including results that may be erroneous due to errors in the code. By using this software, you agree not to hold the author(s) liable for any injury or cost to any party from using this software or relying on its results.</p></body></html>", 0, QApplication::UnicodeUTF8));
        OK->setText(QApplication::translate("AboutDialog", "OK", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class AboutDialog: public Ui_AboutDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_ABOUT_H
