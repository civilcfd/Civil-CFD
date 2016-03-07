/********************************************************************************
** Form generated from reading UI file 'render.ui'
**
** Created by: Qt User Interface Compiler version 4.8.7
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_RENDER_H
#define UI_RENDER_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QDialog>
#include <QtGui/QGridLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QPlainTextEdit>
#include <QtGui/QProgressBar>
#include <QtGui/QPushButton>
#include <QtGui/QSpacerItem>
#include "QVTKWidget.h"

QT_BEGIN_NAMESPACE

class Ui_Render
{
public:
    QGridLayout *gridLayout_2;
    QGridLayout *gridLayout;
    QSpacerItem *verticalSpacer_2;
    QPushButton *Return;
    QSpacerItem *verticalSpacer;
    QPushButton *Stop;
    QVTKWidget *VTKRender;
    QLabel *status;
    QSpacerItem *verticalSpacer_3;
    QProgressBar *progressBar;
    QPlainTextEdit *output;

    void setupUi(QDialog *Render)
    {
        if (Render->objectName().isEmpty())
            Render->setObjectName(QString::fromUtf8("Render"));
        Render->resize(546, 614);
        gridLayout_2 = new QGridLayout(Render);
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        gridLayout = new QGridLayout();
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        verticalSpacer_2 = new QSpacerItem(20, 10, QSizePolicy::Minimum, QSizePolicy::Maximum);

        gridLayout->addItem(verticalSpacer_2, 1, 0, 1, 1);

        Return = new QPushButton(Render);
        Return->setObjectName(QString::fromUtf8("Return"));
        Return->setEnabled(false);
        Return->setMinimumSize(QSize(150, 0));
        Return->setMaximumSize(QSize(150, 16777215));
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/icons/resources/system-log-out.svg"), QSize(), QIcon::Normal, QIcon::Off);
        Return->setIcon(icon);
        Return->setIconSize(QSize(22, 22));

        gridLayout->addWidget(Return, 4, 1, 1, 1, Qt::AlignHCenter);

        verticalSpacer = new QSpacerItem(20, 20, QSizePolicy::Minimum, QSizePolicy::Maximum);

        gridLayout->addItem(verticalSpacer, 5, 0, 1, 1);

        Stop = new QPushButton(Render);
        Stop->setObjectName(QString::fromUtf8("Stop"));
        Stop->setMinimumSize(QSize(150, 0));
        Stop->setMaximumSize(QSize(150, 16777215));
        QIcon icon1;
        icon1.addFile(QString::fromUtf8(":/icons/resources/process-stop.svg"), QSize(), QIcon::Normal, QIcon::Off);
        Stop->setIcon(icon1);
        Stop->setIconSize(QSize(22, 22));

        gridLayout->addWidget(Stop, 4, 0, 1, 1, Qt::AlignHCenter);

        VTKRender = new QVTKWidget(Render);
        VTKRender->setObjectName(QString::fromUtf8("VTKRender"));
        VTKRender->setMinimumSize(QSize(100, 250));

        gridLayout->addWidget(VTKRender, 6, 0, 1, 2);

        status = new QLabel(Render);
        status->setObjectName(QString::fromUtf8("status"));
        status->setMaximumSize(QSize(16777215, 20));
        status->setAlignment(Qt::AlignCenter);

        gridLayout->addWidget(status, 0, 0, 1, 2);

        verticalSpacer_3 = new QSpacerItem(20, 10, QSizePolicy::Minimum, QSizePolicy::Maximum);

        gridLayout->addItem(verticalSpacer_3, 3, 0, 1, 1);

        progressBar = new QProgressBar(Render);
        progressBar->setObjectName(QString::fromUtf8("progressBar"));
        progressBar->setValue(0);

        gridLayout->addWidget(progressBar, 2, 0, 1, 2);

        output = new QPlainTextEdit(Render);
        output->setObjectName(QString::fromUtf8("output"));
        QSizePolicy sizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(output->sizePolicy().hasHeightForWidth());
        output->setSizePolicy(sizePolicy);

        gridLayout->addWidget(output, 7, 0, 1, 2);


        gridLayout_2->addLayout(gridLayout, 0, 0, 1, 1);


        retranslateUi(Render);

        QMetaObject::connectSlotsByName(Render);
    } // setupUi

    void retranslateUi(QDialog *Render)
    {
        Render->setWindowTitle(QApplication::translate("Render", "Render", 0, QApplication::UnicodeUTF8));
        Return->setText(QApplication::translate("Render", "Return", 0, QApplication::UnicodeUTF8));
        Stop->setText(QApplication::translate("Render", "Stop Rendering", 0, QApplication::UnicodeUTF8));
        status->setText(QApplication::translate("Render", "Status", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class Render: public Ui_Render {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_RENDER_H
