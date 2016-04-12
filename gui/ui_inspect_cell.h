/********************************************************************************
** Form generated from reading UI file 'inspect_cell.ui'
**
** Created by: Qt User Interface Compiler version 4.8.7
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_INSPECT_CELL_H
#define UI_INSPECT_CELL_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QDialog>
#include <QtGui/QGridLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QPlainTextEdit>
#include <QtGui/QPushButton>
#include <QtGui/QSpacerItem>

QT_BEGIN_NAMESPACE

class Ui_InspectDialog
{
public:
    QGridLayout *gridLayout_2;
    QGridLayout *gridLayout;
    QPlainTextEdit *output;
    QGridLayout *gridLayout_5;
    QLabel *label_4;
    QPlainTextEdit *z;
    QLabel *maxZ;
    QGridLayout *gridLayout_3;
    QPlainTextEdit *x;
    QLabel *label;
    QLabel *maxX;
    QPushButton *Return;
    QSpacerItem *verticalSpacer;
    QPushButton *InspectCell;
    QLabel *status;
    QGridLayout *gridLayout_4;
    QPlainTextEdit *y;
    QLabel *label_2;
    QLabel *maxY;
    QSpacerItem *verticalSpacer_2;
    QSpacerItem *verticalSpacer_3;
    QGridLayout *gridLayout_6;
    QPlainTextEdit *timestep;
    QLabel *label_3;
    QLabel *label_5;

    void setupUi(QDialog *InspectDialog)
    {
        if (InspectDialog->objectName().isEmpty())
            InspectDialog->setObjectName(QString::fromUtf8("InspectDialog"));
        InspectDialog->resize(546, 500);
        gridLayout_2 = new QGridLayout(InspectDialog);
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        gridLayout = new QGridLayout();
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        output = new QPlainTextEdit(InspectDialog);
        output->setObjectName(QString::fromUtf8("output"));
        QSizePolicy sizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(output->sizePolicy().hasHeightForWidth());
        output->setSizePolicy(sizePolicy);
        output->setMinimumSize(QSize(0, 220));
        output->setReadOnly(true);

        gridLayout->addWidget(output, 16, 0, 1, 2);

        gridLayout_5 = new QGridLayout();
        gridLayout_5->setObjectName(QString::fromUtf8("gridLayout_5"));
        label_4 = new QLabel(InspectDialog);
        label_4->setObjectName(QString::fromUtf8("label_4"));

        gridLayout_5->addWidget(label_4, 0, 0, 1, 1, Qt::AlignRight);

        z = new QPlainTextEdit(InspectDialog);
        z->setObjectName(QString::fromUtf8("z"));
        z->setMaximumSize(QSize(80, 31));
        z->setTabChangesFocus(true);
        z->setLineWrapMode(QPlainTextEdit::NoWrap);

        gridLayout_5->addWidget(z, 0, 1, 1, 1);

        maxZ = new QLabel(InspectDialog);
        maxZ->setObjectName(QString::fromUtf8("maxZ"));

        gridLayout_5->addWidget(maxZ, 0, 2, 1, 1);


        gridLayout->addLayout(gridLayout_5, 9, 0, 1, 2);

        gridLayout_3 = new QGridLayout();
        gridLayout_3->setObjectName(QString::fromUtf8("gridLayout_3"));
        gridLayout_3->setSizeConstraint(QLayout::SetMaximumSize);
        x = new QPlainTextEdit(InspectDialog);
        x->setObjectName(QString::fromUtf8("x"));
        x->setMaximumSize(QSize(81, 31));
        x->setTabChangesFocus(true);
        x->setLineWrapMode(QPlainTextEdit::NoWrap);

        gridLayout_3->addWidget(x, 0, 1, 1, 1);

        label = new QLabel(InspectDialog);
        label->setObjectName(QString::fromUtf8("label"));
        label->setMaximumSize(QSize(16777215, 31));

        gridLayout_3->addWidget(label, 0, 0, 1, 1, Qt::AlignRight);

        maxX = new QLabel(InspectDialog);
        maxX->setObjectName(QString::fromUtf8("maxX"));

        gridLayout_3->addWidget(maxX, 0, 2, 1, 1);


        gridLayout->addLayout(gridLayout_3, 5, 0, 1, 2);

        Return = new QPushButton(InspectDialog);
        Return->setObjectName(QString::fromUtf8("Return"));
        Return->setEnabled(false);
        Return->setMinimumSize(QSize(150, 0));
        Return->setMaximumSize(QSize(150, 16777215));
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/icons/resources/system-log-out.svg"), QSize(), QIcon::Normal, QIcon::Off);
        Return->setIcon(icon);
        Return->setIconSize(QSize(22, 22));

        gridLayout->addWidget(Return, 14, 1, 1, 1, Qt::AlignHCenter);

        verticalSpacer = new QSpacerItem(20, 20, QSizePolicy::Minimum, QSizePolicy::Maximum);

        gridLayout->addItem(verticalSpacer, 15, 0, 1, 1);

        InspectCell = new QPushButton(InspectDialog);
        InspectCell->setObjectName(QString::fromUtf8("InspectCell"));
        InspectCell->setMinimumSize(QSize(150, 0));
        InspectCell->setMaximumSize(QSize(150, 16777215));
        QIcon icon1;
        icon1.addFile(QString::fromUtf8(":/icons/resources/system-search.svg"), QSize(), QIcon::Normal, QIcon::Off);
        InspectCell->setIcon(icon1);
        InspectCell->setIconSize(QSize(22, 22));

        gridLayout->addWidget(InspectCell, 14, 0, 1, 1, Qt::AlignHCenter);

        status = new QLabel(InspectDialog);
        status->setObjectName(QString::fromUtf8("status"));
        status->setMaximumSize(QSize(16777215, 20));
        QFont font;
        font.setBold(true);
        font.setWeight(75);
        status->setFont(font);
        status->setAlignment(Qt::AlignCenter);

        gridLayout->addWidget(status, 0, 0, 1, 2, Qt::AlignLeft);

        gridLayout_4 = new QGridLayout();
        gridLayout_4->setObjectName(QString::fromUtf8("gridLayout_4"));
        y = new QPlainTextEdit(InspectDialog);
        y->setObjectName(QString::fromUtf8("y"));
        y->setMaximumSize(QSize(80, 31));
        y->setTabChangesFocus(true);
        y->setLineWrapMode(QPlainTextEdit::NoWrap);

        gridLayout_4->addWidget(y, 0, 1, 1, 1);

        label_2 = new QLabel(InspectDialog);
        label_2->setObjectName(QString::fromUtf8("label_2"));

        gridLayout_4->addWidget(label_2, 0, 0, 1, 1, Qt::AlignRight);

        maxY = new QLabel(InspectDialog);
        maxY->setObjectName(QString::fromUtf8("maxY"));

        gridLayout_4->addWidget(maxY, 0, 2, 1, 1);


        gridLayout->addLayout(gridLayout_4, 6, 0, 1, 2);

        verticalSpacer_2 = new QSpacerItem(20, 10, QSizePolicy::Minimum, QSizePolicy::Maximum);

        gridLayout->addItem(verticalSpacer_2, 1, 0, 1, 1);

        verticalSpacer_3 = new QSpacerItem(20, 10, QSizePolicy::Minimum, QSizePolicy::Maximum);

        gridLayout->addItem(verticalSpacer_3, 13, 0, 1, 1);

        gridLayout_6 = new QGridLayout();
        gridLayout_6->setObjectName(QString::fromUtf8("gridLayout_6"));
        timestep = new QPlainTextEdit(InspectDialog);
        timestep->setObjectName(QString::fromUtf8("timestep"));
        timestep->setMaximumSize(QSize(81, 32));
        timestep->setReadOnly(true);

        gridLayout_6->addWidget(timestep, 0, 1, 1, 1);

        label_3 = new QLabel(InspectDialog);
        label_3->setObjectName(QString::fromUtf8("label_3"));

        gridLayout_6->addWidget(label_3, 0, 0, 1, 1, Qt::AlignRight);

        label_5 = new QLabel(InspectDialog);
        label_5->setObjectName(QString::fromUtf8("label_5"));

        gridLayout_6->addWidget(label_5, 0, 2, 1, 1);


        gridLayout->addLayout(gridLayout_6, 3, 0, 1, 2);


        gridLayout_2->addLayout(gridLayout, 0, 0, 1, 1);

        QWidget::setTabOrder(x, y);
        QWidget::setTabOrder(y, z);
        QWidget::setTabOrder(z, InspectCell);
        QWidget::setTabOrder(InspectCell, Return);
        QWidget::setTabOrder(Return, output);

        retranslateUi(InspectDialog);

        QMetaObject::connectSlotsByName(InspectDialog);
    } // setupUi

    void retranslateUi(QDialog *InspectDialog)
    {
        InspectDialog->setWindowTitle(QApplication::translate("InspectDialog", "Inspect Cell", 0, QApplication::UnicodeUTF8));
        label_4->setText(QApplication::translate("InspectDialog", "z:", 0, QApplication::UnicodeUTF8));
        maxZ->setText(QApplication::translate("InspectDialog", "of", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("InspectDialog", "x:", 0, QApplication::UnicodeUTF8));
        maxX->setText(QApplication::translate("InspectDialog", "of", 0, QApplication::UnicodeUTF8));
        Return->setText(QApplication::translate("InspectDialog", "Return", 0, QApplication::UnicodeUTF8));
        InspectCell->setText(QApplication::translate("InspectDialog", "Inspect Cell", 0, QApplication::UnicodeUTF8));
        status->setText(QApplication::translate("InspectDialog", "Inspect Cell", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("InspectDialog", "y:", 0, QApplication::UnicodeUTF8));
        maxY->setText(QApplication::translate("InspectDialog", "of", 0, QApplication::UnicodeUTF8));
        label_3->setText(QApplication::translate("InspectDialog", "timestep:", 0, QApplication::UnicodeUTF8));
        label_5->setText(QString());
    } // retranslateUi

};

namespace Ui {
    class InspectDialog: public Ui_InspectDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_INSPECT_CELL_H
