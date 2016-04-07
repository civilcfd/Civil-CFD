/********************************************************************************
** Form generated from reading UI file 'inspect.ui'
**
** Created by: Qt User Interface Compiler version 4.8.7
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_INSPECT_H
#define UI_INSPECT_H

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

class Ui_Inspect
{
public:
    QGridLayout *gridLayout_2;
    QGridLayout *gridLayout;
    QGridLayout *gridLayout_4;
    QPlainTextEdit *y;
    QLabel *label_2;
    QLabel *maxY;
    QSpacerItem *verticalSpacer_2;
    QPushButton *Return;
    QSpacerItem *verticalSpacer;
    QPushButton *InspectCell;
    QLabel *status;
    QSpacerItem *verticalSpacer_3;
    QPlainTextEdit *output;
    QGridLayout *gridLayout_3;
    QPlainTextEdit *x;
    QLabel *label;
    QLabel *maxX;
    QGridLayout *gridLayout_5;
    QLabel *label_4;
    QPlainTextEdit *z;
    QLabel *maxZ;

    void setupUi(QDialog *Inspect)
    {
        if (Inspect->objectName().isEmpty())
            Inspect->setObjectName(QString::fromUtf8("Inspect"));
        Inspect->resize(546, 500);
        gridLayout_2 = new QGridLayout(Inspect);
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        gridLayout = new QGridLayout();
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        gridLayout_4 = new QGridLayout();
        gridLayout_4->setObjectName(QString::fromUtf8("gridLayout_4"));
        y = new QPlainTextEdit(Inspect);
        y->setObjectName(QString::fromUtf8("y"));
        y->setMaximumSize(QSize(80, 31));
        y->setTabChangesFocus(true);
        y->setLineWrapMode(QPlainTextEdit::NoWrap);

        gridLayout_4->addWidget(y, 0, 1, 1, 1);

        label_2 = new QLabel(Inspect);
        label_2->setObjectName(QString::fromUtf8("label_2"));

        gridLayout_4->addWidget(label_2, 0, 0, 1, 1, Qt::AlignRight);

        maxY = new QLabel(Inspect);
        maxY->setObjectName(QString::fromUtf8("maxY"));

        gridLayout_4->addWidget(maxY, 0, 2, 1, 1);


        gridLayout->addLayout(gridLayout_4, 4, 0, 1, 2);

        verticalSpacer_2 = new QSpacerItem(20, 10, QSizePolicy::Minimum, QSizePolicy::Maximum);

        gridLayout->addItem(verticalSpacer_2, 1, 0, 1, 1);

        Return = new QPushButton(Inspect);
        Return->setObjectName(QString::fromUtf8("Return"));
        Return->setEnabled(false);
        Return->setMinimumSize(QSize(150, 0));
        Return->setMaximumSize(QSize(150, 16777215));
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/icons/resources/system-log-out.svg"), QSize(), QIcon::Normal, QIcon::Off);
        Return->setIcon(icon);
        Return->setIconSize(QSize(22, 22));

        gridLayout->addWidget(Return, 12, 1, 1, 1, Qt::AlignHCenter);

        verticalSpacer = new QSpacerItem(20, 20, QSizePolicy::Minimum, QSizePolicy::Maximum);

        gridLayout->addItem(verticalSpacer, 13, 0, 1, 1);

        InspectCell = new QPushButton(Inspect);
        InspectCell->setObjectName(QString::fromUtf8("InspectCell"));
        InspectCell->setMinimumSize(QSize(150, 0));
        InspectCell->setMaximumSize(QSize(150, 16777215));
        QIcon icon1;
        icon1.addFile(QString::fromUtf8(":/icons/resources/system-search.svg"), QSize(), QIcon::Normal, QIcon::Off);
        InspectCell->setIcon(icon1);
        InspectCell->setIconSize(QSize(22, 22));

        gridLayout->addWidget(InspectCell, 12, 0, 1, 1, Qt::AlignHCenter);

        status = new QLabel(Inspect);
        status->setObjectName(QString::fromUtf8("status"));
        status->setMaximumSize(QSize(16777215, 20));
        QFont font;
        font.setBold(true);
        font.setWeight(75);
        status->setFont(font);
        status->setAlignment(Qt::AlignCenter);

        gridLayout->addWidget(status, 0, 0, 1, 2);

        verticalSpacer_3 = new QSpacerItem(20, 10, QSizePolicy::Minimum, QSizePolicy::Maximum);

        gridLayout->addItem(verticalSpacer_3, 11, 0, 1, 1);

        output = new QPlainTextEdit(Inspect);
        output->setObjectName(QString::fromUtf8("output"));
        QSizePolicy sizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(output->sizePolicy().hasHeightForWidth());
        output->setSizePolicy(sizePolicy);
        output->setMinimumSize(QSize(0, 220));

        gridLayout->addWidget(output, 14, 0, 1, 2);

        gridLayout_3 = new QGridLayout();
        gridLayout_3->setObjectName(QString::fromUtf8("gridLayout_3"));
        gridLayout_3->setSizeConstraint(QLayout::SetMaximumSize);
        x = new QPlainTextEdit(Inspect);
        x->setObjectName(QString::fromUtf8("x"));
        x->setMaximumSize(QSize(81, 31));
        x->setTabChangesFocus(true);
        x->setLineWrapMode(QPlainTextEdit::NoWrap);

        gridLayout_3->addWidget(x, 0, 1, 1, 1);

        label = new QLabel(Inspect);
        label->setObjectName(QString::fromUtf8("label"));
        label->setMaximumSize(QSize(16777215, 31));

        gridLayout_3->addWidget(label, 0, 0, 1, 1, Qt::AlignRight);

        maxX = new QLabel(Inspect);
        maxX->setObjectName(QString::fromUtf8("maxX"));

        gridLayout_3->addWidget(maxX, 0, 2, 1, 1);


        gridLayout->addLayout(gridLayout_3, 3, 0, 1, 2);

        gridLayout_5 = new QGridLayout();
        gridLayout_5->setObjectName(QString::fromUtf8("gridLayout_5"));
        label_4 = new QLabel(Inspect);
        label_4->setObjectName(QString::fromUtf8("label_4"));

        gridLayout_5->addWidget(label_4, 0, 0, 1, 1, Qt::AlignRight);

        z = new QPlainTextEdit(Inspect);
        z->setObjectName(QString::fromUtf8("z"));
        z->setMaximumSize(QSize(80, 31));
        z->setTabChangesFocus(true);
        z->setLineWrapMode(QPlainTextEdit::NoWrap);

        gridLayout_5->addWidget(z, 0, 1, 1, 1);

        maxZ = new QLabel(Inspect);
        maxZ->setObjectName(QString::fromUtf8("maxZ"));

        gridLayout_5->addWidget(maxZ, 0, 2, 1, 1);


        gridLayout->addLayout(gridLayout_5, 7, 0, 1, 2);


        gridLayout_2->addLayout(gridLayout, 0, 0, 1, 1);


        retranslateUi(Inspect);

        QMetaObject::connectSlotsByName(Inspect);
    } // setupUi

    void retranslateUi(QDialog *Inspect)
    {
        Inspect->setWindowTitle(QApplication::translate("Inspect", "Render", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("Inspect", "y:", 0, QApplication::UnicodeUTF8));
        maxY->setText(QApplication::translate("Inspect", "of", 0, QApplication::UnicodeUTF8));
        Return->setText(QApplication::translate("Inspect", "Return", 0, QApplication::UnicodeUTF8));
        InspectCell->setText(QApplication::translate("Inspect", "Inspect Cell", 0, QApplication::UnicodeUTF8));
        status->setText(QApplication::translate("Inspect", "Inspect Cell", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("Inspect", "x:", 0, QApplication::UnicodeUTF8));
        maxX->setText(QApplication::translate("Inspect", "of", 0, QApplication::UnicodeUTF8));
        label_4->setText(QApplication::translate("Inspect", "z:", 0, QApplication::UnicodeUTF8));
        maxZ->setText(QApplication::translate("Inspect", "of", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class Inspect: public Ui_Inspect {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_INSPECT_H
