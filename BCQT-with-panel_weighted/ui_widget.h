/********************************************************************************
** Form generated from reading UI file 'widget.ui'
**
** Created by: Qt User Interface Compiler version 5.6.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_WIDGET_H
#define UI_WIDGET_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QListView>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QTableView>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_Widget
{
public:
    QTableView *tableView;
    QListView *listView;
    QRadioButton *lineSmooth;
    QDoubleSpinBox *nodeTran;
    QLabel *label;
    QLabel *label_2;
    QDoubleSpinBox *nodeSize;
    QLabel *label_3;
    QDoubleSpinBox *edgeTran;
    QLabel *label_4;
    QDoubleSpinBox *edgeSize;

    void setupUi(QWidget *Widget)
    {
        if (Widget->objectName().isEmpty())
            Widget->setObjectName(QStringLiteral("Widget"));
        Widget->resize(800, 800);
        tableView = new QTableView(Widget);
        tableView->setObjectName(QStringLiteral("tableView"));
        tableView->setGeometry(QRect(800, 0, 200, 800));
        listView = new QListView(Widget);
        listView->setObjectName(QStringLiteral("listView"));
        listView->setGeometry(QRect(689, 0, 111, 251));
        lineSmooth = new QRadioButton(Widget);
        lineSmooth->setObjectName(QStringLiteral("lineSmooth"));
        lineSmooth->setGeometry(QRect(700, 0, 100, 20));
        nodeTran = new QDoubleSpinBox(Widget);
        nodeTran->setObjectName(QStringLiteral("nodeTran"));
        nodeTran->setGeometry(QRect(700, 50, 67, 24));
        nodeTran->setMaximum(1);
        nodeTran->setSingleStep(0.1);
        nodeTran->setValue(1);
        label = new QLabel(Widget);
        label->setObjectName(QStringLiteral("label"));
        label->setGeometry(QRect(700, 30, 59, 16));
        label_2 = new QLabel(Widget);
        label_2->setObjectName(QStringLiteral("label_2"));
        label_2->setGeometry(QRect(700, 80, 59, 16));
        nodeSize = new QDoubleSpinBox(Widget);
        nodeSize->setObjectName(QStringLiteral("nodeSize"));
        nodeSize->setGeometry(QRect(700, 100, 67, 24));
        nodeSize->setSingleStep(0.1);
        nodeSize->setValue(1);
        label_3 = new QLabel(Widget);
        label_3->setObjectName(QStringLiteral("label_3"));
        label_3->setGeometry(QRect(700, 130, 59, 16));
        edgeTran = new QDoubleSpinBox(Widget);
        edgeTran->setObjectName(QStringLiteral("edgeTran"));
        edgeTran->setGeometry(QRect(700, 150, 67, 24));
        edgeTran->setSingleStep(0.1);
        edgeTran->setValue(1);
        label_4 = new QLabel(Widget);
        label_4->setObjectName(QStringLiteral("label_4"));
        label_4->setGeometry(QRect(700, 180, 59, 16));
        edgeSize = new QDoubleSpinBox(Widget);
        edgeSize->setObjectName(QStringLiteral("edgeSize"));
        edgeSize->setGeometry(QRect(700, 200, 67, 24));
        edgeSize->setSingleStep(0.1);
        edgeSize->setValue(1);

        retranslateUi(Widget);

        QMetaObject::connectSlotsByName(Widget);
    } // setupUi

    void retranslateUi(QWidget *Widget)
    {
        Widget->setWindowTitle(QApplication::translate("Widget", "Widget", 0));
        lineSmooth->setText(QApplication::translate("Widget", "smooth", 0));
        label->setText(QApplication::translate("Widget", "node transparacy", 0));
        label_2->setText(QApplication::translate("Widget", "node size", 0));
        label_3->setText(QApplication::translate("Widget", "edge tran", 0));
        label_4->setText(QApplication::translate("Widget", "edge size", 0));
    } // retranslateUi

};

namespace Ui {
    class Widget: public Ui_Widget {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_WIDGET_H
