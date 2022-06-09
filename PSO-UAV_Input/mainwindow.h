#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <iostream>
#include <fstream>
#include <chrono>
#include <thread>
#include <QProcess>
#include "threatsource.h"
#include "solutionchart.h"

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    std::fstream Input_parameters;
    QProcess process;
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
    threatsource *threat = new threatsource[10];

private slots:
    void on_lineEdit_numberofthreatsources_textChanged(const QString &arg1);

    void on_pushButton_OK_clicked();

    void on_pushButton_loaddefault_clicked();

    void on_MainWindow_destroyed();

private:
    Ui::MainWindow *ui;
};
#endif // MAINWINDOW_H
