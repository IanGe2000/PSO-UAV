#include "mainwindow.h"
#include "./ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_lineEdit_numberofthreatsources_textChanged(const QString &s)
{
    if(ui->scrollAreaWidgetContents->layout()){
        delete ui->scrollAreaWidgetContents->layout();
        delete [] threat;
    }
    int number = s.toInt();
    QVBoxLayout *layout1 = new QVBoxLayout();
    threat = new threatsource[10];
    if (number <= 10){
        for (int i = 0; i < number; i++) {
            layout1->addWidget(&threat[i]);
        }
    }
    ui->scrollAreaWidgetContents->setLayout(layout1);
}

void MainWindow::on_pushButton_OK_clicked()
{
    Input_parameters.open("Input_parameters.txt", std::fstream::out | std::fstream::trunc);
    Input_parameters << ui->lineEdit_G_n->text().toInt() << "\n";
    Input_parameters << ui->lineEdit_startingpoint_x->text().toDouble() << " ";
    Input_parameters << ui->lineEdit_startingpoint_y->text().toDouble() << "\n";
    Input_parameters << ui->lineEdit_endpoint_x->text().toDouble() << " ";
    Input_parameters << ui->lineEdit_endpoint_y->text().toDouble() << "\n";
    int number = ui->lineEdit_numberofthreatsources->text().toInt();
    Input_parameters << number << "\n";
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < number; j++) {
            Input_parameters << threat[j].getvalue(i) << " ";
        }
        Input_parameters << "\n";
    }
    Input_parameters << ui->lineEdit_n->text().toDouble() << "\n";
    Input_parameters << ui->lineEdit_N->text().toDouble() << "\n";
    Input_parameters << ui->lineEdit_thetaTmax->text().toDouble() << " ";
    Input_parameters << ui->lineEdit_thetaCmax->text().toDouble() << "\n";
    Input_parameters << ui->lineEdit_maxgeneration->text().toDouble() << "\n";
    Input_parameters.close();
    process.start(QString::fromStdString("./CXX_simulation"));
    process.waitForFinished();
    SolutionChart *solution = new SolutionChart;
    solution->show();
}

void MainWindow::on_pushButton_loaddefault_clicked()
{
    ui->lineEdit_G_n->setText(QString::number(5));
    ui->lineEdit_startingpoint_x->setText(QString::number(0.0));
    ui->lineEdit_startingpoint_y->setText(QString::number(0.0));
    ui->lineEdit_endpoint_x->setText(QString::number(30.0));
    ui->lineEdit_endpoint_y->setText(QString::number(30.0));
    ui->lineEdit_numberofthreatsources->setText(QString::number(7));
    double threatsouce_value[3][7] = {
        {6.2280,17.7810,15.6810, 6.5280,22.5810,15.0570,21.0360},
        {8.5230, 4.6080,17.2080,13.6290,21.1080,11.8350,15.8460},
        {2.2826, 1.9663, 2.8540, 2.0762, 1.9393, 2.4483, 2.4404}
    };
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 7; j++){
            threat[j].setvalue(i,threatsouce_value[i][j]);
        }
    }
    ui->lineEdit_n->setText(QString::number(10));
    ui->lineEdit_N->setText(QString::number(10));
    ui->lineEdit_thetaTmax->setText(QString::number(60));
    ui->lineEdit_thetaCmax->setText(QString::number(45));
    ui->lineEdit_maxgeneration->setText(QString::number(30));
}

void MainWindow::on_MainWindow_destroyed()
{
    delete ui;
}
