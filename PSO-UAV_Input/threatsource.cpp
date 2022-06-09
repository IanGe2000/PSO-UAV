#include "threatsource.h"
#include "ui_threatsource.h"

threatsource::threatsource(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::threatsource)
{
    ui->setupUi(this);
}

threatsource::~threatsource()
{
    delete ui;
}

double threatsource::getvalue(int index)
{
    double value = 0;
    switch (index) {
    case 0:
        value = ui->lineEdit_x->text().toDouble();
        break;
    case 1:
        value = ui->lineEdit_y->text().toDouble();
        break;
    case 2:
        value = ui->lineEdit_r->text().toDouble();
        break;
    }
    return value;
}

void threatsource::setvalue(int index, double value)
{
    switch (index) {
    case 0:
        ui->lineEdit_x->setText(QString::number(value));
        break;
    case 1:
        ui->lineEdit_y->setText(QString::number(value));
        break;
    case 2:
        ui->lineEdit_r->setText(QString::number(value));
        break;
    }
}
