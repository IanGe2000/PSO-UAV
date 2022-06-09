#ifndef THREATSOURCE_H
#define THREATSOURCE_H

#include <QWidget>

namespace Ui {
class threatsource;
}

class threatsource : public QWidget
{
    Q_OBJECT

public:
    explicit threatsource(QWidget *parent = nullptr);
    ~threatsource();
    double getvalue(int index);
    void setvalue(int index, double value);

private:
    Ui::threatsource *ui;
};

#endif // THREATSOURCE_H
