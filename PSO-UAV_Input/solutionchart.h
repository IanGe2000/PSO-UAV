#ifndef SOLUTIONCHART_H
#define SOLUTIONCHART_H

#include <QWidget>
#include <QtCharts/QLineSeries>
#include <QtCharts/QChartView>
#include <fstream>
#include <sstream>
#include <iostream>

namespace Ui {
class SolutionChart;
}

class SolutionChart : public QWidget
{
    std::fstream Threat_source_bmp;
    std::fstream Solution;
    std::string filecontent;
    int numberofthreatsource;
    int G_n;
    int N;
    std::stringstream ss;
    Q_OBJECT
    QtCharts::QLineSeries *threat_series = new QtCharts::QLineSeries[10];
    QtCharts::QLineSeries *solution_series = new QtCharts::QLineSeries[10];
    QtCharts::QChart *chart = new QtCharts::QChart();
    QtCharts::QChartView *chartView = new QtCharts::QChartView(chart);


public:
    explicit SolutionChart(QWidget *parent = nullptr);
    ~SolutionChart();
    void readThreatSourceBmp();
    void readSolution();

private slots:
    void on_SolutionChart_destroyed();

private:
    Ui::SolutionChart *ui;
};

#endif // SOLUTIONCHART_H
