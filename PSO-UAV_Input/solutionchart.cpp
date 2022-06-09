#include "solutionchart.h"
#include "ui_solutionchart.h"

SolutionChart::SolutionChart(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::SolutionChart)
{
    ui->setupUi(this);
    readThreatSourceBmp();
    readSolution();

    chart->legend()->hide();
    for(int i = 0; i < numberofthreatsource; i++){
        chart->addSeries(&threat_series[i]);
    }
    for(int i = 0; i < G_n; i++){
        chart->addSeries(&solution_series[i]);
    }

    chart->createDefaultAxes();
    chart->setTitle("PSO-UAV Solution");

    chartView->setRenderHint(QPainter::Antialiasing);
    chartView->setParent(ui->horizontalFrame);
}

SolutionChart::~SolutionChart()
{
    delete ui;
}

void SolutionChart::on_SolutionChart_destroyed()
{
    delete ui;
}

void SolutionChart::readThreatSourceBmp()
{
    Threat_source_bmp.open("Threat_source_bmp.txt", std::fstream::in);
    int switcher0 = 0;
    int switcher1 = 0;
    int switcher2 = 0;
    double Threat_source_bmp_temp[2][100];
    while (std::getline(Threat_source_bmp, filecontent)) {
        switch (switcher0) {
        case 0:
            numberofthreatsource = std::stoi(filecontent);
            switcher0++;
            break;
        case 1:
            if (switcher2 >= numberofthreatsource)
                break;
            new (&ss) std::stringstream(filecontent);
            for (int i = 0; i < 100; i++){
                ss >> Threat_source_bmp_temp[switcher1][i];
            }
            if (switcher1 == 0)
                switcher1++;
            else{
                switcher1 = 0;
                for (int i = 0; i < 100; i++) {
                    threat_series[switcher2].append(Threat_source_bmp_temp[0][i], Threat_source_bmp_temp[1][i]);
                }
                switcher2++;
            }
            break;
        }
    }
}

void SolutionChart::readSolution()
{
    Solution.open("Solution.txt", std::fstream::in);
    int switcher0 = 0;
    int switcher1 = 0;
    std::getline(Solution, SolutionChart::filecontent);
    SolutionChart::G_n = std::stoi(SolutionChart::filecontent);
    switcher0++;
    std::getline(Solution, SolutionChart::filecontent);
    SolutionChart::N = std::stoi(SolutionChart::filecontent);
    SolutionChart::solution_series = new QtCharts::QLineSeries[G_n];
    switcher0++;
    double xIntervals[N];
    double solution[N];
    while(std::getline(Solution, filecontent)){
        switch(switcher0) {
        case 2:
            new (&ss) std::stringstream(filecontent);
            for (int i = 0; i < N; i++) {
                ss >> xIntervals[i];
            }
            switcher0++;
            break;
        case 3:
            new (&ss) std::stringstream(filecontent);
            for (int i = 0; i < N; i++) {
                ss >> solution[i];
                solution_series[switcher1].append(xIntervals[i], solution[i]);
            }
            switcher1++;
            if(switcher1 == G_n){
                switcher0++;
            }
            break;
        default:
            break;
        }
    }
}
