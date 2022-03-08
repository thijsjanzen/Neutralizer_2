#ifndef MAINWINDOW_HPP
#define MAINWINDOW_HPP

#include <QMainWindow>
#include "qcustomplot.h"
#include "simulation.h"

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
  Q_OBJECT

public:
  MainWindow(QWidget *parent = nullptr);
  ~MainWindow();

  size_t row_size;
  double spec_rate;
  double migr_rate;
  double disp_range;
  int Jm;
  double theta;

  size_t update_speed;

  bool is_running;
  bool is_paused;

  void update_display();

private slots:
  void on_update_params_clicked();

  void on_button_start_clicked();

  void on_button_stop_clicked();

  void on_speed_slider_actionTriggered(int action);

private:

  QVector<double> x_t;
  QVector<double> y_t; // num species

  Ui::MainWindow *ui;
  QImage image_;

  std::unique_ptr<simulation> sim;

  QCPBars *meta_comm_bars;
  QCPBars *local_comm_bars;

  int max_local_comm_bars;
  int max_rank_abund_rank;

  void replot_graphs();
  void set_resolution(int width, int height);
  void update_plots(double t,
                    const std::vector<double>& sp_area_x,
                    const std::vector<double>& sp_area_y);
  void update_params();
};
#endif // MAINWINDOW_HPP
