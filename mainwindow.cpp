#include "mainwindow.hpp"
#include "ui_mainwindow.h"

#include "simulation.h"
#include <sstream>

#include <memory>

MainWindow::MainWindow(QWidget *parent)
  : QMainWindow(parent)
  , ui(new Ui::MainWindow)
{
  ui->setupUi(this);

   ui->plot_species->addGraph();
   ui->plot_species->graph(0)->setPen(QPen(Qt::black));

  // ui->plot_species->graph(0)->setName("Num Species");

 //  QCPPlotTitle *fst_title = new QCPPlotTitle(ui->plot_species, "Number of species");
  // fst_title->setFont(QFont("sans", 12, QFont::Bold));
  // ui->plot_species->plotLayout()->insertRow(0);
  // ui->plot_species->plotLayout()->addElement(0, 0, fst_title);
   ui->plot_species->xAxis->setLabel("Time");
   ui->plot_species->yAxis->setLabel("Number of species");

   ui->plot_rankabund->addGraph();
   ui->plot_rankabund->graph(0)->setPen(QPen(Qt::black));

 //  ui->plot_rankabund->graph(0)->setName("Num Species");

  // QCPPlotTitle *fst_title2 = new QCPPlotTitle(ui->plot_rankabund, "Rank Abundance");
 //  fst_title2->setFont(QFont("sans", 12, QFont::Bold));
  // ui->plot_rankabund->plotLayout()->insertRow(0);
  // ui->plot_rankabund->plotLayout()->addElement(0, 0, fst_title2);
   ui->plot_rankabund->xAxis->setLabel("Rank");
   ui->plot_rankabund->yAxis->setLabel("Relative abundance");
   max_rank_abund_rank = 0;



   ui->plot_meta_comm->addGraph();

   meta_comm_bars = new QCPBars(ui->plot_meta_comm->xAxis,
                                ui->plot_meta_comm->yAxis);

   meta_comm_bars->setName("Metacommunity");

   meta_comm_bars->setPen(QPen(Qt::black));
   meta_comm_bars->setBrush(QBrush(QColor(0,0,255, static_cast<int>(0.9 * 255))));

 //  ui->plot_meta_comm->graph(0)->setName("Num Species");
 //  meta_comm_bars->setAntialiased(false);
 //  meta_comm_bars->setAntialiasedFill(false);

 //  QCPPlotTitle *fst_title3 = new QCPPlotTitle(ui->plot_meta_comm, "Metacommunity");
 //  fst_title3->setFont(QFont("sans", 12, QFont::Bold));
 //  ui->plot_meta_comm->plotLayout()->insertRow(0);
 //  ui->plot_meta_comm->plotLayout()->addElement(0, 0, fst_title3);
   ui->plot_meta_comm->xAxis->setLabel("Number of Individuals per Species");
   ui->plot_meta_comm->yAxis->setLabel("Number of Species");

   ui->plot_local_comm->addGraph();

   local_comm_bars = new QCPBars(ui->plot_local_comm->xAxis,
                                 ui->plot_local_comm->yAxis);

   local_comm_bars->setName("Local community");

   local_comm_bars->setPen(QPen(Qt::black));
   local_comm_bars->setBrush(QBrush(QColor(0,0,255, static_cast<int>(0.9 * 255))));

 //  ui->plot_local_comm->graph(0)->setName("Local Community");
  // local_comm_bars->setAntialiased(false);
  // local_comm_bars->setAntialiasedFill(false);

  // QCPPlotTitle *fst_title4 = new QCPPlotTitle(ui->plot_local_comm, "Local community");
  // fst_title4->setFont(QFont("sans", 12, QFont::Bold));
 //  ui->plot_local_comm->plotLayout()->insertRow(0);
  // ui->plot_local_comm->plotLayout()->addElement(0, 0, fst_title4);
   ui->plot_local_comm->xAxis->setLabel("Number of Individuals per Species");
   ui->plot_local_comm->yAxis->setLabel("Number of Species");
}

MainWindow::~MainWindow()
{
  delete ui;
}

std::string get_string(std::string s, float v) {
    std::string output = s + " " + std::to_string(v) + "\n";
    return output;
}

void update_preston_plot(QCustomPlot* UI,
                         QCPBars* barplot,
                         const std::vector<int>& octaves,
                         int& running_max_y) {
  QVector<double> m_y(octaves.size());
  QVector<double> m_x;
  double max_y = 0.0;
  for (int i = 0; i < m_y.size(); ++i) {
      m_x.push_back(static_cast<double>(i));
      m_y[i] = static_cast<double>(octaves[i]);
      if (m_y[i] > max_y) max_y = m_y[i];
  }

  if (max_y > running_max_y)  running_max_y = max_y;

  barplot->setData(m_x, m_y);
  UI->xAxis->setRange(-0.5, m_y.size() * 1.01);
  UI->yAxis->setRange(0.0, running_max_y * 1.1);

  UI->replot();
  UI->update();
}

void MainWindow::on_update_params_clicked()
{
   update_speed = ui->speed_slider->value();
   row_size = ui->box_size->value();
   spec_rate = ui->box_spec_rate->value();
   migr_rate = ui->box_migration_rate->value();
   disp_range = ui->box_dispersal->value();
   Jm = 1e6;
   theta = 42;

   bool init_mono_dom = ui->checkBox->checkState();

   set_resolution(row_size, row_size);

    sim = std::make_unique<simulation>(row_size,
                                       spec_rate,
                                       migr_rate,
                                       Jm,
                                       disp_range,
                                       theta,
                                       init_mono_dom);
    auto dummy_max_y = 0;
    update_preston_plot(ui->plot_meta_comm,
                        meta_comm_bars,
                        sim->get_meta_octaves(),
                        dummy_max_y);

    update_plots(0);
    ui->button_start->setText("Start");
    auto s = std::to_string(sim->num_species());
    ui->label_sp->setText(QString::fromStdString(s));
    x_t.clear();
    y_t.clear();

    update_display();

    return;
}

void MainWindow::set_resolution(int width, int height) {
    image_ = QImage(width, height, QImage::Format_RGB32);
}

QRgb convert_color(const std::array<size_t, 3>& input) {
  QColor col = {static_cast<int>(input[0]),
                static_cast<int>(input[1]),
                static_cast<int>(input[2])};
  return col.rgb();
}

void MainWindow::update_display() {
  size_t line_size = row_size;
  size_t num_lines = row_size;

  for(size_t i = 0; i < num_lines; ++i) {
      QRgb* row = (QRgb*) image_.scanLine(i);

      size_t start = i * line_size;
      size_t end = start + line_size;

      for(size_t index = start; index < end; ++index) {
          size_t local_index = index - start;
          auto local_color = sim->get_color(index);
          auto converted_color = convert_color(local_color);


          row[local_index] = converted_color; // convert_color(sim->get_color(index));
      }
  }

  int w = ui->q_label->width();
  int h = ui->q_label->height();

  ui->q_label->setPixmap((QPixmap::fromImage(image_)).scaled(w, h, Qt::KeepAspectRatio));
  ui->q_label->update();
  QApplication::processEvents();
}



void MainWindow::update_plots(double t) {
  x_t.append(t);
  y_t.append(sim->num_species());
  ui->plot_species->graph(0)->clearData();
  ui->plot_species->graph(0)->setData(x_t, y_t);

 // ui->plot_species->rescaleAxes();
  auto max_y = *std::max_element(y_t.begin(), y_t.end());
  auto max_x = x_t.back();
  ui->plot_species->yAxis->setRange(0, max_y);
  ui->plot_species->xAxis->setRange(0, max_x);

  ui->plot_species->xAxis->setAutoTickStep(false);
  double tick_step = max_x / 8;
  ui->plot_species->xAxis->setTickStep(tick_step);
  ui->plot_species->xAxis->setTickLabelRotation(45);


  ui->plot_species->replot();


  ui->plot_rankabund->graph(0)->clearData();

  QVector<double> r_x(sim->rank_abund_curve.size());
  QVector<double> r_y(sim->rank_abund_curve.size());
  for (size_t i = 0; i < sim->rank_abund_curve.size(); ++i) {
      r_x[i] = i;
      r_y[i] = sim->rank_abund_curve[i];
  }
  if (sim->rank_abund_curve.size() > max_rank_abund_rank)
    max_rank_abund_rank = sim->rank_abund_curve.size();
  ui->plot_rankabund->graph(0)->setData(r_x, r_y);

 // ui->plot_rankabund->rescaleAxes();
  ui->plot_rankabund->yAxis->setRange(0, 101);
  ui->plot_rankabund->xAxis->setRange(0, max_rank_abund_rank);
  ui->plot_rankabund->replot();

  update_preston_plot(ui->plot_local_comm,
                      local_comm_bars,
                      sim->get_local_octaves(),
                      max_local_comm_bars);

  auto s = std::to_string(sim->num_species());
  ui->label_sp->setText(QString::fromStdString(s));
}

void MainWindow::on_button_start_clicked()
{
  is_running = true;
  while(true) {
    sim->update();
    if (sim->t % update_speed == 0) {
        sim->update_stats();
        update_display();
        update_plots(sim->t);
        std::stringstream s;
        if(!is_running) break;
    }
  }
}

void MainWindow::on_button_stop_clicked()
{
      is_running = false;
      is_paused = true;
      ui->button_start->setText("Resume");
}

void MainWindow::on_speed_slider_actionTriggered(int action) {
   update_speed = ui->speed_slider->value();
}
