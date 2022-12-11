#include "copter/brkga.h"
#include <cmath>
#include <fstream>
#include <ios>
#include <vector>

constexpr double PI = 3.1415926525;

int main(void){

    auto fit_func = [&](const std::vector<double> &features) {
    double x = features[0];
    double y = features[1];
    double temp = x * x + y * y;
    double cost = 0.5 + (std::pow(std::sin(temp), 2) - 0.5) / std::pow((1 + 0.001 * temp), 2);
    return -cost;
  };

  copter::brkga::BRKGAParam param;
  
  param.n_feature = 2;
  param.pop_size = 50;
  param.elite_rate = 0.20;
  param.mutant_rate = 0.20;
  param.elitism_prob = 0.60;

  param.max_iter = 500;
  param.lower_bound = std::vector<double>(param.n_feature, -1.0);
  param.upper_bound = std::vector<double>(param.n_feature, 5.0);

  copter::brkga::BRKGA<decltype(fit_func)> brkga(fit_func, param);

  // brkga.Optimize(std::cout);
  std::ofstream log_file("brkga.csv", std::fstream::out);
  brkga.Optimize(log_file);
  log_file.close();
  auto best_feature= brkga.GetBestFeature();
  std::cout << "x: " << best_feature[0] << " y: " << best_feature[1] << "\n";

  // brkga with init
  std::vector<double> f1{0.0, 0.0}, f2{0.0, 1.3};
  std::vector<copter::brkga::Chromo> init_pop{brkga.Encode(f1), brkga.Encode(f2)};
  brkga.InitPop(init_pop);
  std::ofstream log_file2("init_brkga.csv", std::fstream::out);
  brkga.Optimize(log_file2);
  best_feature= brkga.GetBestFeature();
  std::cout << "x: " << best_feature[0] << " y: " << best_feature[1] << "\n";

   
  return 0;
}