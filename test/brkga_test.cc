#include "copter/brkga.h"
#include <cmath>
#include <fstream>
#include <vector>

constexpr double PI = 3.1415926525;

int main(void){

    auto fit_func = [&](const std::vector<double> &features) {
    double x = features[0];
    double y = features[1];
    double temp = x * x + y * y;
    double cost = 0.5 + (std::pow(std::sin(temp), 2) - 0.5) / std::pow((1 + 0.001 * temp), 2);
    return  1 / cost;
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

  brkga.Optimize(std::cout);

  auto best_feature= brkga.GetBestFeature();
  std::cout << "x: " << best_feature[0] << " y: " << best_feature[1] << "\n";

   
  return 0;
}