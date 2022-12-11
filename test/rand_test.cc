#include "copter/random.h"
#include <iostream>
#include <random>
template <typename ContainerT>
void PrintContainer(const ContainerT &container) {
  for (const auto &elem : container) {
    std::cout << elem << " ";
  }
  std::cout << "\n";
}

int main(void) {

  size_t test_cnt = 10;

  for (size_t i = 0; i < test_cnt; ++i) {

    constexpr std::size_t size = 5;

    std::cout << "-------------test no." << i << "----------------\n";

    // float random vector
    std::uniform_real_distribution<double> uniform_f_dist(0, 1);
    auto rand_float_vec =
        copter::random::GenRandomVec<float>(10, uniform_f_dist);
    PrintContainer(rand_float_vec);

    // int random vector
    std::uniform_int_distribution<int> uniform_i_dist(0, 10);
    auto rand_int_vec = copter::random::GenRandomVec<int>(size, uniform_i_dist);
    PrintContainer(rand_int_vec);
  }

  return 0;
}