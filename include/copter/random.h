#pragma once

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <random>
#include <vector>

namespace copter {
namespace random {

template <typename ElemT, typename DistT>
ElemT GenRandomNum(DistT& dist) {
  std::random_device seed;
  std::ranlux48 engine(seed());
  return dist(engine);
}

template <typename ElemT, typename DistT>
std::vector<ElemT> GenRandomVec(std::size_t size, DistT dist) {
  std::random_device seed;
  std::ranlux48 engine(seed());
  auto generator = [&dist, &engine] { return dist(engine); };
  std::vector<ElemT> rand_vec(size);
  std::generate(rand_vec.begin(), rand_vec.end(), generator);
  return rand_vec;
}

template <typename ElemT> ElemT RandomChoice(const std::vector<ElemT> &arr) {
  auto n_elem = arr.size();
  std::uniform_int_distribution<int> dist(0, n_elem);
  auto rand_idx = GenRandomNum<int>(n_elem);
  return arr.at(rand_idx);
}

template <typename ElemT> ElemT RandomChoice(std::vector<ElemT> &arr) {
  auto n_elem = arr.size();
  std::uniform_int_distribution<int> uniform_dist(0, n_elem-1);
  auto rand_idx = GenRandomNum<int>(uniform_dist);
  return arr.at(rand_idx);
}

} // namespace random
} // namespace copter
