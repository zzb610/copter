#include "copter/random.h"
#include <algorithm>
#include <cstddef>
#include <numeric>
#include <ostream>
#include <utility>
#include <vector>

namespace copter {
namespace brkga {

struct SearchProccess {};

struct BRKGAParam {
  std::size_t max_iter;
  std::size_t early_stop;

  std::size_t pop_size;
  std::size_t n_feature;

  double elite_rate;
  double mutant_rate;
  double elitism_prob;

  std::vector<double> lower_bound;
  std::vector<double> upper_bound;
};

struct Chromo {
  Chromo()=default;
  Chromo(const std::vector<double> &g, double fit) : genes(g), fitness(fit) {}
  std::vector<double> genes;
  double fitness;
};

template <typename FitFuncT> class BRKGA {

public:
  BRKGA(FitFuncT fit_func, BRKGAParam param)
      : fit_func_(std::move(fit_func)), param_(std::move(param)){};

  Chromo Encode(const std::vector<double> &feature) const {
    std::vector<double> gene(param_.n_feature, 1.0);
    for (std::size_t i = 0; i < gene.size(); ++i) {
      gene[i] *= (feature[i] - param_.lower_bound[i]) / param_.upper_bound[i];
    }
    auto fit = fit_func_(feature);
    return Chromo(gene, fit);
  }

  std::vector<double> Decode(const std::vector<double> &genes) const {
    std::vector<double> feature;
    for (std::size_t i = 0; i < genes.size(); ++i) {
      auto range_len = param_.upper_bound[i] - param_.lower_bound[i];
      feature.push_back(genes[i] * range_len + param_.lower_bound[i]);
    }
    return feature;
  }

  double GetFitness(const std::vector<double> &gene) const {
    auto feature = Decode(gene);
    return fit_func_(feature);
  }

  void InitPop(std::vector<Chromo> init_pop) {
    population_ = std::move(init_pop);
    auto init_num = init_pop.size();

    std::uniform_real_distribution<double> dist(0, 1);
    for (std::size_t i = 0; i < param_.pop_size - init_num; ++i) {
      auto rand_gene =
          copter::random::GenRandomVec<double>(param_.n_feature, dist);
      auto fitness = GetFitness(rand_gene);
      population_.push_back(Chromo(rand_gene, fitness));
    }
  }

  std::vector<double> GetBestFeature() const {
    return Decode(best_chromo_.genes);
  }

  void Optimize(std::ostream &log) {
    // init population
    if (population_.empty()) {
      InitPopRandom();
    }

    auto chromo_less = [](const Chromo &lhs, const Chromo &rhs) {
      return lhs.fitness < rhs.fitness;
    };

    best_chromo_ = *std::max_element(population_.cbegin(),
                                          population_.cend(), chromo_less);

    std::vector<std::size_t> elite_idxes, ord_idxes;
    std::cout << "generation_best,best\n";
    for (std::size_t i = 0; i < param_.max_iter; ++i) {

      // next generation
      std::vector<Chromo> new_pop;

      DivPop(elite_idxes, ord_idxes);
      // copy elites
      for (const auto &idx : elite_idxes) {
        new_pop.push_back(population_[idx]);
      }

      std::size_t mutant_num =
          static_cast<std::size_t>(param_.pop_size * param_.mutant_rate);
      auto cross_num = param_.pop_size - elite_idxes.size() - mutant_num;

      // crossover
      for (std::size_t i = 0; i < cross_num; ++i) {
        auto elite_idx = random::RandomChoice(elite_idxes);
        auto ord_idx = random::RandomChoice(ord_idxes);
        auto child = CrossOver(population_[elite_idx], population_[ord_idx]);
        new_pop.push_back(std::move(child));
      }

      // generate mutants
      std::uniform_real_distribution<double> dist(0, 1);
      for (std::size_t i = 0; i < mutant_num; ++i) {
        auto mutant_genes =
            random::GenRandomVec<double>(param_.n_feature, dist);
        auto fit = GetFitness(mutant_genes);
        new_pop.push_back(Chromo(mutant_genes, fit));
      }

      // evolve to next generation
      population_ = std::move(new_pop);

      auto gen_best_chromo_it = std::max_element(
          population_.cbegin(), population_.cend(), chromo_less);
      if (gen_best_chromo_it->fitness > best_chromo_.fitness) {
        best_chromo_ = *gen_best_chromo_it;
      }
      log << gen_best_chromo_it->fitness << "," << best_chromo_.fitness << "\n";
    }

  }

  Chromo CrossOver(const Chromo &elite_p, const Chromo &ord_p) {

    auto n_genes = elite_p.genes.size();
    std::uniform_real_distribution<double> dist(0, 1);
    auto mask = random::GenRandomVec<double>(n_genes, dist);

    std::vector<double> child_genes(n_genes);
    for (std::size_t i = 0; i < n_genes; ++i) {
      if (mask[i] < param_.elitism_prob) {
        child_genes[i] = elite_p.genes[i];
      } else {
        child_genes[i] = ord_p.genes[i];
      }
    }
    auto child_fit = GetFitness(child_genes);
    return Chromo(child_genes, child_fit);
  }

  void InitPopRandom() {
    std::uniform_real_distribution<double> dist(0, 1);
    for (std::size_t i = 0; i < param_.pop_size; ++i) {
      auto rand_gene =
          copter::random::GenRandomVec<double>(param_.n_feature, dist);
      auto fitness = GetFitness(rand_gene);
      population_.push_back(Chromo(rand_gene, fitness));
    }
  }

  void DivPop(std::vector<std::size_t> &elite_idx,
              std::vector<std::size_t> &ord_idx) {

    // argsort
    std::vector<std::size_t> idx(param_.pop_size);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(),
              [this](std::size_t lhs, const std::size_t rhs) {
                return population_[lhs].fitness > population_[rhs].fitness;
              });

    // elites and ordinaries(non-elites)
    auto elite_num =
        static_cast<std::size_t>(param_.pop_size * param_.elite_rate);
    elite_idx.clear();
    ord_idx.clear();
    elite_idx.assign(idx.begin(), idx.begin() + elite_num);
    ord_idx.assign(idx.begin() + elite_num, idx.end());
  }

  FitFuncT fit_func_;
  BRKGAParam param_;
  std::vector<Chromo> population_;
  Chromo best_chromo_;
};

} // namespace brkga
} // namespace copter