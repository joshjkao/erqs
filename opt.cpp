#include "common.hpp"
#include "hamiltonian.hpp"
#include "operations.hpp"
#include "optimization.hpp"
#include "quantumstate.hpp"
#include "validation.hpp"
#include <argparse/argparse.hpp>
#include <array>
#include <chrono>
#include <fstream>
#include <iostream>
#include <memory>
#include <print>
#include <random>

long time_in_ms(auto func) {
  auto start = std::chrono::high_resolution_clock::now();
  func();
  auto end = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  return duration.count();
}

void random_update(std::shared_ptr<SumState> root, std::mt19937 &gen) {
  std::uniform_int_distribution<int> dist(0, 2);
  int rand = dist(gen);
  if (rand == 0) {
    std::println("attempting a random term addition");
    AddRandomTerm(root, 3, gen);
  } else if (rand == 1) {
    std::println("attempting a random unification");
    UnifyRandom(root, gen);
  } else {
    std::println("attempting to remove a random term");
    RemoveRandomTerm(root, gen);
  }
}

int main(int argc, char **argv) {
  auto parser = argparse::ArgumentParser{"optimizer", "0.1"};
  parser.add_argument("-H", "--hamiltonian");
  parser.add_argument("-U", "--updates").scan<'d', size_t>();
  parser.add_argument("-L", "--lambda").scan<'g', double>();
  parser.add_argument("-S", "--start-temperature").scan<'g', double>();
  parser.add_argument("-E", "--end-temperature").scan<'g', double>();
  parser.add_argument("-O", "--output").default_value("out.csv");

  parser.parse_args(argc, argv);

  std::string h_file = parser.get<std::string>("--hamiltonian");
  size_t num_updates = parser.get<size_t>("--updates");
  double lambda = parser.get<double>("--lambda");
  std::string out_file = parser.get<std::string>("--output");
  const double T_start = parser.get<double>("--start-temperature");
  const double T_end = parser.get<double>("--end-temperature");

  std::cout << "using hamiltonian " << h_file << "\n";
  std::cout << "performing " << num_updates << "updates\n";
  std::cout << "starting temperature " << T_start << "\n";
  std::cout << "ending temperature " << T_end << "\n";
  std::cout << "lambda is " << lambda << "\n";

  std::random_device rd;
  auto seed = rd();
  std::mt19937 gen{seed};
  auto root = ZeroOneTensor(ONES);
  // auto root = RandomSumState(3, 2, 2, ONES, gen);

  Simplify(root);
  Validate(root, CHECK_ALL);

  Print(root);
  std::cout << std::endl;

  PauliHamiltonian H = ReadFromFile(h_file);

  Print(H);

  auto cost_function = [](std::shared_ptr<SumState> state, double energy,
                          double lambda_) -> double {
    return lambda_ * static_cast<double>(CountNodes(state)) + energy;
  };

  double minf_local = OptimizeCoefficients(H, root);

  double current_cost = cost_function(root, minf_local, lambda);
  std::cout << "localized min " << minf_local << "\n";

  std::vector<double> energies(num_updates, 0);
  std::vector<size_t> n_nodes(num_updates, 0);
  std::vector<double> free_energy(num_updates, 0);

  std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);

  // Track the best state seen so far, since uphill MC moves can temporarily
  // worsen current_cost
  auto best_root = Clone(root);
  double best_cost = current_cost;

  const double max_step =
      static_cast<double>(num_updates > 1uz ? num_updates - 1uz : 1uz);

  for (size_t i = 0; i < num_updates; ++i) {
    // Exponential cooling schedule (set T_start == T_end for fixed-temperature
    // MC)
    double progress = static_cast<double>(i) / max_step;
    double T_mc = T_start * std::pow(T_end / T_start, progress);

    auto root_cpy = Clone(root);
    random_update(root_cpy, gen);
    Simplify(root_cpy);
    Validate(root_cpy, CHECK_ALL);
    Print(root_cpy);

    double new_minf_local = OptimizeCoefficients(H, root_cpy);
    std::cout << "localized min " << new_minf_local << "\n";

    Normalize(root_cpy);

    // cost_function evaluates F = E + lambda * Memory(root_cpy)
    double new_cost = cost_function(root_cpy, new_minf_local, lambda);
    double delta_cost = new_cost - current_cost;

    bool accept = false;
    if (delta_cost <= 0.0) {
      accept = true;
    } else if (T_mc > 0.0) {
      double acceptance_prob = std::exp(-delta_cost / T_mc);
      accept = (uniform_dist(gen) < acceptance_prob);
    }

    if (accept) {
      root = std::move(root_cpy);
      minf_local = new_minf_local;
      current_cost = new_cost;
      std::cout << "change accepted (dF = " << delta_cost << ", T = " << T_mc
                << ")\n";

      if (current_cost < best_cost) {
        best_root = Clone(root);
        best_cost = current_cost;
      }
    } else {
      std::cout << "change rejected (dF = " << delta_cost << ", T = " << T_mc
                << ")\n";
    }

    energies[i] = minf_local;
    n_nodes[i] = CountNodes(root);
    free_energy[i] = current_cost;
  }

  Print(root);
  std::cout << "final minimum energy: " << minf_local << "\n";

  std::ofstream ofile(out_file);
  ofile << "energy,memory,free_energy\n";
  for (size_t i = 0; i < num_updates; ++i) {
    ofile << energies[i] << "," << n_nodes[i] << "," << free_energy[i] << "\n";
  }
  ofile.close();
}
