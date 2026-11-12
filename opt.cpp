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
  parser.add_argument("-O", "--output").default_value("out.csv");

  parser.parse_args(argc, argv);

  std::string h_file = parser.get<std::string>("--hamiltonian");
  size_t num_updates = parser.get<size_t>("--updates");
  double lambda = parser.get<double>("--lambda");
  std::string out_file = parser.get<std::string>("--output");

  std::cout << "using hamiltonian " << h_file << "\n";
  std::cout << "performing " << num_updates << "updates\n";

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

  auto cost_function = [=](std::shared_ptr<SumState> state,
                          double energy) -> double {
    return lambda * static_cast<double>(CountNodes(state)) + energy;
  };

  double minf_local = OptimizeCoefficients(H, root);

  double current_cost = cost_function(root, minf_local);
  std::cout << "localized min " << minf_local << "\n";

  std::vector<double> energies(num_updates, 0);
  std::vector<size_t> n_nodes(num_updates, 0);
  std::vector<double> free_energy(num_updates, 0);

  for (size_t i = 0; i < num_updates; ++i) {
    auto root_cpy = Clone(root);
    random_update(root_cpy, gen);
    Simplify(root_cpy);
    Validate(root_cpy, CHECK_ALL);
    Print(root_cpy);

    // InvokeBottomUp(visitor, root_cpy);
    double new_minf_local = OptimizeCoefficients(H, root_cpy);
    std::cout << "localized min " << new_minf_local << "\n";

    Normalize(root_cpy);

    double new_cost = cost_function(root_cpy, new_minf_local);

    if (new_cost <= current_cost) {
      root = root_cpy;
      minf_local = new_minf_local;
      current_cost = new_cost;
      std::cout << "change accepted\n";
    } else {
      std::cout << "change rejected\n";
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
