#pragma once
#include "quantumstate.h"
#include <memory>

void SetCoefficients(const std::vector<complex> &coeffs,
                     std::shared_ptr<SumState> &state);
void SetCoefficients(const std::vector<double> &coeffs,
                     std::shared_ptr<SumState> &state);
size_t CountCoefficients(const std::shared_ptr<SumState> &state);
size_t CountNodes(const std::shared_ptr<SumState> &root);

// ---- COEFFICIENT OPTIMIZATION ---- //
double OptimizeCoefficients(const PauliHamiltonian &H,
                            std::shared_ptr<SumState> &state);

// ---- STRUCTURE CHANGING ADJUSTMENTS ---- //
void Prune(std::shared_ptr<SumState> &root, double threshold);
void AddRandomTerm(std::shared_ptr<SumState> &root, std::mt19937 &gen);
void UnifyRandom(std::shared_ptr<SumState> &root, std::mt19937 &gen);
void Normalize_slow(std::shared_ptr<SumState> &root);
void Normalize(std::shared_ptr<SumState> &root);
void RemoveSingles(std::shared_ptr<SumState> &root);
void SquashPures(std::shared_ptr<SumState> &root);
// void Factorize(std::shared_ptr<SumState> &root);

void Simplify(ptr_variant &var);
void Simplify(std::shared_ptr<SumState> &ptr);
