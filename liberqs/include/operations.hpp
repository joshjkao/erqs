#pragma once
#include "common.hpp"
#include "hamiltonian.hpp"
#include "quantumstate.hpp"
#include "skewoperator.hpp"

// ---- INNER PRODUCTS ---- //
SkewOperator Inner(const std::shared_ptr<PureState> &p1,
                   const std::shared_ptr<PureState> &p2);
SkewOperator Inner(const std::shared_ptr<ProductState> &p1,
                   const std::shared_ptr<ProductState> &p2);
SkewOperator Inner(const std::shared_ptr<ProductState> &p1,
                   const std::shared_ptr<PureState> &p2);
SkewOperator Inner(const std::shared_ptr<PureState> &p1,
                   const std::shared_ptr<ProductState> &p2);
SkewOperator Inner(const ptr_variant &p1, const ptr_variant &p2);
SkewOperator Inner(const std::shared_ptr<SumState> &p1,
                   const std::shared_ptr<SumState> &p2);
SkewOperator CompressConstants(const SkewOperator &op);
SkewOperator Simplify(const SkewOperator &op);

void Add(SkewOperator &o1, const SkewOperator &o2);
void Multiply(SkewOperator &o, const complex &c);
SkewOperator Multiply(const SkewOperator &o1, const SkewOperator &o2);

ptr_variant Tensor(const ptr_variant &p1, const ptr_variant &p2);
std::shared_ptr<SumState> Tensor(const std::shared_ptr<SumState> &p1,
                                 const std::shared_ptr<SumState> &p2);

// ---- OPERATIONS ---- //
std::shared_ptr<SumState> Operate(const PauliOperator &pauli,
                                  const std::shared_ptr<SumState> &state);
std::shared_ptr<SumState> Operate(const PauliHamiltonian &h,
                                  const std::shared_ptr<SumState> &state);
double Norm(const std::shared_ptr<SumState> &state);
double ExpectedValue(const PauliOperator &pauli,
                     const std::shared_ptr<SumState> &state);
double ExpectedValue(const PauliHamiltonian &h,
                     const std::shared_ptr<SumState> &state);
double ExpectedValue_slow(const PauliHamiltonian &h,
                          const std::shared_ptr<SumState> &state);
