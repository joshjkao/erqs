#pragma once
#include <bitset>
#include <complex>
#include <memory>
#include <variant>

constexpr size_t NQUBITS = 16;
using complex = std::complex<double>;
using QSpace = std::bitset<NQUBITS>;
using BitString = std::bitset<NQUBITS>;

struct PureState;
struct ProductState;
struct SumState;

using ptr_variant =
    std::variant<std::shared_ptr<PureState>, std::shared_ptr<ProductState>>;
