#pragma once
#include <bitset>
#include <complex>
#include <memory>
#include <variant>

constexpr size_t NQUBITS = 4;
using complex = std::complex<double>;
using QSpace = std::bitset<NQUBITS>;
using BitString = std::bitset<NQUBITS>;

constexpr QSpace ONES{~QSpace{0}};
constexpr QSpace ZEROS{0};

struct PureState;
struct ProductState;
struct SumState;

using ptr_variant =
    std::variant<std::shared_ptr<PureState>, std::shared_ptr<ProductState>>;
