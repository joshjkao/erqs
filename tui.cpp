// Copyright 2020 Arthur Sonzogni. All rights reserved.
// Use of this source code is governed by the MIT license that can be found in
// the LICENSE file.
#include <random>
#include <sstream>
#include <string> // for operator+, to_string

#include "optimization.h"
#include "quantumstate.h"

#include "ftxui/component/component.hpp" // for Button, Horizontal, Renderer
#include "ftxui/component/component_base.hpp"     // for ComponentBase
#include "ftxui/component/screen_interactive.hpp" // for ScreenInteractive
#include "ftxui/dom/elements.hpp" // for separator, gauge, text, Element, operator|, vbox, border

using namespace ftxui;

// Z local terms
PauliOperator op1{
    .x = BitString{"0000"}, .y = BitString{"0000"}, .z = BitString{"1000"}};
PauliOperator op2{
    .x = BitString{"0000"}, .y = BitString{"0000"}, .z = BitString{"0100"}};
PauliOperator op3{
    .x = BitString{"0000"}, .y = BitString{"0000"}, .z = BitString{"0010"}};
PauliOperator op4{
    .x = BitString{"0000"}, .y = BitString{"0000"}, .z = BitString{"0001"}};

// XX interaction terms
PauliOperator int1{
    .x = BitString{"1100"}, .y = BitString{"0000"}, .z = BitString{"0000"}};
PauliOperator int2{
    .x = BitString{"0110"}, .y = BitString{"0000"}, .z = BitString{"0000"}};
PauliOperator int3{
    .x = BitString{"0011"}, .y = BitString{"0000"}, .z = BitString{"0000"}};
PauliOperator int4{
    .x = BitString{"1001"}, .y = BitString{"0000"}, .z = BitString{"0000"}};

PauliHamiltonian H{{0, 0, 0, 0, -0.5, -0.5, -0.5, -0.5},
                   {op1, op2, op3, op4, int1, int2, int3, int4}};

ButtonOption Style() {
  auto option = ButtonOption::Simple();
  option.transform = [](const EntryState &s) {
    auto element = text(s.label);
    if (s.focused) {
      element |= bold;
    }
    return element | center | borderEmpty | flex;
  };
  return option;
}

int main() {
  std::mt19937 gen;
  auto root = RandomSumState(3, 2, 2, ~QSpace{0}, gen);
  auto root_prev = Clone(root);

  auto btn_make_rand = Button(
      "make random",
      [&] {
        root_prev = Clone(root);
        root = RandomSumState(3, 2, 2, ~QSpace{0}, gen);
      },
      Style());
  auto btn_random_term = Button(
      "add random term",
      [&] {
        root_prev = Clone(root);
        AddRandomTerm(root, gen);
      },
      Style());
  auto btn_unify = Button(
      "unify random",
      [&] {
        root_prev = Clone(root);
        UnifyRandom(root, gen);
      },
      Style());
  auto btn_clean = Button(
      "cleanup",
      [&] {
        root_prev = Clone(root);
        Prune(root, 1e-3);
        RemoveSingles(root);
        SquashPures(root);
      },
      Style());
  auto btn_opt = Button(
      "optimize these coeffs",
      [&] {
        root_prev = Clone(root);
        Normalize(root);
        Prune(root, 1e-3);
        OptimizeCoefficients(H, root);
      },
      Style());

  // The tree of components. This defines how to navigate using the keyboard.
  // The selected `row` is shared to get a grid layout.
  int row = 0;
  auto buttons = Container::Vertical({
      Container::Horizontal({btn_make_rand, btn_random_term}, &row) | flex,
      Container::Horizontal({btn_unify, btn_clean}, &row) | flex,
      Container::Horizontal({btn_opt}, &row) | flex,
  });

  // Modify the way to render them on screen:
  auto component = Renderer(buttons, [&] {
    std::stringstream ss;
    PrintToStream(ss, root);
    std::stringstream ss_prev;
    PrintToStream(ss_prev, root_prev);
    std::string s = ss.str();
    std::string s_prev = ss_prev.str();
    return vbox({hbox({buttons->Render(), paragraph(s_prev), paragraph(s)})}) |
           flex | border;
  });

  auto screen = ScreenInteractive::FitComponent();
  screen.Loop(component);
  return 0;
}
