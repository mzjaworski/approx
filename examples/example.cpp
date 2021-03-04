//          Copyright Mateusz Jaworski 2021 - 2021.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE.md or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#include "../include/approx/approx.hpp"
#include <iostream>

int main() {

    const std::function function_1 = [](const int x) -> double { return x+1; };
    std::cout << "Result: " << mz::approx::riemann::approximate<mz::approx::riemann::method::mid_point>
            (function_1, {{0, 10, 5}}) << std::endl;

    const std::function function_2 = [](const double x) -> double { return x+1; };
    std::cout << "Result: " << mz::approx::riemann::approximate<mz::approx::riemann::method::left_point>
            (function_2, {{0, 10, 10000}}) << std::endl;

    const std::function function_3 = [](const double x, const double y) -> double { return sin(x)+cos(y); };
    std::cout << "Result: " << mz::approx::riemann::approximate<mz::approx::riemann::method::right_point>
            (function_3, {{0, 10, 100}, {0, 10, 100}}) << std::endl;

    const std::function function_4 = [](const double x, const int y) -> double { return sin(x)+y*50; };
    std::cout << "Result: " << mz::approx::riemann::approximate<mz::approx::riemann::method::mid_point>
            (function_4, {{0, 10, 100}, {0, 10, 5}}) << std::endl;

    const std::function function_5 = [](const double x, const int y, const char z, const double v) -> double { return sin(x)+y*50+cos(z*v); };
    std::cout << "Result: " << mz::approx::riemann::approximate<mz::approx::riemann::method::mid_point>
            (function_5, {{0, 10, 100}, {0, 20, 10}, {0, 100, 50}, {0, 10, 10}}) << std::endl;

    return 0;
}
