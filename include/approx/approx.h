//          Copyright Mateusz Jaworski 2021 - 2021.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE.md or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#ifndef APPROX_APPROX_H
#define APPROX_APPROX_H

#include "../../src/internals/internals.h"
#include "../../src/riemann/riemann.h"
#include "../../src/trapezoidal/trapezoidal.h"

namespace mz::approx {

    // approximate using riemann sums
    namespace riemann {

        namespace method{

            struct left_point;

            struct mid_point;

            struct right_point;

        }

        template <typename Type>
        struct variable_integration_info;

        // used to approximate functions
        template <typename method, typename Arg, typename ...Args, std::enable_if_t<mz::approx::internals::all_types_are_arithmetic<Arg, Args...>(), bool>>
        constexpr double approximate(const std::function<double(Arg, Args...)>& function,
                                     const mz::approx::internals::make_tuple_of<mz::approx::riemann::variable_integration_info, Arg,Args...>& info);

        // used to approximate the area under a curve given by a vector of inputs->output tuples
        // points should be sorted in an ascending order in regards to their inputs, for now this function accepts only 1D inputs
        // I will expand it in the future, (it already contains some code to handle them)
        template <typename method, typename Input, typename Output, std::enable_if_t<mz::approx::internals::all_types_are_arithmetic<Input, Output>(), bool>>
        constexpr double approximate(const std::vector<std::tuple<std::tuple<Input>, Output>>& points);

    }

    // approximate using trapezoidal rule
    namespace trapezoidal {

        template <typename Type>
        struct variable_integration_info;

        // used to approximate functions
        template <typename Arg, typename ...Args, std::enable_if_t<mz::approx::internals::all_types_are_arithmetic<Arg, Args...>(), bool>>
        constexpr double approximate(const std::function<double(Arg, Args...)>& function,
                                     const mz::approx::internals::make_tuple_of<mz::approx::trapezoidal::variable_integration_info, Arg,Args...>& info);

    }

    // TODO add different methods

}



#endif
