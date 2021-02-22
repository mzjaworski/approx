//          Copyright Mateusz Jaworski 2021 - 2021.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE.md or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#ifndef APPROX_APPROX_H
#define APPROX_APPROX_H

#include "../../src/riemann/riemann.h"

namespace mz::approx {

    // approximate using riemann sums
    namespace riemann {

        namespace method{

            template <typename T>
            struct left_point;

            template <typename T>
            struct mid_point;

            template <typename T>
            struct right_point;

        }

        template <typename Type>
        struct variable_integration_info;

        template <template <typename T> typename method, typename Arg, typename ...Args, std::enable_if_t<mz::approx::internals::all_types_are_arithmetic<Arg, Args...>(), bool>>
        constexpr double approximate(const std::function<double(Arg, Args...)>& function,
                                     const mz::approx::internals::make_tuple_of<mz::approx::riemann::variable_integration_info, Arg,Args...>& info);

    }

    // TODO add different methods

}

#endif
