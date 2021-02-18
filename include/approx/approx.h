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

            template <typename I, typename R>
            struct left_point;

            template <typename I, typename R>
            struct mid_point;

            template <typename I, typename R>
            struct right_point;

        }

        template <template <typename I, typename R> typename method, typename Arg, typename ...Args, std::enable_if_t<mz::approx::internals::check_if_all_are_arithmetic<Arg, Args...>::value, bool>>
        constexpr double approximate(const std::function<double(Arg, Args...)>& function,
                                     const typename mz::approx::internals::make_tuple_of<mz::approx::internals::variable_integration_info, Arg,Args...>::type& info);

    }

    // TODO add different methods

}



#endif
