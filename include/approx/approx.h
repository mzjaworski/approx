//          Copyright Mateusz Jaworski 2021 - 2021.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE.md or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#ifndef APPROX_APPROX_H
#define APPROX_APPROX_H

#include "../../src/approx/internals.h"

namespace mz::approx {
    using namespace mz::approx::internals;

    template <typename Arg, typename ...Args, std::enable_if_t<check_if_all_are_arithmetic<Arg, Args...>::value, bool> = true>
    constexpr double integrate_midpoint_riemann(const std::function<double(Arg, Args...)>& function, const typename make_tuple_of<variable_integration_info, Arg,Args...>::type& info){

        const constexpr auto calculate_total_number_of_points = [](auto&& ...info_struct) {
            return (info_struct.steps * ...);
        };

        const auto initialize_dimension_data = [](auto&& dimension_data, const auto info_struct){
            auto& [current_coordinate, starting_position, stop_at, step_size] = dimension_data;
            const auto& [from, to, steps] = info_struct;
            stop_at = to;
            step_size = (to - from) / steps;
            current_coordinate = starting_position = step_size / 2;
            return dimension_data;
        };

        const auto calculate_delta = [](auto&& ...dimension_data){
            return (static_cast<double>(dimension_data.step_size) * ...);
        };

        const unsigned int total_number_of_points = std::apply(calculate_total_number_of_points , info);
        using point_data_type = typename make_tuple_of<dimension_data, Arg, Args...>::type;

        // pick entry from each tuple and input them into the initialize_dimension_data function
        point_data_type point_data = [&]<size_t... Is>(std::index_sequence<Is...>){
            return std::make_tuple(initialize_dimension_data(typename std::tuple_element<Is, point_data_type>::type(), std::get<Is>(info))... );
        }(std::make_index_sequence<std::tuple_size_v<point_data_type>>());

        // calculate n-dimensional delta size
        const double delta_size = std::apply(calculate_delta, point_data);

        double result = 0;
        for (size_t current_point = 0; current_point < total_number_of_points; current_point++){

            // bind point coordinates to the given function and evaluate it
            auto output = std::apply(bind_function_parameters<decltype(function), Arg, Args...>{function}, point_data)();

            // add output to the result
            result += output * delta_size;

            // advance point coordinates
            std::apply(advance_to_next_point<Arg, Args...>(), point_data);
        }

        return result;
    }

}


#endif
