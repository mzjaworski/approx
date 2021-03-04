//          Copyright Mateusz Jaworski 2021 - 2021.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE.md or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#ifndef TEST_APPROX_TRAPEZOIDAL_HPP
#define TEST_APPROX_TRAPEZOIDAL_HPP

#include <numeric>

namespace mz::approx::trapezoidal{

    template <typename Type>
    struct variable_integration_info{
        Type from = 0;
        Type to = 0;
        const unsigned long int steps = 0;
    };

    template <typename Arg, typename ...Args, std::enable_if_t<mz::approx::internals::all_types_are_arithmetic<Arg, Args...>(), bool> = true>
    constexpr double approximate(const std::function<double(Arg, Args...)>& function,
                                 const mz::approx::internals::make_tuple_of<mz::approx::trapezoidal::variable_integration_info, Arg,Args...>& info){

        constexpr auto calculate_total_number_of_points = [](auto&& ...info_struct) {
            return (info_struct.steps * ...);
        };

        constexpr auto get_number_of_steps_in_first_dimension = [](auto&& head, auto&& ...tail){
            return head.steps;
        };

        constexpr auto initialize_dimension_data = []<typename T>(auto&& dimension_data, const mz::approx::trapezoidal::variable_integration_info<T>& info_struct){
            auto& [current_coordinate, starting_position, stop_at, step_size, compensation] = dimension_data;
            auto [from, to, steps] = info_struct;

            // check whenever we have to swap integration range
            if (from > to)
                std::swap(from, to);

            stop_at = to;
            step_size = (to - from) / steps;
            current_coordinate = starting_position = from;
            return dimension_data;
        };

        constexpr auto calculate_delta = [](auto&& ...dimension_data){
            return (static_cast<double>(dimension_data.step_size) * ...);
        };

        const unsigned long int first_dimension_steps = std::apply( get_number_of_steps_in_first_dimension, info);
        const unsigned long int total_number_of_points = std::apply(calculate_total_number_of_points , info);

        using point_data_type = mz::approx::internals::make_tuple_of<mz::approx::internals::dimension_data, Arg, Args...>;

        // pick each entry from info tuple and input them into the initialize_dimension_data function together with a default constructed dimension_data struct
        point_data_type point_data = [&]<size_t... Is>(std::index_sequence<Is...>){
            return std::make_tuple(initialize_dimension_data(typename std::tuple_element<Is, point_data_type>::type(), std::get<Is>(info))... );
        }(std::make_index_sequence<std::tuple_size_v<point_data_type>>());

        // calculate n-dimensional delta
        const double delta = std::apply(calculate_delta, point_data);

        double result = 0.0;
        for (unsigned long int current_point = 0; current_point < total_number_of_points;){

            // get first output
            const auto output_1 = std::apply(mz::approx::internals::bind_function_parameters<decltype(function), Arg, Args...>{function}, point_data)();

            // advance point coordinates
            std::apply(mz::approx::internals::advance_to_next_point<mz::approx::internals::greather_than, Arg, Args...>, point_data);

            // get second output
            const auto output_2 = std::apply(mz::approx::internals::bind_function_parameters<decltype(function), Arg, Args...>{function}, point_data)();

            // check whenever we have advanced to the last point for the first variable if yes
            // advance again to point to the first position again
            if (!(++current_point % first_dimension_steps))
                std::apply(mz::approx::internals::advance_to_next_point<mz::approx::internals::greather_than, Arg, Args...>, point_data);

            // add outputs to the result
            result += ((output_1 + output_2) / 2) * delta;

        }

        return result;
    }

}

#endif
