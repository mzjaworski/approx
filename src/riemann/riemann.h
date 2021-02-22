//          Copyright Mateusz Jaworski 2021 - 2021.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE.md or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#ifndef APPROX_RIEMANN_H
#define APPROX_RIEMANN_H

#include "../internals/internals.h"

namespace mz::approx::riemann {

    namespace method {

        template <typename T>
        struct left_point{

            static constexpr T init(const T& from, const T& step_size){
                return from;
            }

        };

        template <typename T>
        struct mid_point{

            static constexpr T init(const T& from, const T& step_size){
                return from + (step_size / 2);
            }

        };

        template <typename T>
        struct right_point{

            static constexpr T init(const T& from, const T& step_size){
                return from + step_size;
            }

        };
    }

    template <typename Type>
    struct variable_integration_info{
        Type from = 0;
        Type to = 0;
        const unsigned int steps = 0;
    };

    template <template <typename T> typename method, typename Arg, typename ...Args, std::enable_if_t<mz::approx::internals::all_types_are_arithmetic<Arg, Args...>(), bool> = true>
    constexpr double approximate(const std::function<double(Arg, Args...)>& function,
                                 const mz::approx::internals::make_tuple_of<mz::approx::riemann::variable_integration_info, Arg,Args...>& info){

        constexpr auto calculate_total_number_of_points = [](auto&& ...info_struct) {
            return (info_struct.steps * ...);
        };

        constexpr auto initialize_dimension_data = []<typename T>(auto&& dimension_data, const mz::approx::riemann::variable_integration_info<T>& info_struct){
            auto& [current_coordinate, starting_position, stop_at, step_size, compensation] = dimension_data;
            auto [from, to, steps] = info_struct;

            // check whenever we have to swap integration range
            if (from > to)
                std::swap(from, to);

            stop_at = to;
            step_size = (to - from) / steps;
            current_coordinate = starting_position = method<T>::init(from, step_size);
            return dimension_data;
        };

        constexpr auto calculate_delta = [](auto&& ...dimension_data){
            return (static_cast<double>(dimension_data.step_size) * ...);
        };

        const unsigned int total_number_of_points = std::apply(calculate_total_number_of_points , info);
        using point_data_type = mz::approx::internals::make_tuple_of<mz::approx::internals::dimension_data, Arg, Args...>;

        // pick each entry from info tuple and input them into the initialize_dimension_data function together with a default constructed dimension_data struct
        point_data_type point_data = [&]<size_t... Is>(std::index_sequence<Is...>){
            return std::make_tuple(initialize_dimension_data(typename std::tuple_element<Is, point_data_type>::type(), std::get<Is>(info))... );
        }(std::make_index_sequence<std::tuple_size_v<point_data_type>>());

        // calculate n-dimensional delta
        const double delta = std::apply(calculate_delta, point_data);

        double result = 0;
        for (size_t current_point = 0; current_point < total_number_of_points; current_point++){

            // bind point coordinates to the given function and evaluate it
            const auto output = std::apply(mz::approx::internals::bind_function_parameters<decltype(function), Arg, Args...>{function}, point_data)();

            // add slice area to the result
            result += output * delta;

            // advance point coordinates
            std::apply(mz::approx::internals::advance_to_next_point<Arg, Args...>, point_data);
        }

        return result;
    }

}

#endif
