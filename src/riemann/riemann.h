//          Copyright Mateusz Jaworski 2021 - 2021.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE.md or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#ifndef APPROX_RIEMANN_H
#define APPROX_RIEMANN_H

#include <numeric>
#include <span>

namespace mz::approx::riemann {

    namespace method {

        struct left_point: public mz::approx::internals::equal_greather_than {

            template <typename T>
            static constexpr T init(const T& from, const T& step_size){
                return from;
            }

            template <typename Output, typename Input, typename ...Inputs>
            static constexpr Output estimate_area(const std::tuple<std::tuple<Input, Inputs...>, Output>& centre,
                                                  const std::tuple<std::tuple<Input, Inputs...>, Output>& right){

                // bind parameters for convenience
                const auto& [c_inputs, c_output] = centre;
                const auto& [r_inputs, _] = right;

                // calculate difference between points
                const auto diff = mz::approx::internals::difference(r_inputs, c_inputs);

                // if difference of a first dimension is not equal to zero estimate area and return it
                return mz::approx::internals::first_entry_equals_to_zero(diff) ? 0 :
                    c_output * std::apply(mz::approx::internals::calculate_delta<Output, Input, Inputs...>, diff);
            }

        };

        struct mid_point: public mz::approx::internals::greather_than{

            template <typename T>
            static constexpr T init(const T& from, const T& step_size){
                return from + (step_size / 2);
            }

            template <typename Output, typename Input, typename ...Inputs>
            static constexpr Output estimate_area(const std::tuple<std::tuple<Input, Inputs...>, Output>& centre,
                                                  const std::tuple<std::tuple<Input, Inputs...>, Output>& right){

                // bind parameters for convenience
                const auto& [c_inputs, c_output] = centre;
                const auto& [r_inputs, r_output] = right;

                // calculate difference between points
                const auto diff = mz::approx::internals::difference(r_inputs, c_inputs);

                // estimate area and return it
                // TODO add some actual interpolation methods and a way to select them
                return mz::approx::internals::first_entry_equals_to_zero(diff) ? 0 : std::midpoint(c_output, r_output)
                    * std::apply(mz::approx::internals::calculate_delta<Output, Input, Inputs...>, diff);
            }

        };

        struct right_point: public mz::approx::internals::greather_than{

            template <typename T>
            static constexpr T init(const T& from, const T& step_size){
                return from + step_size;
            }

            template <typename Output, typename Input, typename ...Inputs>
            static constexpr Output estimate_area(const std::tuple<std::tuple<Input, Inputs...>, Output>& centre,
                                                  const std::tuple<std::tuple<Input, Inputs...>, Output>& right){

                // bind parameters for convenience
                const auto& [c_inputs, _] = centre;
                const auto& [r_inputs, r_output] = right;

                // calculate difference between points
                const auto diff = mz::approx::internals::difference(r_inputs, c_inputs);

                // estimate area and return it
                return mz::approx::internals::first_entry_equals_to_zero(diff) ? 0 : r_output
                    * std::apply(mz::approx::internals::calculate_delta<Output, Input, Inputs...>, diff);
            }

        };
    }

    template <typename Type>
    struct variable_integration_info{
        Type from = 0;
        Type to = 0;
        const unsigned long int steps = 0;
    };

    template <typename Method, typename Arg, typename ...Args, std::enable_if_t<mz::approx::internals::all_types_are_arithmetic<Arg, Args...>(), bool> = true>
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
            current_coordinate = starting_position = Method::init(from, step_size);
            return dimension_data;
        };

        constexpr auto calculate_delta = [](auto&& ...dimension_data){
            return (static_cast<double>(dimension_data.step_size) * ...);
        };

        const unsigned long int total_number_of_points = std::apply(calculate_total_number_of_points , info);
        using point_data_type = mz::approx::internals::make_tuple_of<mz::approx::internals::dimension_data, Arg, Args...>;

        // pick each entry from info tuple and input them into the initialize_dimension_data function together with a default constructed dimension_data struct
        point_data_type point_data = [&]<size_t... Is>(std::index_sequence<Is...>){
            return std::make_tuple(initialize_dimension_data(typename std::tuple_element<Is, point_data_type>::type(), std::get<Is>(info))... );
        }(std::make_index_sequence<std::tuple_size_v<point_data_type>>());

        // calculate n-dimensional delta
        const double delta = std::apply(calculate_delta, point_data);

        double result = 0;
        for (unsigned long int current_point = 0; current_point < total_number_of_points; current_point++){

            // bind point coordinates to the given function and evaluate it
            const auto output = std::apply(mz::approx::internals::bind_function_parameters<decltype(function), Arg, Args...>{function}, point_data)();

            // add slice area to the result
            result += output * delta;

            // advance point coordinates
            std::apply(mz::approx::internals::advance_to_next_point<Method, Arg, Args...>, point_data);
        }

        return result;
    }

    template <typename method, typename Input, typename Output, std::enable_if_t<mz::approx::internals::all_types_are_arithmetic<Input, Output>(), bool> = true>
    constexpr double approximate(const std::vector<std::tuple<std::tuple<Input>, Output>>& points){

        using point_type = std::tuple<std::tuple<Input>, Output>;

        // create a local copy of the points
        std::vector<point_type> local_data;
        local_data.reserve(points.size() + 2);

        // append copy of the first and the last element to ensure coherence of the algorithm
        local_data.insert(local_data.end(), points.begin(), points.end());
        local_data.push_back(points.back());

        // create vector views to avoid any further unnecessary copying
        auto equal_range = std::span(local_data.begin(), std::prev(local_data.end()));
        auto front_iter = std::next(local_data.begin(), 1);

        const auto estimate_area = [&front_iter](const auto& result, const auto& val){

            // if points are adjacent estimate area under the curve connecting them
            const auto output = mz::approx::internals::points_are_adjacent(val, *front_iter)
                    ? method::estimate_area(val, *front_iter) : 0;

            // advance secondary iterators
            front_iter++;

            // add output to the current result
            return result + output;
        };

        // calculate and return results
        return std::accumulate(equal_range.begin(), equal_range.end(), Output(0), estimate_area);
    }

}

#endif
