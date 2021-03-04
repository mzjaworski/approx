//          Copyright Mateusz Jaworski 2021 - 2021.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE.md or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#ifndef APPROX_INTERNALS_HPP
#define APPROX_INTERNALS_HPP

#include <type_traits>
#include <functional>
#include <iostream>
#include <numeric>
#include <cmath>
#include <tuple>

namespace mz::approx::internals {

    template <typename ...Args>
    void debug_print(Args&& ...args){
        ((std::cout << args << " "), ...) << std::endl;
    }

    template <typename Type>
    struct dimension_data{
        Type current_coordinate = 0;
        Type starting_position = 0;
        Type stop_at = 0;
        Type step_size = 0;
        Type compensation = 0;
    };

    template <typename ...Tail>
    constexpr bool all_types_are_arithmetic(){
        return (std::is_arithmetic_v<Tail> && ...);
    }

    template <typename ...>
    struct dummy{};

    template <template <typename ...> class C, template <typename> class W, typename T, typename ...P>
    struct make_internal;

    template <template <typename ...> class C, template <typename> class W, typename ...E, typename H>
    struct make_internal<C, W, dummy<E...>, H>{
        using type = C<E..., W<H>>;
    };

    template <template <typename ...> class C, template <typename> class W, typename ...E>
    struct make_internal<C, W, dummy<E...>>{
        using type = C<E...>;
    };

    template <template <typename ...> class C, template <typename> class W, typename ...E, typename H, typename ...T>
    struct make_internal<C, W, dummy<E...>, H, T...>{
        using type = typename make_internal<C, W, dummy<E..., W<H>>, T...>::type;
    };

    // Kahan sum function usable as a binary predicate in eg. accumulate
    template <typename T, typename I>
    constexpr std::tuple<T,T> kahan_sum(const std::tuple<T,T>& sum_compensation, const I& next){
        auto [sum, compensation] = sum_compensation;
        const auto y = static_cast<T>(next) - compensation;
        const auto t = sum + y;
        compensation = (t - sum) - y;
        sum = t;
        return {sum, compensation};
    }

    template <typename T, std::enable_if_t<std::is_arithmetic_v<std::remove_cvref_t<T>>, bool> = true>
    constexpr bool eq(const T& lhs, const T& rhs){
        return std::fabs(lhs - rhs) < std::numeric_limits<std::remove_cvref_t<T>>::epsilon();
    }

    template <typename T, std::enable_if_t<std::is_floating_point_v<std::remove_cvref_t<T>>, bool> = true>
    constexpr bool gt(const T& lhs, const T& rhs){
        return !eq(lhs,rhs) ? lhs > rhs : false;
    }

    template <typename T, std::enable_if_t<std::negation_v<std::is_floating_point<std::remove_cvref_t<T>>>, bool> = true>
    constexpr bool gt(const T& lhs, const T& rhs){
        return lhs > rhs;
    }

    template <typename T, std::enable_if_t<std::is_floating_point_v<std::remove_cvref_t<T>>, bool> = true>
    constexpr bool egt(const T& lhs, const T& rhs){
        return eq(lhs,rhs) ? true : lhs > rhs;
    }

    struct greather_than{

            template <typename T>
            static constexpr auto compare = gt<T>;

    };

    struct equal_greather_than{

        template <typename T>
        static constexpr auto compare = egt<T>;

    };

    template <template <typename> class W, typename ...T>
    using make_tuple_of = typename make_internal<std::tuple, W, dummy<>, T...>::type;

    // Apply normal summation for non floating point types
    template <typename T, std::enable_if_t<std::negation_v<std::is_floating_point<T>>, bool> = true>
    constexpr void advance_coordinate(dimension_data<T>& input){
        input.current_coordinate += input.step_size;
    }

    // Apply Kahan summation algorithm for the floating point types
    template <typename T, std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
    constexpr void advance_coordinate(dimension_data<T>& input){
        const auto y = input.step_size - input.compensation;
        const auto t = input.current_coordinate + y;
        input.compensation = (t - input.current_coordinate) - y;
        input.current_coordinate = t;
    }

    template <typename Comparator>
    constexpr void advance_to_next_point(){}

    template <typename Comparator, typename Head, typename ...Tail>
    constexpr void advance_to_next_point(dimension_data<Head>& head, dimension_data<Tail>& ...tail){
        advance_coordinate(head);

        if (Comparator::template compare<Head>(head.current_coordinate, head.stop_at)){
            head.current_coordinate = head.starting_position;
            head.compensation = 0;

            advance_to_next_point<Comparator>(tail...);
        }
    }

    template <class F, typename Head, typename ...Tail>
    struct bind_function_parameters{
        F f;

        constexpr auto operator()(const dimension_data<Head>& head, const dimension_data<Tail>& ...tail){
            return std::bind_front(f, head.current_coordinate, (tail.current_coordinate)...);
        }
    };

    template <typename F, size_t Index>
    constexpr auto zip_with_helper(F&& f, const auto&... i){
        return f(std::get<Index>(i)...);
    }

    template <typename F, template <typename ...> typename ...T, typename ...I>
    constexpr auto zip_with(F&& f, const T<I...>& ...args){
        return [&]<size_t ...V>(std::integer_sequence<size_t, V...> ...){
            return std::make_tuple(zip_with_helper<F, V>(std::forward<F>(f), args...)...);
        }(std::make_index_sequence<std::tuple_size_v<T<I...>>>()...);
    }

    template <typename ...T>
    constexpr std::tuple<T...> difference(const std::tuple<T...>& lhs, const std::tuple<T...>& rhs){
        return zip_with(std::minus<void>(), lhs, rhs);
    }

    template <typename ...T>
    constexpr std::tuple<T...> midpoint(const std::tuple<T...>& lhs, const std::tuple<T...>& rhs){

        const auto calculate_midpoint = [](auto&& lhs, auto&& rhs){
            return std::midpoint(lhs,rhs);
        };

        return zip_with(calculate_midpoint, lhs, rhs);
    }

    template <typename ...T>
    constexpr bool has_negative_entry(const std::tuple<T...>& input){

        const auto is_negative = [](auto&& ...args){
            return ((args < 0) || ...);
        };

        return std::apply(is_negative, input);
    }

    template <typename Output, typename ...Inputs>
    constexpr Output calculate_delta(Inputs ...inputs){
        return ((gt(static_cast<Output>(inputs), static_cast<Output>(0)) ? static_cast<Output>(inputs) : 1) * ...);
    }

    template <typename Output, typename Input, typename ...Inputs>
    constexpr bool points_are_adjacent(const std::tuple<std::tuple<Input, Inputs...>, Output>& lhs,
                                       const std::tuple<std::tuple<Input, Inputs...>, Output>& rhs){

        const auto& [lhs_i, _] = lhs;
        const auto& [rhs_i, __] = rhs;

        return !has_negative_entry(difference(rhs_i, lhs_i));
    }

    template <typename ...T>
    constexpr bool first_entry_equals_to_zero(const std::tuple<T...>& input){

        const auto is_zero = []<typename I>(I&& head, auto&& ...tail){
            return eq(head, static_cast<I>(0));
        };

        return std::apply(is_zero, input);
    }

}

#endif
