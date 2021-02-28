//          Copyright Mateusz Jaworski 2021 - 2021.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE.md or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#ifndef APPROX_INTERNALS_H
#define APPROX_INTERNALS_H

#include <type_traits>
#include <functional>
#include <iostream>
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

    // Kahan sum function usable as a binary predicate
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
    constexpr bool eq(T&& lhs, T&& rhs){
        return std::fabs(lhs - rhs) < std::numeric_limits<std::remove_cvref_t<T>>::epsilon();
    }

    template <typename T, std::enable_if_t<std::is_floating_point_v<std::remove_cvref_t<T>>, bool> = true>
    constexpr bool gt(T&& lhs, T&& rhs){
        return !eq(lhs,rhs) ? lhs > rhs : false;
    }

    template <typename T, std::enable_if_t<std::negation_v<std::is_floating_point<std::remove_cvref_t<T>>>, bool> = true>
    constexpr bool gt(T&& lhs, T&& rhs){
        return lhs > rhs;
    }

    template <typename Head>
    constexpr void advance_to_next_point(dimension_data<Head>& head){
        advance_coordinate(head);

        if (gt(head.current_coordinate, head.stop_at)){
            head.current_coordinate = head.starting_position;
            head.compensation = 0;
        }
    }

    template <typename Head, typename ...Tail>
    constexpr void advance_to_next_point(dimension_data<Head>& head, dimension_data<Tail>& ...tail){
        advance_coordinate(head);

        if (gt(head.current_coordinate, head.stop_at)){
            head.current_coordinate = head.starting_position;
            head.compensation = 0;

            advance_to_next_point(tail...);
        }
    }

    template <class F, typename Head, typename ...Tail>
    struct bind_function_parameters{
        F f;

        constexpr auto operator()(const dimension_data<Head>& head, const dimension_data<Tail>& ...tail){
            return std::bind_front(f, head.current_coordinate, (tail.current_coordinate)...);
        }
    };

}

#endif
