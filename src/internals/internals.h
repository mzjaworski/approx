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
    struct dimension_data {
        Type current_coordinate = 0;
        Type starting_position = 0;
        Type stop_at = 0;
        Type step_size = 0;
        Type compensation = 0;
    };

    template <typename Type, std::enable_if_t<std::is_arithmetic_v<Type>, bool> = true>
    struct variable_integration_info {
        Type from = 0;
        Type to = 0;
        const unsigned int steps = 0;
    };

    template <typename ...Args>
    struct type_list_size;

    template <typename Head>
    struct type_list_size<Head>{
        static constexpr size_t value = 1;
    };

    template <typename Head, typename ...Tail>
    struct type_list_size<Head, Tail...>{
        static constexpr size_t value = 1 + type_list_size<Tail...>::value;
    };


    template <typename ...P>
    struct dummy {};

    template <template <typename ...> class Obj, typename T, typename ...P>
    struct internal;

    template <template <typename ...> class Obj, typename ...P1, typename T, typename L>
    struct internal<Obj, dummy<P1...>, T, L>
    {
        using type = Obj<P1..., T>;
    };

    template <template <typename ...> class Obj, typename ...P1, typename T>
    struct internal<Obj, dummy<P1...>, T>
    {
        using type = Obj<P1...>;
    };

    template <template <typename ...> class Obj, typename ...P1, typename T, typename ...P2>
    struct internal<Obj, dummy<P1...>, T, P2...>
    {
        using type = typename internal<Obj, dummy<P1..., T>, P2...>::type;
    };

    template <template <typename ...> class T, typename ...P>
    struct subst_all_but_last
    {
        using type = typename internal<T, dummy<>, P...>::type;
    };


    template <typename ...Args>
    struct check_if_all_are_arithmetic;

    template <typename Head>
    struct check_if_all_are_arithmetic<Head>{
        static constexpr bool value = std::is_arithmetic_v<Head>;
    };

    template <typename Head, typename ...Tail>
    struct check_if_all_are_arithmetic<Head, Tail...>{
        static constexpr bool value = std::conjunction_v<std::is_arithmetic<Head>, check_if_all_are_arithmetic<Tail...>>;
    };


    template <template <typename ...> class C, template <typename I> class W, typename T, typename ...P>
    struct make_internal;

    template <template <typename ...> class C, template <typename I> class W, typename ...E, typename H>
    struct make_internal<C, W, dummy<E...>, H>{
        using type = C<E..., W<H>>;
    };

    template <template <typename ...> class C, template <typename I> class W, typename ...E>
    struct make_internal<C, W, dummy<E...>>{
        using type = C<E...>;
    };

    template <template <typename ...> class C, template <typename I> class W, typename ...E, typename H, typename ...T>
    struct make_internal<C, W, dummy<E...>, H, T...>{
        using type = typename make_internal<C, W, dummy<E..., W<H>>, T...>::type;
    };

    template <template <typename I> class W, typename ...T>
    struct make_tuple_of{
        using type = typename make_internal<std::tuple, W, dummy<>, T...>::type;
    };


    // Apply normal summation for non floating point types
    template <typename T, std::enable_if_t<std::negation_v<std::is_floating_point<T>>, bool> = true>
    auto advance_coordinate(dimension_data<T>& input){
        input.current_coordinate += input.step_size;
    }

    // Apply Kahan summation algorithm for the floating point types
    template <typename T, std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
    auto advance_coordinate(dimension_data<T>& input){
        auto y = input.step_size - input.compensation;
        auto t = input.current_coordinate + y;
        input.compensation = (t - input.current_coordinate) - y;
        input.current_coordinate = t;
    }


    template <typename ...Args>
    struct advance_to_next_point;

    template <typename Head>
    struct advance_to_next_point<Head>{

        auto operator()(dimension_data<Head>& head){
            advance_coordinate(head);

            if (head.current_coordinate > head.stop_at){
                head.current_coordinate = head.starting_position;
                head.compensation = 0;
            }
        }
    };

    template <typename Head, typename ...Tail>
    struct advance_to_next_point<Head, Tail...>{

        auto operator()(dimension_data<Head>& head, dimension_data<Tail>& ...tail){
            advance_coordinate(head);

            if (head.current_coordinate > head.stop_at){
                head.current_coordinate = head.starting_position;
                head.compensation = 0;

                advance_to_next_point<Tail...>()(tail...);
            }
        }
    };


    template <class F, typename Head, typename ...Tail>
    struct bind_function_parameters{
        F f;

        auto operator()(const dimension_data<Head>& head, const dimension_data<Tail>& ...tail){
            return std::bind_front(std::forward<F>(f), head.current_coordinate, (tail.current_coordinate)...);
        }
    };

}

#endif
