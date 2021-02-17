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

    template <typename Type>
    struct dimension_data {
        Type current_coordinate;
        Type starting_position;
        Type stop_at;
        Type step_size;
    };

    template <typename Type, std::enable_if_t<std::is_arithmetic_v<Type>, bool> = true>
    struct variable_integration_info {
        Type from;
        Type to;
        const unsigned int steps;
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


    template <typename ...Args>
    struct advance_to_next_point;

    template <typename Head>
    struct advance_to_next_point<Head>{

        auto operator()(dimension_data<Head>& head){
            head.current_coordinate += head.step_size;

            if (head.current_coordinate >= head.stop_at){
                head.current_coordinate = head.starting_position;
            }
        }
    };

    template <typename Head, typename ...Tail>
    struct advance_to_next_point<Head, Tail...>{

        auto operator()(dimension_data<Head>& head, dimension_data<Tail>& ...tail){
            head.current_coordinate += head.step_size;

            if (head.current_coordinate >= head.stop_at){
                head.current_coordinate = head.starting_position;
                advance_to_next_point<Tail...>()(tail...);
            }
        }
    };


    template <typename ...Args>
    struct bind_function_parameters;

    template <class F>
    struct bind_function_parameters<F>{
        F f;

        auto operator()(){
            return f;
        }
    };

    template <class F, typename Head, typename ...Tail>
    struct bind_function_parameters<F, Head, Tail...>{
        F f;

        auto operator()(const dimension_data<Head>& head, const dimension_data<Tail>& ...tail){
            auto temp = std::bind_front(f, head.current_coordinate);
            return bind_function_parameters<decltype(temp), Tail...>{temp}(tail...);
        }
    };

}

#endif
