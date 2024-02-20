#pragma once

namespace configor
{

enum class token_type
{
    uninitialized,

    literal_true,
    literal_false,
    literal_null,

    value_string,
    value_integer,
    value_float,

    begin_array,
    end_array,

    begin_object,
    end_object,

    name_separator,
    value_separator,

    end_of_input
};

inline const char* to_string(token_type token)
{
    switch (token)
    {
    case token_type::uninitialized:
        return "uninitialized";
    case token_type::literal_true:
        return "literal_true";
    case token_type::literal_false:
        return "literal_false";
    case token_type::literal_null:
        return "literal_null";
    case token_type::value_string:
        return "value_string";
    case token_type::value_integer:
        return "value_integer";
    case token_type::value_float:
        return "value_float";
    case token_type::begin_array:
        return "begin_array";
    case token_type::end_array:
        return "end_array";
    case token_type::begin_object:
        return "begin_object";
    case token_type::end_object:
        return "end_object";
    case token_type::name_separator:
        return "name_separator";
    case token_type::value_separator:
        return "value_separator";
    case token_type::end_of_input:
        return "end_of_input";
    }
    return "unknown";
}

}  // namespace configor
