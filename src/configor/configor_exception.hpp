
#pragma once
#include <stdexcept>  // std::runtime_error, std::exception_ptr, std::rethrow_exception
#include <string>     // std::string

#ifndef CONFIGOR_ASSERT
#include <cassert>  // assert
#define CONFIGOR_ASSERT(...) assert(__VA_ARGS__)
#endif

namespace configor
{

//
// exceptions
//

class configor_exception : public std::runtime_error
{
public:
    explicit configor_exception(const char* message)
        : std::runtime_error(message)
    {
    }

    explicit configor_exception(const std::string& message)
        : std::runtime_error(message)
    {
    }
};

class configor_type_error : public configor_exception
{
public:
    explicit configor_type_error(const std::string& message)
        : configor_exception("config type error: " + message)
    {
    }
};

class configor_invalid_key : public configor_exception
{
public:
    explicit configor_invalid_key(const std::string& message)
        : configor_exception("invalid config key error: " + message)
    {
    }
};

class configor_invalid_iterator : public configor_exception
{
public:
    explicit configor_invalid_iterator(const std::string& message)
        : configor_exception("invalid config iterator error: " + message)
    {
    }
};

class configor_deserialization_error : public configor_exception
{
public:
    explicit configor_deserialization_error(const std::string& message)
        : configor_exception("config deserialization error: " + message)
    {
    }
};

class configor_serialization_error : public configor_exception
{
public:
    explicit configor_serialization_error(const std::string& message)
        : configor_exception("config serialization error: " + message)
    {
    }
};

//
// error handler
//

enum class error_policy
{
    strict = 1,  // throw exceptions
    record = 2,  // record error message
    ignore = 3,  // ignore all errors
};

class error_handler
{
public:
    virtual void handle(std::exception_ptr eptr) = 0;
};

template <error_policy _Policy>
class error_handler_with;

template <>
class error_handler_with<error_policy::strict> : public error_handler
{
public:
    virtual void handle(std::exception_ptr eptr) override
    {
        std::rethrow_exception(eptr);
    }
};

template <>
class error_handler_with<error_policy::ignore> : public error_handler
{
public:
    virtual void handle(std::exception_ptr eptr) override
    {
        // DO NOTHING
    }
};

template <>
class error_handler_with<error_policy::record> : public error_handler
{
public:
    virtual void handle(std::exception_ptr eptr) override
    {
        try
        {
            if (eptr)
            {
                std::rethrow_exception(eptr);
            }
        }
        catch (const configor_exception& e)
        {
            this->error = e.what();
        }
    }

    std::string error;
};

}  // namespace configor
