#pragma once
#include "configor_basic.hpp"

namespace configor
{

using config  = basic_config<>;
using wconfig = basic_config<wconfig_args>;

template <typename _ConfTy, typename = typename std::enable_if<is_config<_ConfTy>::value>::type>
inline void swap(_ConfTy& lhs, _ConfTy& rhs)
{
    lhs.swap(rhs);
}

}  // namespace configor

namespace std
{
template <typename _Args>
struct hash<::configor::basic_config<_Args>>
{
    using argument_type = ::configor::basic_config<_Args>;
    using result_type   = size_t;

    result_type operator()(argument_type const& config) const
    {
        return hash<typename argument_type::string_type>{}(config.dump());
    }
};
}  // namespace std
