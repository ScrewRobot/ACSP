#ifndef ACSP_SAFE_CHECK_HPP
#define ACSP_SAFE_CHECK_HPP

#include <cmath>
#include<limits>

namespace ACSP::math
{
    template <typename T>
    inline auto is_positive_infinity(T val) -> bool
    {
        return val == std::numeric_limits<T>::infinity();
    }

    template <typename T>
    inline auto is_negative_infinity(T val) -> bool
    {
        return val == -std::numeric_limits<T>::infinity();
    }

    template <typename T>
    inline auto is_nan(T val) -> bool
    {
        return std::isnan(val);
    }

    template <typename T>
    inline auto is_safe(T val) -> bool
    {
        bool safe = true;
        safe &= !is_positive_infinity(val);
        safe &= !is_negative_infinity(val);
        safe &= !is_nan(val);

        return safe;

    }
}



#endif //ACSP_SAFE_CHECK_HPP
