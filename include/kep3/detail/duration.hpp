#ifndef KEP3_DETAIL_DURATION_HPP
#define KEP3_DETAIL_DURATION_HPP

#include <chrono>

namespace kep3::detail
{

template <typename T>
struct is_duration : std::false_type {
};

template <typename Rep, typename Period>
struct is_duration<std::chrono::duration<Rep, Period>> : std::true_type {
};

} // namespace kep3::detail

#endif // KEP3_DETAIL_DURATION_HPP
