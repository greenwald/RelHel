#ifndef __MathUtils__h
#define __MathUtils__h

#include <cstdlib>
#include <string>
#include <vector>

constexpr const bool is_odd(int v) noexcept
{ return v & 0x1; }

constexpr const bool is_even(int v) noexcept
{ return !is_odd(v); }

inline std::string parity_to_string(int p)
{ return (p > 0) ? "+" : "-"; }

/// \return vector of spins from |s1-s2| to (s1 + s2)
const std::vector<unsigned> triangle(unsigned s1, unsigned s2);

/// \return sign of argument
template <typename T>
constexpr const T sign_of(T t) noexcept
{ return ((T(0) < t) - (t < T(0))); }

/// \return sign_of(val) * sqrt(abs(val))
constexpr const double signed_sqrt(double val) noexcept
{ return sign_of(val) * std::sqrt(std::abs(val)); }

/// \return exponent string := "" (n == 0), "s" (n == 1), "s^n" (otherwise)
std::string exponential_string(std::string s, int n)
{ return (n == 0) ? "" : (s + (n == 1) ? "" : "^" + std::to_string(n)); }

#endif
