#ifndef __MathUtils__h
#define __MathUtils__h

#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

namespace relhel {

/// \return whether value is odd
constexpr const bool is_odd(int v) noexcept
{ return v & 0x1; }

/// \return whether value is even
constexpr const bool is_even(int v) noexcept
{ return !is_odd(v); }

/// (-1)^n
constexpr int pow_negative_one(int n)
{ return is_odd(n) ? -1 : +1; }

/// \return '+' or '-' depending on sign of argument
inline std::string parity_to_string(int p)
{ return (p > 0) ? "+" : "-"; }

/// \return vector of spins from |s1-s2| to (s1 + s2)
const std::vector<unsigned> triangle(unsigned s1, unsigned s2);

/// \return vector of spin projections from -j to j
const std::vector<int> projections(unsigned j);

/// \return vector of spin projections of a collection of particles
const std::vector<std::vector<int> > projections(const std::vector<unsigned>& J);

/// \return sign of argument
template <typename T>
constexpr const T sign_of(T t) noexcept
{ return ((T(0) < t) - (t < T(0))); }

/// \return sign_of(val) * sqrt(abs(val))
constexpr const double signed_sqrt(double val) noexcept
{ return sign_of(val) * std::sqrt(std::abs(val)); }

/// \return exponent string := "" (n == 0), "s" (n == 1), "s^n" (otherwise)
inline std::string exponential_string(std::string s, int n)
{ return n == 0 ? "" : (s + (n == 1 ? "" : "^" + std::to_string(n))); }

/// \return spin as string
inline std::string spin_to_string(int s)
{ return is_even(s) ? std::to_string(s / 2) : (std::to_string(s) + "/2"); }

}

#endif
