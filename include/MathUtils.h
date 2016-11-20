#ifndef relhel__MathUtils_h
#define relhel__MathUtils_h

#include "Parity.h"

#include <cmath>
#include <cstdlib>
#include <stdexcept>
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

/// \return whether three spins fulfil triangle relationship
inline const bool triangle(unsigned two_j1, unsigned two_j2, unsigned two_j3)
{ return two_j3 >= std::abs((int)two_j1 - (int)two_j2) and two_j3 <= (two_j1 + two_j2); }

/// \return vector of (twice) spins from |two_s1-two_s2| to (two_s1 + two_s2)
const std::vector<unsigned> triangle(unsigned two_s1, unsigned two_s2);

/// \return vector of (twice) spin projections from -two_j to two_j
const std::vector<int> projections(unsigned two_j);

/// \return vector of spin projections of a collection of particles
const std::vector<std::vector<int> > projections(const std::vector<unsigned>& two_J);

/// \return sign of argument
template <typename T>
constexpr const T sign_of(T t) noexcept
{ return ((T(0) < t) - (t < T(0))); }

/// \return sign_of(val) * sqrt(abs(val))
inline const double signed_sqrt(double val) noexcept
{ return sign_of(val) * std::sqrt(std::abs(val)); }

/// \return exponent string := "" (n == 0), "s" (n == 1), "s^n" (otherwise)
inline std::string exponential_string(std::string s, int n)
{ return n == 0 ? "" : (s + (n == 1 ? "" : "^" + std::to_string(n))); }

/// \return spin as string
inline std::string spin_to_string(int two_s)
{ return is_even(two_s) ? std::to_string(two_s / 2) : (std::to_string(two_s) + "/2"); }

}

#endif
