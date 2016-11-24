/// \file
/// \author Daniel Greenwald

#ifndef RationalNumber_h
#define RationalNumber_h

#include "MathUtils.h"

#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <string>
#include <vector>

namespace relhel {

/// \class PrimeFactors
/// \brief Prime factor decomposition of an unsigned integer
/// \author Daniel Greenwald
class PrimeFactors
{
 public:

    /// maps prime factor to exponent
    using map_type = std::map<unsigned, unsigned>;
    
    /// access operator
    map_type::mapped_type& operator[](const map_type::key_type& key)
    { return Factors_[key]; }

    /// access begin
    map_type::iterator begin()
    { return Factors_.begin(); }

    /// access begin
    map_type::const_iterator begin() const
    { return Factors_.begin(); }

    /// access end
    map_type::iterator end()
    { return Factors_.end(); }

    /// access end
    map_type::const_iterator end() const
    { return Factors_.end(); }

    const bool empty() const
    { return Factors_.empty(); }
    
    /// erase
    map_type::iterator erase(map_type::const_iterator pos)
    { return Factors_.erase(pos); }

    map_type::iterator find(const map_type::key_type& key)
    { return Factors_.find(key); }
    
    /// explicit conversion to unsigned
    explicit operator unsigned() const;

    /// const access to factors map
    const map_type& factors() const
    { return Factors_; }
    
 private:

    map_type Factors_;
};

/// check if equal to unity
inline const bool is_one(const PrimeFactors& pf)
{ return pf.empty(); }

/// convert to string
std::string to_string(const PrimeFactors& pf);

/// \return prime number decomposition of product of arguments
PrimeFactors decompose(const std::vector<unsigned>& N);

/// \return prime number decomposition of product of arguments
template <typename ... Types>
PrimeFactors decompose(unsigned n, Types ... other_n)
{ std::vector<unsigned> V({n, other_n...}); return decompose(V); }

/// \return prime number decomposition of product of factorials of arguments
PrimeFactors decompose_factorial(const std::vector<unsigned>& N);

/// \return prime number decomposition of product of factorials of arguments
template <typename ... Types>
PrimeFactors decompose_factorial(unsigned n, Types ... other_n)
{ std::vector<unsigned> V({n, other_n...}); return decompose_factorial(V); }

/// equality operator
inline const bool operator==(const PrimeFactors& lhs, const PrimeFactors& rhs)
{ return lhs.factors() == rhs.factors(); }

/// remove factors with exponent = 0
PrimeFactors& remove_zeroes(PrimeFactors& PF);

/// remove common factors from arguments;
/// remove_zeroes is not called.
void remove_common(PrimeFactors& A, PrimeFactors& B);

/// remove common factors from arguments, storing into return value;
/// remove_zeroes is not called.
PrimeFactors common(PrimeFactors& A, PrimeFactors& B);

/// multiplication assignment
inline PrimeFactors& operator*=(PrimeFactors& lhs, const PrimeFactors& rhs)
{ for (const auto& p_e : rhs) lhs[p_e.first] += p_e.second; return lhs; }

/// multiplication
inline PrimeFactors operator*(PrimeFactors A, const PrimeFactors& B)
{ return A *= B; }

/// \return exponentiation of factorized unsigned
inline PrimeFactors pow(PrimeFactors pf, unsigned n)
{ for (auto& p_e : pf) p_e.second *= n; return pf; }

/// \return square root of factorized unsigned, throws if sqrt is not unsigned.
PrimeFactors sqrt(PrimeFactors pf);

/// factorize PrimeFactors into sqrt-able and non-sqrt-able portions;
/// note, no sqrt is applied
std::array<PrimeFactors, 2> factorize_sqrt(PrimeFactors pf);

/// \return sign factor
double sign_factor(int n, int d);

/// \class RationalNumber
/// \brief Ractional number represented by prime decompositions of numerator and denominator
/// \author Daniel Greenwald
class RationalNumber
{
public:

    /// constructor
    /// \param N exponents of the numerator's prime numbers
    /// \param D exponents of the denominator's prime numbers
    /// \param s sign of number (converted to +-1, 0, inf, or nan)
    RationalNumber(const PrimeFactors& N, const PrimeFactors& D, double s = 1)
        : Numerator_(N), Denominator_(D), Sign_(std::isfinite(s) ? sign_of(s) : s)
    {
        remove_common(Numerator_, Denominator_);
        remove_zeroes(Numerator_);
        remove_zeroes(Denominator_);
    }

    /// int constructor
    /// \param N numerator
    /// \param D denominator
    explicit RationalNumber(int N, int D = 1)
        : RationalNumber(decompose(abs(N)), decompose(abs(D)), sign_factor(N, D)) {}

    /// \return Numerator_
    const PrimeFactors& numerator() const
    { return Numerator_; }

    /// \return Denominator_
    const PrimeFactors& denominator() const
    { return Denominator_; }

    /// \return Sign_
    const double sign() const
    { return Sign_; }

    /// convert to double
    explicit operator double() const
    { return Sign_ * static_cast<unsigned>(Numerator_) / static_cast<unsigned>(Denominator_); }

    /// flip the sign
    void negate()
    { Sign_ *= -1; }
    
private:
    
    /// prime factorization of numerator
    PrimeFactors Numerator_;

    /// prime factorization of denominator
    PrimeFactors Denominator_;
    
    // Prefactor, including sign
    double Sign_{1};

};

/// convert to string as [sign][N]/[D]
std::string to_string(const RationalNumber& f);

/// convert to string as {N, D}
std::string to_header_string(const RationalNumber& f);

/// convert to detailed string
std::string to_detailed_string(const RationalNumber& f);

/// equality operator
inline const bool operator==(const RationalNumber& lhs, const RationalNumber& rhs)
{ return lhs.sign() == rhs.sign() and lhs.numerator() == rhs.numerator() and lhs.denominator() == rhs.denominator(); }

/// inequality operator
inline const bool operator!=(const RationalNumber& lhs, const RationalNumber& rhs)
{ return !(lhs == rhs); }

/// less than operator
const bool operator<(const RationalNumber& lhs, const RationalNumber& rhs);

/// \return whether RationalNumber is zero
inline const bool is_zero(const RationalNumber& f)
{ return f.sign() == 0; }

/// \return whether RationalNumber is unity
inline const bool is_one(const RationalNumber& f)
{ return is_one(f.numerator()) and is_one(f.denominator()) and f.sign() == 1; }

/// remove common factors from arguments;
void remove_common(RationalNumber& A, RationalNumber& B);

/// remove common factors from arguments, storing into return value;
RationalNumber common(RationalNumber& A, RationalNumber& B);

/// unary minus
inline RationalNumber operator-(const RationalNumber& f)
{ return RationalNumber(f.numerator(), f.denominator(), -f.sign()); }

/// addition assignment
RationalNumber& operator+=(RationalNumber& lhs, const RationalNumber& rhs);

/// addition
inline RationalNumber operator+(RationalNumber lhs, const RationalNumber& rhs)
{ return lhs += rhs; }

/// multiplication assignment
RationalNumber& operator*=(RationalNumber& lhs, const RationalNumber& rhs);

/// multiplication
inline RationalNumber operator*(RationalNumber lhs, const RationalNumber& rhs)
{ return lhs *= rhs; }

/// multiplication assignment
inline RationalNumber& operator*=(RationalNumber& lhs, const PrimeFactors& rhs)
{ return lhs = RationalNumber(lhs.numerator() * rhs, lhs.denominator(), lhs.sign()); }

/// multiplication
inline RationalNumber operator*(RationalNumber lhs, const PrimeFactors& rhs)
{ return lhs *= rhs; }

/// multiplication
inline RationalNumber operator*(const PrimeFactors& lhs, RationalNumber rhs)
{ return rhs *= lhs; }

/// division of two PrimeFactors
inline RationalNumber operator/(const PrimeFactors& N, const PrimeFactors& D)
{ return RationalNumber(N, D); }

/// multiple sign factors
double multiply_sign_factors(double s1, double s2);

/// invert sign
double invert_sign_factor(double s);

/// invert RationalNumber
inline RationalNumber invert(const RationalNumber& f)
{ return RationalNumber(f.denominator(), f.numerator(), invert_sign_factor(f.sign())); }

/// \return square root of absolute value of number with sign
/// preserved, potentially throws if sqrt is not int
inline RationalNumber sqrt(const RationalNumber& f)
{ return RationalNumber(sqrt(f.numerator()), sqrt(f.denominator()), f.sign()); }

/// factorize rational number into cleanly sqrt-able and
/// non-cleanly-sqrt-able portions; note the sqrt is not taken of
/// either and the sign is transfered to the sqrt-able portion.
/// \return array := [sqrt-able, non-sqrt-able]
std::array<RationalNumber, 2> factorize_sqrt(const RationalNumber& f);

/// convert to sqrt string := [sign][N_0]/[D_0] sqrt([N_1][D_1]);
std::string to_sqrt_string(const std::array<RationalNumber, 2>& f);

/// convert to sqrt string := [sign][N]/[D] sqrt([n][d]);
/// where clean sqrt is pulled out
inline std::string to_sqrt_string(const RationalNumber& f)
{ return to_sqrt_string(factorize_sqrt(f)); }


/// \return exponentiated rational number with sign preserved
inline RationalNumber pow(const RationalNumber& f, unsigned n)
{ return RationalNumber(pow(f.numerator(), n), pow(f.denominator(), n), f.sign()); }

/// \return absolute value of number
inline RationalNumber abs(const RationalNumber& f)
{ return RationalNumber(f.numerator(), f.denominator(), std::abs(f.sign())); }

}
 
#endif
