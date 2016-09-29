#ifndef TFracNum_h
#define TFracNum_h

#include <string>
#include <vector>

const char IOUTSTRING[3] = "%d";

/// \class PrimeFactors
/// \brief Prime factor decomposition of an unsigned integer
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

    /// erase
    map_type::iterator erase(map_type::const_iterator pos)
    { return Factors_.erase(pos); }

    /// explicit conversion to unsigned
    explicit operator unsigned() const;

    /// const access to factors map
    const map_type& factors() const
    { return Factors_; }
    
 private:

    map_type Factors_;
};

/// convert to string
std::string to_string(const PrimeFactors& pf);

/// \return prime number decomposition of product of arguments
PrimeFactors decompose(const std::vector<unsigned>& N);

/// \return prime number decomposition of product of arguments
template <typename ... Types>
PrimeFactors decompose(unsigned n, Types ... other_n)
{ return decompose({n, other_n...}); }

/// \return prime number decomposition of product of factorials of arguments
PrimeFactors decompose_factorial(const std::vector<unsigned>& N);

/// \return prime number decomposition of product of factorials of arguments
template <typename ... Types>
PrimeFactors decompose_factorial(unsigned n, Types ... other_n)
{ return decompose_factorial({n, other_n...}); }

/// equality operator
const bool operator==(const PrimeFactors& lhs, const PrimeFactors& rhs)
{ return lhs.factors() == rhs.factors(); }

/// remove factors with exponent = 0
PrimeFactors& remove_zeroes(PrimeFactors& PF);

/// remove common factors from arguments;
/// remove_zeroes is not called.
PrimeFactors remove_common(PrimeFactors& A, PrimeFactors& B);

/// remove common factors from arguments, storing into return value;
/// remove_zeroes is not called.
PrimeFactors common(PrimeFactors& A, PrimeFactors& B);

/// multiplication assignment
PrimeFactors& operator*=(PrimeFactors& lhs, const PrimeFactors& rhs)
{ for (const auto& p_e : rhs) lhs[p_e.first] += p_e.second; return lhs; }

/// multiplication
PrimeFactors operator*(PrimeFactors A, const PrimeFactors& B)
{ return A *= B; }

/// \return exponentiation of factorized unsigned
PrimeFactors pow(PrimeFactors pf, unsigned n)
{ for (auto& p_e : pf) p_e.second *= n; return pf; }

/// \return square root of factorized unsigned, throws if sqrt is not unsigned.
PrimeFactors sqrt(PrimeFactors pf);

/// factorize PrimeFactors into sqrt-able and non-sqrt-able portions;
/// note, no sqrt is applied
std::array<PrimeFactors, 2> factorize_sqrt(PrimeFactors pf);

/// \return sign factor
int sign_factor(int n, int d);

/// \class TFracNum
/// \brief Ractional number represented by prime decompositions of numerator and denominator
/// \author Jan.Friedrich@ph.tum.de, Daniel Greenwald
class TFracNum
{
public:

    /// constructor
    /// \param N exponents of the numerator's prime numbers
    /// \param D exponents of the denominator's prime numbers
    /// \param s sign of number (converted to +-1, 0, inf, or nan)
    TFracNum(const PrimeFactors& N, const PrimeFactors& D, double s = 1)
        : Numerator_(N), Denominator_(D), Sign_(std::isfinite(s) ? sign_of(s) : s)
    {
        remove_common(Numerator_, Denominator_);
        remove_zeroes(Numerator_);
        remove_zeroes(Denominator_);
    }

    /// int constructor
    /// \param N numerator
    /// \param D denominator
    TFracNum(int N, int D == 1)
        : TFracNum(decompose(abs(N)), decompose(abs(D)), sign_factor(N, D)) {}

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
    
    static TFracNum am0_to_J(unsigned J, int m, int m0)
    { return TFracNum(decompose(m0) * decompose_factorial(J + m, J - m), decompose_factorial(2 * J)); }

    static TFracNum c_sub_l(unsigned l)
    { return TFracNum(decompose(ell) * pow(decompose_factorial(ell), 2), decompose_factorial(2 * ell)); }
    
    static TFracNum cm_sub_l(unsigned l, int m)
    { return (l == 0) ? TFracNum(1) : TFracNum(decompose((l + m / 2)) * pow(decompose_factorial(l), 2), decompose_factorial(2 * l)); }
    
    static TFracNum cm_sub_l_2(unsigned l, int m)
    { return (l == 0) ? TFracNum(1) : TFracNum(decompose(l + m) * pow(decompose_factorial(l), 4), pow(decompose_factorial(2 * l), 2)); }

private:
    
    /// prime factorization of numerator
    PrimeFactors Numerator_;

    /// prime factorization of denominator
    PrimeFactors Denominator_;
    
    // Prefactor, including sign
    double Sign_{1};

public:

    const static TFracNum Zero;
    const static TFracNum One;
    const static TFracNum Two;
    const static TFracNum mTwo;
    const static TFracNum Quarter;


};

/// convert to string as [sign][N]/[D]
std::string to_string(const TFracNum& f);

/// convert to string as {N, D}
std::string to_header_string(const TFracNum& f);

/// convert to detailed string
std::string to_detailed_string(const TFracNum& f);

/// equality operator
const bool operator==(const TFracNum& lhs, const TFracNum& rhs)
{ return lhs.sign() == rhs.sign() and lhs.numerator() == rhs.numerator() and lhs.denominator() == rhs.denominator(); }

/// inequality operator
const bool operator!=(const TFracNum& lhs, const TFracNum& rhs)
{ return !(lhs == rhs); }

/// less than operator
const bool operator<(const TFracNum& lhs, const TFracNum& rhs);

/// \return whether TFracNum is zero
const bool is_zero(const TFracNum& f)
{ return f.sign() == 0; }

/// unary minus
TFracNum operator-(const TFracNum& f)
{ return TFracNum(f.numerator(), f.denominator(), -f.sign()); }

/// addition assignment
TFracNum& operator+=(TFracNum& lhs, const TFracNum& rhs);

/// addition
TFracNum operator+(TFracNum lhs, const TFracNum& rhs)
{ return lhs += rhs; }

/// multiplication assignment
TFracNum& operator*=(TFracNum& lhs, const TFracNum& rhs);

/// multiplication
TFracNum operator*(TFracNum lhs, const TFracNum& rhs)
{ return lhs *= rhs; }

/// multiple sign factors
double multiply_sign_factors(double s1, double s2);

/// invert sign
double invert_sign_factor(double s);

/// invert TFracNum
TFracNum invert(const TFracNum& f)
{ return TFracNum(f.denominator(), f.numerator(), invert_sign_factor(f.sign())); }

/// \return square root of absolute value of number with sign
/// preserved, potentially throws if sqrt is not int
TFracNum sqrt(const TFracNum& f)
{ return TFracNum(sqrt(f.numerator()), sqrt(f.denominator()), f.sign()); }

/// factorize rational number into cleanly sqrt-able and
/// non-cleanly-sqrt-able portions; note the sqrt is not taken of
/// either and the sign is transfered to the sqrt-able portion.
/// \return array := [sqrt-able, non-sqrt-able]
std::array<TFracNum, 2> factorize_sqrt(const TFracNum& f);

/// convert to sqrt string := [sign][N]/[D] sqrt([n][d]);
/// where clean sqrt is pulled out
std::string to_sqrt_string(const TFracNum& f)
{ return to_sqrt_string(factorize_sqrt(f)); }

/// convert to sqrt string := [sign][N_0]/[D_0] sqrt([N_1][D_1]);
std::string to_sqrt_string(const std::array<TFracNum, 2>& f);

/// \return exponentiated rational number with sign preserved
TFracNum pow(const TFracNum& f, unsigned n)
{ return TFracNum(pow(f.numerator(), n), pow(f.denominator(), n), f.sign()); }

/// \return absolute value of number
TFracNum abs(const TFracNum& f)
{ return TFracNum(f.numerator(), f.denominator(), std::abs(f.sign())); }

TFracNum sum_roots(const TFracNum& A, const TFracNum& B)
{ return pow(sqrt(A) + sqrt(B), 2); }

#endif
