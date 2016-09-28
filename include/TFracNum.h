#ifndef TFracNum_h
#define TFracNum_h

#include <string>
#include <vector>

const char IOUTSTRING[3] = "%d";

/// \typedef PrimeFactors
/// maps prime factor to exponent
using PrimeFactors = std::map<unsigned, unsigned>;

/// \return sign factor
int sign_factor(int n, int d);

/// \return prime number decomposition of unsigned integer
PrimeFactors decompose(unsigned n);

/// \return prime number decomposition of factorial of unsigned integer
PrimeFactors decompose_factorial(unsigned n);

/// \return multiplication assignment
PrimeFactors& operator*=(PrimeFactors& A, const PrimeFactors& B);

/// \class TFracNum
/// \brief Fractional Number.
///
/// Fractional number with numerator and denominator represented by
/// their prime number decomposition. Some arithmetical operations are
/// included but \b not complete.
///
/// \author Jan.Friedrich@ph.tum.de, Daniel Greenwald
class TFracNum
{

public:

    /// default constructor
    TFracNum() = default;

    /// vector constructor
    /// \param N exponents of the numerator's prime numbers
    /// \param D exponents of the denominator's prime numbers
    /// \param s sign of number
    TFracNum(const PrimeFactors& N, const PrimeFactors& D, double s);

    /// int constructor
    /// \param N numerator
    /// \param D denominator
    TFracNum(int N, int D)
        : TFracNum(decompose(abs(N)), decompose(abs(D)), sign_factor(N, D)) {}

    /// construct TFracNum from numerator and demoninator specified by factorials
    /// \param N numerator = N!
    /// \param D denimnator  = D!
    static TFracNum factorial_TFracNum(int N, int D)
    { return TFracNum(decompose_factorial(abs(N)), decompose_factorial(abs(D)), sign_factor(N, D)); }

    /// \return Numerator_
    const PrimeFactors& numerator() const
    { return Numerator_; }

    /// \return Denominator_
    const PrimeFactors& denominator() const
    { return Denominator_; }

    double sign() const
    { return Sign_; }



    
    //! Largest common divisor of numerator and denominator
    long DenomCommonDivisor(const TFracNum& rhs) const;

    //! The return value c satisfies ssqrt(c)=ssqrt(a)+ssqrt(b)
    //! Here the signed square root function ssqrt(n)=sign(n)*sqrt(abs(n)) is used
    TFracNum SumSignedRoots(const TFracNum& rhs) const;

    //! String in the form Num/Den. If Den=1, only Num is given
    const char* FracString() const;

    //! String of the square-root in the form <tt>Num/Den#RNum/RDen</tt>
    //! All numbers after the '#' are to be square-rooted, so Num/Den#RNum/RDen=Num/Den*sqrt(RNum/RDen)
    const char* FracStringSqrt() const;

    //! Complete information about the fractional number is put to cout
    std::ostream& Print(std::ostream& out) const;

    //! Return the double-precision real value
    const double& Dval(const bool& squareRootTheResult = false) const;

    //! Try square root operation
    /*! In case of success, return true. In case this does not lead to a
     fractional number, the number is left untouched and return value is false*/
    bool Sqrt();

    //! Force sign to plus
    bool Abs();

    //! Inversion of number, so nominator and denumerator are flipped
    bool Invert();

    //! Return sign as +1 (also for zero or undefined number) or -1
    long GetSign() const { return (Dval() < 0) ? -1 : 1; }

    //! Return numerator
    const long& GetNumerator() const;

    //! Return denominator
    const long& GetDenominator() const;

    //! Output some comparative values for two fractional numbers
    bool PrintDifference(const TFracNum&) const;

    //! String containing NOM and DEN
    const char* HeaderString() const;

    //! Check whether two fractional numbers are equal
    bool operator==(const TFracNum& rhs) const;
    //! Check whether two fractional numbers are not equal
    bool operator!=(const TFracNum& rhs) const { return not (*this == rhs); }
    //! Check whether left-hand number is greater than right-hand number
    bool operator>(const TFracNum& rhs) const;
    //! Multiply two fractional numbers
    TFracNum& operator*=(const TFracNum& rhs);
    //! Add two fractional numbers
    TFracNum& operator+=(const TFracNum& rhs);

    static TFracNum am0_to_J(const long& J, const long& m, const long& m0);
    static TFracNum c_sub_ell(const long& ell);
    static TFracNum cm0_sub_ell(const long& ell, const long& m0);
    static TFracNum cm0_sub_ell_2(const long& ell, const long& m0);

private:
    //
    // since Num is appearing as short form of "Number",
    // nom/NOM is taken when the numerator is meant
    //

    PrimeFactors Numerator_;
    PrimeFactors Denominator_;
    
    // Prime number decomposition of numerator. Field length is maxPrimNom,
    //  NOM[0] is the exponent of 2, NOM[1] of 3, and so on.
    std::vector<int> NumeratorExponents_;

    // Prime number decomposition of denominator, analogue to NOM
    std::vector<int> DenominatorExponents_;

    // Prefactor, including sign
    double SignPrefactor_{1};

    // Integers of numerator and denominator
    mutable int   _numerator_{0};
    mutable bool   _nomCacheRebuildRequired{true};
    mutable int   _denominator{1};
    mutable bool   _denCacheRebuildRequired{true};
    mutable double _value{0};
    mutable bool   _valueCacheRebuildRequired{true};

    void resetAllCaches() const;

    static int getNumberFromFactorization(const std::vector<int>& vector);

    static bool Debug_;

public:

    const static TFracNum Zero;
    const static TFracNum One;
    const static TFracNum Two;
    const static TFracNum mTwo;
    const static TFracNum Quarter;


};

/// unary minus
TFracNum operator-(const TFracNum& f)
{ return TFracNum(f.numerator(), f.denominator(), -f.sign()); }

/// multiple sign factors
double multiply_sign_factors(double s1, double s2);

/// addition assignment
TFracNum& operator+=(TFracNum& lhs, const TFracNum& rhs)
{
    return lhs = TFracNum(lhs.numerator() * rhs.denominator() + rhs.numerator() * lhs.denominator(),
                          lhs.denominator() * rhs.denominator(), multiply_sign_factors(lhs.sign(), rhs.sign()));
}

/// invert sign
double invert_sign_factor(double s);

/// invert TFracNum
TFracNum invert(const TFracNum& f)
{ return TFracNum(f.denominator(), f.numerator(), invert_sign_factor(f.sign())); }

inline
std::ostream&
operator <<(std::ostream&            out,
            const TFracNum&          fracNum)
{
    return fracNum.Print(out);
}


inline
TFracNum
operator *(TFracNum lhs, const TFracNum& rhs)
{
    lhs *= rhs;
    return lhs;
}


inline
TFracNum
operator +(TFracNum lhs, const TFracNum& rhs)
{
    lhs += rhs;
    return lhs;
}


#endif
