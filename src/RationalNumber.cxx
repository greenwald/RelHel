#include "RationalNumber.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string>

namespace relhel {

namespace prime_cache {

    // initialize with first primes up to 1000
    static std::vector<unsigned> primes_ = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53,
                                            59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
                                            127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181,
                                            191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251,
                                            257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317,
                                            331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397,
                                            401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463,
                                            467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557,
                                            563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619,
                                            631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701,
                                            709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787,
                                            797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863,
                                            877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953,
                                            967, 971, 977, 983, 991, 997};

    //-------------------------
    bool is_prime(unsigned p)
    {
        if (primes_.empty() and p == 2)
            return true;

        // check for upper limit of search range: p / 2
        auto last = std::upper_bound(primes_.begin(), primes_.end(), 0.5 * p);

        // search if any prime divides p
        for (auto it = primes_.begin(); it != last; ++it)
            if (p % *it == 0)
                return false;

        // if couldn't search high enough range, through
        if (last == primes_.end())
            throw;

        // else p is prime
        return true;
    }

    //-------------------------
    std::vector<unsigned>::iterator next(std::vector<unsigned>::iterator it)
    {
        // increment to next prime
        ++it;
        // if prime available, return it
        if (it != primes_.end())
            return it;
        // else add to cache, brute force
        // initialize p to next value
        unsigned p = primes_.empty() ? 2 : primes_.back() + 1;
        // increase until prime is found
        while (!is_prime(p)) ++p;
        // add to cache
        primes_.push_back(p);
        return --primes_.end();
    }
}

//-------------------------
PrimeFactors::operator unsigned() const
{
    return std::accumulate(Factors_.begin(), Factors_.end(), 1u,
                           [](unsigned a, const map_type::value_type& p_e)
                           { return a * std::pow(p_e.first, p_e.second); });
}

//-------------------------
std::string to_string(const PrimeFactors& pf)
{
    return is_one(pf) ? "1" : 
        std::accumulate(pf.begin(), pf.end(), std::string(),
                        [](std::string& s, const PrimeFactors::map_type::value_type& p_e)
                        { return s += " * " + exponential_string(std::to_string(p_e.first), p_e.second); }
            ).erase(0, 3);
}

//-------------------------
PrimeFactors decompose(const std::vector<unsigned>& N)
{
    PrimeFactors PF;
    for (unsigned n : N) {
        if (n == 0)
            PF[0] = 1;
        else {
            // check if in prime list
            auto pit = std::lower_bound(prime_cache::primes_.begin(), prime_cache::primes_.end(), n);
            if (pit != prime_cache::primes_.end() and *pit == n)
                PF[*pit] = 1;
            else {
                for (auto it = prime_cache::primes_.begin(); n > 1; it = prime_cache::next(it)) {
                    unsigned m = 0;
                    while (n % *it == 0) {
                        ++m;
                        n /= *it;
                    }
                    if (m > 0)
                        PF[*it] = m;
                }
            }
        }
    }
    return PF;
}

//-------------------------
PrimeFactors decompose_factorial(const std::vector<unsigned>& F)
{
    std::vector<unsigned> N;
    N.reserve(std::max(std::accumulate(F.begin(), F.end(), size_t(0)), F.size()) - F.size());
    for (unsigned f : F)
        while (f > 1)
            N.push_back(f--);
    return decompose(N);
}

//-------------------------
PrimeFactors& remove_zeroes(PrimeFactors& PF)
{
    for (auto it = PF.begin(); it != PF.end();)
        if (it->second == 0)
            it = PF.erase(it);
        else
            ++it;
    return PF;
}

//-------------------------
void remove_common(PrimeFactors& A, PrimeFactors& B)
{    
    for (auto& p_e : A) {
        auto it = B.find(p_e.first);
        if (it != B.end()) {
            unsigned c = std::min(p_e.second, it->second);
            p_e.second -= c;
            it->second -= c;
        }
    }
}

//-------------------------
PrimeFactors common(PrimeFactors& A, PrimeFactors& B)
{
    PrimeFactors C;
    for (auto& p_e : A) {
        auto it = B.find(p_e.first);
        if (it != B.end()) {
            unsigned c = std::min(p_e.second, it->second);
            C[p_e.first] = c;
            p_e.second -= c;
            it->second -= c;
        }
    }
    return C;
}

//-------------------------
PrimeFactors sqrt(PrimeFactors pf)
{
    for (auto& p_e : pf)
        if (is_even(p_e.second))
            p_e.second /= 2;
        else
            throw std::runtime_error("sqrt is not integral");
    return pf;
}

//-------------------------
std::array<PrimeFactors, 2> factorize_sqrt(PrimeFactors pf)
{
    PrimeFactors pf2;
    for (auto& p_e : pf)
        if (is_odd(p_e.second)) {
            pf2[p_e.first] = 1;
            --p_e.second;
        }
    return {pf, pf2};
}

//-------------------------
std::string to_string(const std::vector<long>& V)
{
    return std::accumulate(V.begin(), V.end(), std::string(),
                           [](std::string& s, long v)
                           {return s += " " + std::to_string(v);}).erase(0, 1);
}

//-------------------------
double sign_factor(int n, int d)
{
    if (n != 0 and d != 0)
        return sign_of(n * d);
    if (d != 0)
        return 0;
    if (n == 0)
        return std::numeric_limits<double>::quiet_NaN();
    return sign_of(n) * std::numeric_limits<double>::infinity();
}

//-------------------------
void remove_common(RationalNumber& A, RationalNumber& B)
{
    auto nA = A.numerator();
    auto nB = B.numerator();
    auto dA = A.denominator();
    auto dB = B.denominator();
    remove_common(nA, nB);
    remove_common(dA, dB);
    A = RationalNumber(nA, dA, A.sign());
    B = RationalNumber(nB, dB, B.sign());
}

//-------------------------
RationalNumber common(RationalNumber& A, RationalNumber& B)
{
    auto nA = A.numerator();
    auto nB = B.numerator();
    auto dA = A.denominator();
    auto dB = B.denominator();
    auto nC = common(nA, nB);
    auto dC = common(dA, dB);
    A = RationalNumber(nA, dA, A.sign());
    B = RationalNumber(nB, dB, B.sign());
    return RationalNumber(nC, dC);
}

//-------------------------
const bool operator<(const RationalNumber& lhs, const RationalNumber& rhs)
{
    // check by sign first
    if (lhs.sign() < rhs.sign())
        return true;

    if (lhs.sign() > rhs.sign())
        return false;

    // signs must be the same now:
    if (is_zero(lhs)) // therefore if both are zero
        return false;

    // remove common factors and compare uncommon factors
    auto lhs_temp = lhs;
    auto rhs_temp = rhs;
    remove_common(lhs_temp, rhs_temp);
    return static_cast<double>(lhs_temp) < static_cast<double>(rhs_temp);
}

//-------------------------
RationalNumber& operator+=(RationalNumber& lhs, const RationalNumber& rhs)
{
    if (rhs.sign() == 0 or !std::isfinite(lhs.sign()))
        return lhs;
    if (lhs.sign() == 0 or !std::isfinite(rhs.sign()))
        return lhs = rhs;

    // form new numerator
    auto LR = lhs.numerator() * rhs.denominator();
    auto RL = rhs.numerator() * lhs.denominator();
    // remove common factors
    auto C = common(LR, RL);

    // calculate sum of noncommon factors, including signs
    int num = lhs.sign() * static_cast<unsigned>(LR) + rhs.sign() * static_cast<unsigned>(RL);

    return lhs = RationalNumber(C * decompose(std::abs(num)), lhs.denominator() * rhs.denominator(), static_cast<double>(num));
}

//-------------------------
RationalNumber& operator*=(RationalNumber& lhs, const RationalNumber& rhs)
{
    return lhs = RationalNumber(lhs.numerator() * rhs.numerator(),
                          lhs.denominator() * rhs.denominator(),
                          multiply_sign_factors(lhs.sign(), rhs.sign()));
}

//-------------------------
double multiply_sign_factors(double s1, double s2)
{
    if ((s1 == 0 and std::isinf(s2)) or (s2 == 0 and std::isinf(s1)))
        return std::numeric_limits<double>::quiet_NaN();
    return s1 * s2;
}

//-------------------------
double invert_sign_factor(double s)
{
    if (s == 0)
        return std::numeric_limits<double>::infinity();
    if (std::isinf(s))
        return 0;
    return s;
}

//-------------------------
std::string to_string(const RationalNumber& f)
{
    if (f.sign() == 0)
        return "0";
    if (!std::isfinite(f.sign()))
        return std::to_string(f.sign());
    auto N = static_cast<unsigned>(f.numerator());
    auto D = static_cast<unsigned>(f.denominator());
    return (f.sign() > 0 ? "" : "-") + std::to_string(N) + (D == 1 ? "" : " / " + std::to_string(D));
}

std::string to_header_string(const RationalNumber& f)
{
    if (f.sign() == 0)
        return "{0,1}";
    if (std::isinf(f.sign()))
        return "{1,0}";
    if (std::isnan(f.sign()))
        return "{0,0}";
    return "{" + std::to_string(static_cast<unsigned>(f.numerator()) * f.sign()) 
        + "," + std::to_string(static_cast<unsigned>(f.denominator())) + "}";
}

//-------------------------
std::string to_detailed_string(const RationalNumber& f)
{
    return (std::isfinite(f.sign()) ? (f.sign() > 0 ? "" : "-") : (std::to_string(f.sign()) + " * "))
                    + "(" + to_string(f.numerator()) + ") / (" + to_string(f.denominator()) + ")";
}

//-------------------------
std::array<RationalNumber, 2> factorize_sqrt(const RationalNumber& f)
{
    // factorize numerator and denominator
    auto N = factorize_sqrt(f.numerator());
    auto D = factorize_sqrt(f.denominator());
    // and create new rational numbers
    return {RationalNumber(N[0], D[0], f.sign()), RationalNumber(N[1], D[1], 1)};
}

//-------------------------
std::string to_sqrt_string(const std::array<RationalNumber, 2>& f)
{
    if (is_zero(f[0]) or is_zero(f[1]))
        return "0";

    if (is_one(f[1]))
        return to_string(sqrt(f[0]));

    std::string s;
    if (is_one(f[0].numerator()) and is_one(f[0].denominator()))
        s = static_cast<double>(sqrt(f[0])) < 0 ? "-" : "";
    else
        s = to_string(sqrt(f[0])) + " * ";
    
    return s + "sqrt("
        + std::to_string(static_cast<unsigned>(f[1].numerator())) + " / "
        + std::to_string(static_cast<unsigned>(f[1].denominator())) + ")";
}


}
