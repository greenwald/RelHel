#include "TFracNum.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <stdio.h>
#include <string>

namespace prime_cache {

    // initialize with first primes up to 1000
    static std::vector<int> primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
                                      31, 37, 41, 43, 47, 53, 59, 61, 67,
                                      71, 73, 79, 83, 89, 97, 101, 103, 107,
                                      109, 113, 127, 131, 137, 139, 149, 151,
                                      157, 163, 167, 173, 179, 181, 191, 193,
                                      197, 199, 211, 223, 227, 229, 233, 239,
                                      241, 251, 257, 263, 269, 271, 277, 281,
                                      283, 293, 307, 311, 313, 317, 331, 337,
                                      347, 349, 353, 359, 367, 373, 379, 383,
                                      389, 397, 401, 409, 419, 421, 431, 433,
                                      439, 443, 449, 457, 461, 463, 467, 479,
                                      487, 491, 499, 503, 509, 521, 523, 541,
                                      547, 557, 563, 569, 571, 577, 587, 593,
                                      599, 601, 607, 613, 617, 619, 631, 641,
                                      643, 647, 653, 659, 661, 673, 677, 683,
                                      691, 701, 709, 719, 727, 733, 739, 743,
                                      751, 757, 761, 769, 773, 787, 797, 809,
                                      811, 821, 823, 827, 829, 839, 853, 857,
                                      859, 863, 877, 881, 883, 887, 907, 911,
                                      919, 929, 937, 941, 947, 953, 967, 971,
                                      977, 983, 991, 997};

    //-------------------------
    bool is_prime(unsigned p)
    {
        // check for upper limit of search range: p / 2
        auto last = std::upper_bound(primes.begin(), primes.end(), 0.5 * p);
        // search if any prime divides p
        for (auto it = primes.begin(); it != last; ++it)
            if (p % *it == 0)
                return false;
        // if couldn't search high enough range, through
        if (last == primes.end())
            throw;
        // else p is prime
        return true;
    }

    //-------------------------
    std::vector<int>::iterator next(std::vector<int>::iterator it)
    {
        // increment to next prime
        ++it;
        // if prime available, return it
        if (it != primes.end())
            return it;
        // else add to cache, brute force
        // initialize p to next value
        int p = primes.back() + 1;
        // increase until prime is found
        while (!is_prime(p)) ++p;
        // add to cache
        primes.push_back(p);
        return --primes.end();
    }
}

//-------------------------
PrimeFactors decompose(unsigned N)
{
    PrimeFactors PF;
    for (auto it = prime_cache::primes.begin(); N > 1; it = prime_cache::next(it)) {
        unsigned n = 0;
        while (N % *it == 0) {
            ++n;
            N /= *it;
        }
        if (n > 0)
            PrimeFactors[*it] == n;
    }
    return PF;
}

//-------------------------
PrimeFactors decompose_factorial(unsigned n)
{
    std::vector<PrimeFactors> PFs;
    PFs.reserve(n - 1);
    while (n > 1) PFs.push_back(decompose(n--));
    return std::accumulate(PFs.begin(), PFs.end(), PrimeFactors(),
                           [](PrimeFactors& PF, const PrimeFactors& pf)
                           { return PF *= pf; });
}

//-------------------------
PrimeFactors& operator*=(PrimeFactors& A, const PrimeFactors& B)
{
    for (const auto& p_e : B)
        A[p_e.first] += p_e.second;
    return A;
}

bool TFracNum::Debug_ = false;

const TFracNum TFracNum::Zero = TFracNum(0, 1);
const TFracNum TFracNum::One = TFracNum(1, 1);
const TFracNum TFracNum::Two = TFracNum(2, 1);
const TFracNum TFracNum::mTwo = TFracNum(-2, 1);
const TFracNum TFracNum::Quarter = TFracNum(1, 4);

const char* SQUAREROOT_CHAR = "#";

//-------------------------
std::string to_string(const std::vector<long>& V)
{
    return std::accumulate(V.begin(), V.end(), std::string(),
                           [](std::string& s, long v)
                           {return s += " " + std::to_string(v);}).erase(0, 1);
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
TFracNum::TFracNum(const PrimeFactors& N, const PrimeFactors& D, double s)
    : NumeratorExponents_(N), DenominatorExponents_(D), Sign_(s)
{
    // cancel common terms
    for (auto& p_e : NumeratorExponents_) {
        auto it = DenominatorExponents_.find(p_e.first);
        if (it != DenominatorExponents_.end()) {
            auto common = std::min(it->second, p_e.second);
            p_e.second -= common;
            it->second -= common;
        }
    }
    remove_zeroes(NumeratorExponents_);
    remove_zeroes(DenominatorExponents_);
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
    return std::numeric_limits<double>::infinity();
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

const long& TFracNum::GetNumerator() const
{
    if (_nomCacheRebuildRequired) {
        if (_signPrefac == 0) {
            _numerator = 0;
        } else {
            _numerator = getNumberFromFactorization(_NOM);
        }
        _nomCacheRebuildRequired = false;
    }
    return _numerator;
}


const long& TFracNum::GetDenominator() const
{
    if (_denCacheRebuildRequired) {
        if (_signPrefac == 0) {
            _denominator = 1;
        } else {
            _denominator = getNumberFromFactorization(_DEN);
        }
        _denCacheRebuildRequired = false;
    }
    return _denominator;
}


const double& TFracNum::Dval(const bool& squareRootTheResult) const
{
    if (_valueCacheRebuildRequired) {
        if (_signPrefac == 0) {
            _value = 0.;
        } else {
            double absVal = ((double)GetNumerator()) / ((double)GetDenominator());
            if (squareRootTheResult) {
                absVal = std::sqrt(absVal);
            }
            _value = ((double)_signPrefac) * absVal;
        }
        _valueCacheRebuildRequired = false;
    }
    return _value;
}


void TFracNum::resetAllCaches() const
{
    _numerator   = 0;
    _denominator = 0;
    _value       = 0.;
    _nomCacheRebuildRequired   = true;
    _denCacheRebuildRequired   = true;
    _valueCacheRebuildRequired = true;
}


long TFracNum::DenomCommonDivisor(const TFracNum& rhs) const
{
    size_t minPD = min(_DEN.size(), rhs._DEN.size());
    rpwa::primeNumbers::entryType comdiv = 1;
    for (size_t i = 0; i < minPD; i++) {
        long ppot = _DEN[i];
        if (rhs._DEN[i] < ppot) {
            ppot = rhs._DEN[i];
        }
        while (ppot-- > 0) {
            comdiv *= rpwa::primeNumbers::instance().primeNumber(i);
        }
    }
    return comdiv;
}


TFracNum
TFracNum::SumSignedRoots(const TFracNum& rhs) const
{
    TFracNum mixed = (*this) * rhs;
    TFracNum aa = *this;
    TFracNum bb = rhs;
    if (mixed.Sqrt()) {
        bool flipsign = (aa.Dval() + bb.Dval() < 0);
        aa.Abs();
        bb.Abs();
        TFracNum res = aa + bb + TFracNum::Two * mixed;
        if (flipsign) {
            res.FlipSign();
        }
        return res;
    }
    std::cout << "Error in TFracNum::SumSignedRoots()" << endl << "this:" << Dval() << endl
              << *this << endl
              << "b:" << rhs.Dval() << endl
              << rhs;
    throw;
}


bool TFracNum::Sqrt()
{
    if (_signPrefac == 0 or GetNumerator() == 0) {
        return true;
    }
    // TODO: move this to the loops further down (or remove it completely?)
    if (_debug) {
        long sqrt_ok = 1;
        for (size_t i = 0; i < _NOM.size(); i++) {
            if (_NOM[i] % 2) {
                sqrt_ok = 0;
                break;
            }
        }
        if (sqrt_ok == 1) {
            for (size_t i = 0; i < _DEN.size(); i++) {
                if (_DEN[i] % 2) {
                    sqrt_ok = 0;
                    break;
                }
            }
        }
        if (sqrt_ok == 0) {
            cout << "square root not possible for this fracnum :(" << endl;
            cout << *this;
        }
    }
    for (size_t i = 0; i < _NOM.size(); i++) {
        if (_NOM[i] % 2) {
            return false;
        }
    }
    for (size_t i = 0; i < _DEN.size(); i++) {
        if (_DEN[i] % 2) {
            return false;
        }
    }
    for (size_t i = 0; i < _NOM.size(); i++) {
        _NOM[i] /= 2;
    }
    for (size_t i = 0; i < _DEN.size(); i++) {
        _DEN[i] /= 2;
    }
    resetAllCaches();
    return true;
}

bool TFracNum::Abs()
{
    if (_signPrefac == 0) {
        return true;
    }
    _signPrefac = 1;
    if (not _nomCacheRebuildRequired) {
        _numerator = labs(_numerator);
    }
    if (not _valueCacheRebuildRequired) {
        _value     = fabs(_value);
    }
    return true;
}


bool TFracNum::Invert()
{
    if (_signPrefac == -7777) {
        _NOM = vector<long>();
        _DEN = vector<long>();
        _signPrefac = -6666;
        resetAllCaches();
        return false;
    }
    if (GetNumerator() == 0) {
        _NOM = vector<long>();
        _DEN = vector<long>();
        _signPrefac = -7777;
        resetAllCaches();
        return false;
    }
    vector<long> oldNOM = _NOM;
    _NOM = _DEN;
    _DEN = oldNOM;
    resetAllCaches();
    return true;
}


bool TFracNum::operator==(const TFracNum& b) const
{
    if (_signPrefac == 0 && b._signPrefac == 0) {
        return true;
    }
    if (_signPrefac != b._signPrefac) {
        return false;
    }
    if (_NOM.size() != b._NOM.size()) {
        return false;
    }
    if (_DEN.size() != b._DEN.size()) {
        return false;
    }
    for (size_t i = 0; i < _NOM.size(); i++) {
        if (_NOM[i] != b._NOM[i]) {
            return false;
        }
    }
    for (size_t i = 0; i < _DEN.size(); i++) {
        if (_DEN[i] != b._DEN[i]) {
            return false;
        }
    }
    return true;
}


bool TFracNum::PrintDifference(const TFracNum& b) const
{
    if (_signPrefac == 0 && b._signPrefac == 0) {
        cout << "Both zero, they are equal." << endl;
        return true;
    }
    if (_signPrefac != b._signPrefac) {
        cout << "Different sign: " << _signPrefac << "!=" << b._signPrefac
             << endl;
        return false;
    }
    if (_NOM.size() != b._NOM.size()) {
        cout << "Different maxPrimNom: " << _NOM.size() << "!=" << b._NOM.size()
             << endl;
        return false;
    }
    if (_DEN.size() != b._DEN.size()) {
        cout << "Different maxPrimDen: " << _DEN.size() << "!=" << b._DEN.size()
             << endl;
        return false;
    }
    for (size_t i = 0; i < _NOM.size(); i++) {
        if (_NOM[i] != b._NOM[i]) {
            cout << "Different numerator contribution at prime " << i << ": "
                 << _NOM[i] << "!=" << b._NOM[i] << endl;
            return false;
        }
    }
    for (size_t i = 0; i < _DEN.size(); i++) {
        if (_DEN[i] != b._DEN[i]) {
            cout << "Different denominator contribution at prime " << i << ": "
                 << _DEN[i] << "!=" << b._DEN[i] << endl;
            return false;
        }
    }
    cout << "Well, they're simply equal!" << endl;
    return true;
}


const char*
TFracNum::HeaderString() const
{
    char* hstr = new char[30];
    if (_signPrefac == 0) {
        sprintf(hstr, "{0,1}");
        return hstr;
    }
    if (_signPrefac == -7777) {
        sprintf(hstr, "{1,0}");
        return hstr;
    }
    if (_signPrefac == -6666) {
        sprintf(hstr, "{0,0}");
        return hstr;
    }
    if (_signPrefac == 1) {
        sprintf(hstr, "{%ld,%ld}", GetNumerator(), GetDenominator());
    } else {
        sprintf(hstr, "{%ld,%ld}", -GetNumerator(), GetDenominator());
    }
    return hstr;
}


bool TFracNum::operator>(const TFracNum& rhs) const
{
    if (Dval() > rhs.Dval()) {
        return true;
    }
    return false;
}



TFracNum& TFracNum::operator+=(const TFracNum& rhs)
{
    long den_cdiv = DenomCommonDivisor(rhs);
    long bdc = rhs.GetDenominator() / den_cdiv;
    long adc = GetDenominator() / den_cdiv;
    *this = TFracNum(_signPrefac * GetNumerator() * bdc + rhs._signPrefac * rhs.GetNumerator() * adc, GetDenominator() * bdc);
    return *this;
}


TFracNum& TFracNum::operator*=(const TFracNum& rhs)
{

    // if one of the two numbers is undetermined,
    // the product is also undetermined
    if (_signPrefac == -6666 or rhs._signPrefac == -6666) {
        *this = TFracNum(vector<long>(), vector<long>(), -6666);
        return *this;
    }

    // if one of the two numbers contains division by zero,
    // and the other nominator is zero, the product is undetermined
    if ((_signPrefac == -7777 and rhs._signPrefac == 0) or
            (_signPrefac == 0     and rhs._signPrefac == -7777)) {
        *this = TFracNum(vector<long>(), vector<long>(), -6666);
        return *this;
    }

    // other cases with division by zero; product is also infinity
    if ((_signPrefac == -7777 or rhs._signPrefac == -7777)) {
        *this = TFracNum(vector<long>(), vector<long>(), -7777);
        return *this;
    }

    if (_signPrefac * rhs._signPrefac == 0) {
        *this = TFracNum::Zero;
        return *this;
    }

    if (_NOM.size() < rhs._NOM.size()) {
        _NOM.resize(rhs._NOM.size(), 0);
    }
    for (size_t i = 0; i < rhs._NOM.size(); ++i) {
        _NOM[i] += rhs._NOM[i];
    }

    if (_DEN.size() < rhs._DEN.size()) {
        _DEN.resize(rhs._DEN.size(), 0);
    }
    for (size_t i = 0; i < rhs._DEN.size(); ++i) {
        _DEN[i] += rhs._DEN[i];
    }

    for (size_t i = 0; i < min(_NOM.size(), _DEN.size()); ++i) {
        if (_NOM[i] != 0 and _DEN[i] != 0) {
            if (_DEN[i] > _NOM[i]) {
                _DEN[i] -= _NOM[i];
                _NOM[i] = 0;
            } else {
                _NOM[i] -= _DEN[i];
                _DEN[i] = 0;
            }
        }
    }

    _signPrefac *= rhs._signPrefac;
    resetAllCaches();
    return *this;
}

std::ostream& TFracNum::Print(std::ostream& out) const
{
    if (_debug) {
        out << "nom prime list: " << _NOM.size() << ",pointer " << to_string(_NOM) << endl;
        out << "den prime list: " << _DEN.size() << ",pointer " << to_string(_DEN) << endl;
        out << "NOM:";
        for (size_t i = 0; i < _NOM.size(); i++) {
            out << _NOM[i] << ",";
        }
        out << endl;
        out << "DEN:";
        for (size_t i = 0; i < _DEN.size(); i++) {
            out << _DEN[i] << ",";
        }
        out << endl;
    }
    out << "sign_prefac=" << _signPrefac << endl;

    bool integrity = true;
    for (size_t i = 0; i < _NOM.size(); i++) {
        if (_NOM[i] < 0 or _NOM[i] > 1000) {
            integrity = false;
        }
    }
    for (size_t i = 0; i < _DEN.size(); i++) {
        if (_DEN[i] < 0 or _DEN[i] > 1000) {
            integrity = false;
        }
    }
    if (not integrity) {
        return out;
    }

    long nom = 1;
    size_t ipn = 0;
    if (_signPrefac < 0) {
        out << "-NOM = ";
    } else {
        out << " NOM = ";
    }
    long FirstTerm = 1;
    while (ipn < _NOM.size()) {
        if (_NOM[ipn] != 0) {
            out << rpwa::primeNumbers::instance().primeNumber(ipn) << "^" << _NOM[ipn];
            FirstTerm = 0;
        }
        for (long jj = 0; jj < _NOM[ipn]; jj++) {
            nom *= rpwa::primeNumbers::instance().primeNumber(ipn);
        }
        ipn++;
        if (!FirstTerm and ipn < _NOM.size() and _NOM[ipn] != 0) {
            out << " * ";
        }
    }
    out << " = " << nom << endl;

    long den = 1;
    size_t ipd = 0;
    out << " DEN = ";
    FirstTerm = 1;
    while (ipd < _DEN.size()) {
        if (_DEN[ipd] != 0) {
            out << rpwa::primeNumbers::instance().primeNumber(ipd) << "^" << _DEN[ipd];
            FirstTerm = 0;
        }
        for (long jj = 0; jj < _DEN[ipd]; jj++)
            den *= rpwa::primeNumbers::instance().primeNumber(ipd);
        ipd++;
        if (!FirstTerm and ipd < _DEN.size() and _DEN[ipd] != 0) {
            out << " * ";
        }
    }
    out << " = " << den << endl;
    out << "NOM_INT=" << GetNumerator() << endl;
    out << "DEN_INT=" << GetDenominator() << endl;
    out << "dvalue=" << Dval() << endl;
    return out;
}


//TODO: fix this horrible sprintf mess
const char*
TFracNum::FracString() const
{
    char* formstr = new char[50];
    char* fstr = new char[100];
    const long& numerator = GetNumerator();
    const long& denominator = GetDenominator();
    if (numerator == 0) {
        sprintf(fstr, "0");
    } else if (denominator == 1) {
        sprintf(formstr, "%%c%s", IOUTSTRING);
        sprintf(fstr, formstr, _signPrefac < 0 ? '-' : '+', numerator);
    } else {
        sprintf(formstr, "%%c%s/%s", IOUTSTRING, IOUTSTRING);
        sprintf(fstr, formstr, _signPrefac < 0 ? '-' : '+', numerator, denominator);
    }
    return fstr;
}


const char NULLSTRING[1] = ""; // workaround CINT warning when sprintf(s,"");


//TODO: fix this horrible sprintf mess
const char*
TFracNum::FracStringSqrt() const
{
    char* formstr = new char[50];
    char* fstr = new char[200];
    if (GetNumerator() == 0) {
        sprintf(fstr, "0");
        return fstr;
    }
    size_t ipn = 0;
    rpwa::primeNumbers::entryType SQRT_NOM_INT = 1;
    rpwa::primeNumbers::entryType NOM_INT_REST = 1;
    while (ipn < _NOM.size()) {
        for (long jj = 0; jj < _NOM[ipn] / 2; jj++) {
            SQRT_NOM_INT *= rpwa::primeNumbers::instance().primeNumber(ipn);
        }
        if (_NOM[ipn] % 2) {
            NOM_INT_REST *= rpwa::primeNumbers::instance().primeNumber(ipn);
        }
        ipn++;
    }
    size_t ipd = 0;
    long SQRT_DEN_INT = 1;
    long DEN_INT_REST = 1;
    while (ipd < _DEN.size()) {
        for (long jj = 0; jj < _DEN[ipd] / 2; jj++) {
            SQRT_DEN_INT *= rpwa::primeNumbers::instance().primeNumber(ipd);
        }
        if (_DEN[ipd] % 2) {
            DEN_INT_REST *= rpwa::primeNumbers::instance().primeNumber(ipd);
        }
        ipd++;
    }

    char* sqrtstr = new char[100];
    bool one1 = false;
    bool one2 = false;
    if (SQRT_DEN_INT == 1) {
        if (SQRT_NOM_INT == 1) {
            sprintf(sqrtstr, "%s", NULLSTRING);
            one1 = true;
        } else {
            sprintf(formstr, "%s", IOUTSTRING);
            sprintf(sqrtstr, formstr, SQRT_NOM_INT);
        }
    } else {
        sprintf(formstr, "%s/%s", IOUTSTRING, IOUTSTRING);
        sprintf(sqrtstr, formstr, SQRT_NOM_INT, SQRT_DEN_INT);
    }

    char* reststr = new char[100];
    if (DEN_INT_REST == 1) {
        if (NOM_INT_REST == 1) {
            sprintf(reststr, "%s", NULLSTRING);
            one2 = true;
        } else {
            sprintf(formstr, "%s%s", SQUAREROOT_CHAR, IOUTSTRING);
            sprintf(reststr, formstr, NOM_INT_REST);
        }
    } else {
        sprintf(formstr, "%s%s/%s", SQUAREROOT_CHAR, IOUTSTRING, IOUTSTRING);
        sprintf(reststr, formstr, NOM_INT_REST, DEN_INT_REST);
    }

    if (one1 && one2) {
        sprintf(sqrtstr, "1");
    }
    sprintf(fstr, "%c%s%s", _signPrefac < 0 ? '-' : '+', sqrtstr, reststr);
    return fstr;
}


void TFracNum::removeZerosFromVector(vector<long>& vector)
{
    if (vector.empty()) {
        return;
    }
    size_t neededSize = vector.size();
    for ( ; neededSize > 0; --neededSize) {
        if (vector[neededSize - 1] != 0) {
            break;
        }
    }
    vector.resize(neededSize);
}


long TFracNum::getNumberFromFactorization(const vector<long>& vector)
{
    long retval = 1;
    for (size_t i = 0; i < vector.size(); ++i) {
        for (long j = 0; j < vector[i]; ++j) {
            retval *= rpwa::primeNumbers::instance().primeNumber(i);
        }
    }
    return retval;
}


// TODO: check if this can be deleted
#if(0)
TFracNum TFracNum::a_to_J(long J, long m)
{
    long kappa = (J - m) % 2;
    cout << "kappa=" << kappa << endl;
    long nom_ptr[1] = { 1 };
    TFracNum twofac(kappa, 0, nom_ptr, 0, 1);
    TFracNum fac1(J + m, 2 * J, "factorial");
    TFracNum fac2(J - m, 1, "factorial");
    return twofac * fac1 * fac2;
}
#endif

TFracNum TFracNum::am0_to_J(const long& J, const long& m, const long& m0)
{
    TFracNum twofac(vector<long>(1, m0),
                    vector<long>(),
                    1);
    TFracNum fac1(J + m, 2 * J, "factorial");
    TFracNum fac2(J - m, 1, "factorial");
    return twofac * fac1 * fac2;
}


TFracNum TFracNum::c_sub_ell(const long& ell)
{
    if (ell == 0) {
        return TFracNum::One;
    }
    TFracNum two_to_ell(vector<long>(1, ell),
                        vector<long>(),
                        1);
    TFracNum fac1(ell, 1, "factorial");
    TFracNum fac2(ell, 2 * ell, "factorial");
    return two_to_ell * fac1 * fac2;
}


TFracNum TFracNum::cm0_sub_ell(const long& ell, const long& m0)
{
    if (ell == 0) {
        return TFracNum::One;
    }
    TFracNum two_to_ell(vector<long>(1, (ell + m0) / 2),
                        vector<long>(),
                        1);
    TFracNum fac1(ell, 1, "factorial");
    TFracNum fac2(ell, 2 * ell, "factorial");
    return two_to_ell * fac1 * fac2;
}


TFracNum TFracNum::cm0_sub_ell_2(const long& ell, const long& m0)
{
    //return  am0_to_J(ell, 0, m0);
    if (ell == 0) {
        return TFracNum::One;
    }
    TFracNum two_to_ell(vector<long>(1, (ell + m0) ),
                        vector<long>(),
                        1);
    TFracNum fac1a(ell, 1, "factorial");
    TFracNum fac1b(ell, 1, "factorial");
    TFracNum fac2a(ell, 2 * ell, "factorial");
    TFracNum fac2b(ell, 2 * ell, "factorial");
    return two_to_ell * fac1a * fac1b * fac2a * fac2b;
}
