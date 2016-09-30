#ifndef TJWFTENSOR_HH
#define TJWFTENSOR_HH

#include <iostream>
#include <vector>

#include "RationalNumber.h"

class TensorTerm
{
public:

    /// \enum index_type
    /// index type for components
    enum class index : unsigned { omega, epsilon, chi, phi };

    /// constructor
    TensorTerm(index i, const std::vector<int>& pzm, const RationalNumber& prefac)
        : Prefactor_(prefac)
    { PZM_[i] = pzm; }

    /// copy-ish constructor
    TensorTerm(const TensorTerm& S, const TensorTerm& L,
               unsigned contractions, int o_contractions, index cp_contraction)

    /// \return Prefactor_
    const RationalNumber& prefactor() const
    { return Prefactor_; }

    /// \return GammaPot_
    const std::array<int, 2>& gammaPot() const
    { return GammaPot_; }
    
    /// access operator
    const std::vector<int>& operator[](index_type i) const
    { return PZM_[i]; }
    
    long LJContraction       (const long& ncon, const bool& even);
// TODO: optimize this call
    void Multiply            (const char& name, const std::vector<long>& pzm_field, const TFracNum& prefac);
    long SpinInnerContraction(const long& cPsiInt);
    bool SameStructure       (const TTensorTerm& rhs) const;
    void AddTwoTerms         (const TTensorTerm& rhs);
    bool IsNonZero           () const { return (_prefac != TFracNum::Zero); }

    std::ostream& Print(const char& flag = 'n', std::ostream& out = std::cout) const;

    const TFracNum& GetPreFac() const { return _prefac; }
    const long&     GetGamS()   const { return _gam_s_pot; }
    const long&     GetGamSig() const { return _gam_sig_pot; }

private:

    void shrinkVectors(const size_t& rOme,
                       const size_t& rEps,
                       const size_t& rChi,
                       const size_t& rPhi);


    std::array<std::vector<int>, 4> PZM_;
    std::array<int, 2> GammaPot_({0, 0});

    /* std::vector<long> _ome_pzm; */
    /* std::vector<long> _eps_pzm; */
    /* std::vector<long> _chi_pzm; */
    /* std::vector<long> _phi_pzm; */

    long _gam_s_pot{0};
    long _gam_sig_pot{0};

    RationalNumber Prefactor_{RationalNumber(1)};

};

using TTensorTerms = std::vector<TTensorTerm>;

inline
std::ostream&
operator <<(std::ostream&            out,
            const TTensorTerm&       tensorTerm)
{
    return tensorTerm.Print('n', out);
}


class TTensorSum
{


public:
    TTensorSum() : _terms() { }

    void AddTerm(const TTensorTerm& addt);
    size_t SpinInnerContraction(const long& cPsiInt);

    TTensorSum LSContraction(const TTensorSum& L,
                             const long& contr,
                             const long& co,
                             const long& ce,
                             const char& conType) const;

    TTensorSum LJContraction(const long& cChiPhi, const bool& even);

    void removeZeroTerms();
    size_t GetNterms() const { return _terms.size(); }

    void Print(char = 'n') const;

    const TTensorTerm& GetTerm(const size_t& i) const { return _terms[i]; }

    const TTensorTerms& terms() const
    { return _terms; }
    
private:

    TTensorTerms _terms;

};

#endif
