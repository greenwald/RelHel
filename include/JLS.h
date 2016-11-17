#ifndef __JLS__h
#define __JLS__h

#include <MathUtils.h>

#include <string>

namespace relhel {

/// \class JLS
/// \brief Holds J, L, S for objects that need them
/// \author Daniel Greenwald
class JLS
{
public:
    /// constructor
    /// \param two_j (twice) the parent spin
    /// \param l orbital angular momentum
    /// \param two_s (twice) the spin angular momentum
    JLS(unsigned two_j, unsigned l, unsigned two_s)
        : TwoJ_(two_j), L_(l), TwoS_(two_s) {}

    const unsigned twoJ() const
    { return TwoJ_; }

    unsigned& twoJ()
    { return TwoJ_; }

    const unsigned L() const
    { return L_; }
    
    unsigned& L()
    { return L_; }

    const unsigned twoS() const
    { return TwoS_; }
    
    unsigned twoS()
    { return TwoS_; }
    
private:

    unsigned TwoJ_;
    unsigned L_;
    unsigned TwoS_;

};

/// convert to string
inline std::string to_string(const JLS& jls)
{ return "(" + spin_to_string(jls.twoJ()) + ";" + std::to_string(jls.L()) + "," + spin_to_string(jls.twoS()) + ")"; }

/// equality operator
const bool operator==(const JLS& lhs, const JLS& rhs)
{ return lhs.twoJ() == rhs.twoJ() and lhs.L() == rhs.L() and lhs.twoS() == rhs.twoS(); }

/// inequality operator
const bool operator!=(const JLS& lhs, const JLS& rhs)
{ return !(lhs == rhs); }

}

#endif
