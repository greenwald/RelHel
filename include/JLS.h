#ifndef __JLS__h
#define __JLS__h

#include <string>

/// \class JLS
/// \brief Holds J, L, S for objects that need them
/// \author Daniel Greenwald
class JLS
{
public:
    /// constructor
    JLS(unsigned j, unsigned l, unsigned s)
        : J_(j), L_(l), S_(s) {}

    const unsigned J() const { return J_; }
    unsigned& J() { return J_; }

    const unsigned L() const { return L_; }
    unsigned& L() { return L_; }

    const unsigned S() const { return S_; }
    unsigned S() { return S_; }
    
private:

    unsigned J_;
    unsigned L_;
    unsigned S_;

};

/// convert to string
inline std::string to_string(const JLS& jls)
{ return "(" + std::to_string(jls.J()) + ";" + std::to_string(jls.L()) + "," + std::to_string(jls.S()) + ")"; }

/// equality operator
const bool operator==(const JLS& lhs, const JLS& rhs)
{ return lhs.J() == rhs.J() and lhs.L() == rhs.L() and lhs.S() == rhs.S(); }

/// inequality operator
const bool operator!=(const JLS& lhs, const JLS& rhs)
{ return !(lhs == rhs); }

#endif
