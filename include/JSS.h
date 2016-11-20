#ifndef relhel__JSS_h
#define relhel__JSS_h

/* #include "TFhh.h" */
/* #include "TLSAmpl.h" */

#include "QuantumNumbers.h"

#include <vector>

namespace relhel {

/// 
class JSS
{
public:

    /// contsructor
    JSS(const QuantumNumbers& parent, const std::vector<QuantumNumbers>& daughters);

    /// contsructor
    JSS(const QuantumNumbers& parent, const QuantumNumbers& d1, const QuantumNumbers& d2)
        : JSS(parent, {d1, d2}) {}

    /* std::vector<TFhh>& fhh() */
    /* { return FhhAmpl_; } */

    void CalcAmpl();

private:

    /// Parent QuantumNumbers
    QuantumNumbers Parent_;

    /// Daughter QuantumNumbers
    std::vector<QuantumNumbers> Daughters_;

    /* std::vector<LSAmplitude> LSAmplitudes_; */
    /* std::vector<TFhh>    FhhAmpl_; */
    /* std::vector<TFhh>    FhhIdAmpl_; */

};

}

#endif
