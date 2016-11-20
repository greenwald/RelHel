#ifndef relhel__LSAmpl_h
#define relhel__LSAmpl_h

#include "JLSIndices.h"
#include "QuantumNumbers.h"

#include "TJwfTensor.h"

#include <vector>

namespace relhel {

/// \brief Relativistic LS-coupling amplitudes
/// \author Jan.Friedrich@ph.tum.de, Daniel Greenwald
class LSAmpl : public JLSIndices
{
public:

    /// Constructor
    LSAmpl(const JLSIndices& jlsi, const std::vector<QuantumNumbers>& daughters);
    
    const TTensorSum& tensorSum() const
    { return TSScalar_; }

    /* const TTensorTerm& GetTerm(const long& i) const */
    /* { return _TSScalar.GetTerm(i); } */

private:

    TTensorSum TSScalar_;
};

}

#endif
