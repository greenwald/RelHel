#include "LSAmpl.h"

#include "SpinWaveFunction.h"

#include <stdexcept>

namespace relhel {

//-------------------------
LSAmpl::LSAmpl(const JLSIndices& jlsi, const std::vector<QuantumNumbers>& daughters)
    : JLSIndices(jlsi)
{
    
    auto total_rank = (spin_sum(daughters) + twoJ()) / 2 + L();
    auto even_contraction = is_even(total_rank);

    if (2 * (psiInternal() + psiChi() + chiPhi() + psiPhi()) != (total_rank - (even_contraction ? 0 : 3)))
        throw std::runtime_error("Invalid contraction occurred.");

    TSpinWaveFunction WFS1(daughters[0], 's');
    TSpinWaveFunction WFS2(daughters[1], 's');
    TTensorSum TSS = WFS1.GetSpinCoupledTensorSum(WFS2, delta(), S());

    CoupledWaveFunctions psi(daughters[0].twoJ(), daughters[1].twoJ());
    
    if (psiInternal() > 0 and TSS.SpinInnerContraction(psiInternal()) == 0)
        return;

    TSpinWaveFunction WFL(L(), 'l');
    TTensorSum TSL = WFL.GetTensorSum('c', 0);

    TTensorSum TSLS = TSS.LSContraction(TSL, psiChi(), chiOmega(), chiEps(), 'c');
    if (TSLS.GetNterms() == 0)
        return;

    TSpinWaveFunction WFJ(J(), 'c');
    TTensorSum TSJ = WFJ.GetTensorSum('p', delta());

    TTensorSum TSLSJ = TSLS.LSContraction(TSJ, psiPhi(), phiOmega(), phiEps(), 'p');

    if (TSLSJ.GetNterms() == 0)
        return;

    TSScalar_ = TSLSJ.LJContraction(chiPhi(), even_contraction);
}

}
