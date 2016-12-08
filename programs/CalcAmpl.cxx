// #include "TFhh.h"
// #include "TJSS.h"

#include "Contraction.h"
#include "WaveFunction.h"

#include <iostream>
#include <stdexcept>
#include <vector>

using namespace relhel;

//-------------------------
// printing function
std::string to_string(unsigned count, const std::vector<Contraction::Term::index>& I)
// { if (count == 0) return ""; return ((count > 1) ? (std::to_string(count) + " x ") : "") + "(" + to_string(I) + "); "; }
{ return std::to_string(count) + " x (" + to_string(I) + "); "; }

int main(int narg, char** args)
{

    if (narg < 6) {
        std::cout << "Usage: " << args[0] << " (Mother 2J) (Daughter 2J) (Daughter 2J) (L) (2S) [2M]" << std::endl;
        return 1;
    }

    WaveFunction phi(std::atoi(args[1]));
    WaveFunction psi1(std::atoi(args[2]));
    WaveFunction psi2(std::atoi(args[3]));

    if (spin(phi) != 0 and narg < 6)
        throw std::invalid_argument("M is ambiguous, you must state it explicitly; main");
    
    unsigned L = std::atoi(args[4]);
    unsigned two_S = std::atoi(args[5]);
    unsigned two_M = std::atoi(args[6]);

    std::cout << "(" + spin_to_string(spin(phi)) + ") -> (" + spin_to_string(spin(psi1)) + ") + (" + spin_to_string(spin(psi2)) + ")"
              << ", L = " << L
              << ", S = " << spin_to_string(two_S)
              << ", M = " << spin_to_string(two_M)
              << std::endl;

    if (!triangle(spin(psi1), spin(psi2), two_S))
        throw std::invalid_argument("j1:j2:S triangle broken; main");

    if (!triangle(spin(phi), 2 * L, two_S))
        throw std::invalid_argument("j0:L:S triangle broken; main");
    
    CoupledWaveFunctions psi(psi1, psi2, two_S, two_M);

    OrbitalAngularMomentumWaveFunction chi(L);

    Contraction C(psi, chi, phi);

    GammaPolynomial P;
    
    for (unsigned psi_internal = 0; psi_internal <= std::min(rank(psi1), rank(psi2)); ++psi_internal) {

        auto K1 = C;
        
        K1.contract(contractions::psi1_psi2, psi_internal);
        if (K1.terms().empty()) continue;
        
        for (unsigned psi_chi = 0; psi_chi <= std::min(rank(psi), rank(chi)); ++psi_chi) {
            for (unsigned psi1_chi = 0; psi1_chi <= (rank(psi1) - psi_internal) and psi1_chi <= psi_chi; ++psi1_chi) {
                auto K2 = K1;
                
                K2.contract(contractions::psi1_chi, psi1_chi);
                if (K2.terms().empty()) continue;
                
                // limit psi2_chi to remaining rank of psi2 after contraction with psi1
                unsigned psi2_chi = std::min(psi_chi - psi1_chi, rank(psi2) - psi_internal);
                K2.contract(contractions::psi2_chi, psi2_chi);
                if (K2.terms().empty()) continue;

                for (unsigned psi_phi = 0; psi_phi <= std::min(rank(psi), rank(phi)); ++psi_phi) {
                    for (unsigned psi1_phi = 0; psi1_phi <= (rank(psi1) - psi_internal - psi1_chi) and psi1_phi <= psi_phi; ++psi1_phi) {
                        auto K3 = K2;
                        
                        K3.contract(contractions::psi1_phi, psi1_phi);
                        if (K3.terms().empty()) continue;

                        // limit psi2_phi to remaining rank of psi2 after contraction with psi1 and chi
                        unsigned psi2_phi = std::min(psi_phi - psi1_phi, rank(psi2) - psi_internal - psi2_chi);
                        K3.contract(contractions::psi2_phi, psi2_phi);
                        if (K3.terms().empty()) continue;

                        for (unsigned chi_phi = 0; chi_phi <= std::min(rank(phi) - psi1_phi - psi2_phi, rank(chi) - psi1_chi - psi2_chi); ++chi_phi) {
                            auto K4 = K3;

                            K4.contract(contractions::chi_phi, chi_phi);
                            if (K4.terms().empty()) continue;

                            auto I3 = triple_contraction(K4);
                            if (I3.size() == 3) {
                                K4.contract(I3);
                                if (K4.terms().empty()) continue;
                            }

                            if (rank(K4) == 0) {

                                // std::cout << "contractions: "
                                //           << to_string(psi_internal, contractions::psi1_psi2)
                                //           << to_string(psi1_chi,     contractions::psi1_chi)
                                //           << to_string(psi2_chi,     contractions::psi2_chi)
                                //           << to_string(psi1_phi,     contractions::psi1_phi)
                                //           << to_string(psi2_phi,     contractions::psi2_phi)
                                //           << to_string(chi_phi,      contractions::chi_phi)
                                //           << to_string(I3.size() == 3, I3)
                                //           << "produce " << K4.terms().size() << " terms" << std::endl;

                                P += gamma_polynomial(K4.terms());
                            }
                        }
                    }
                }
            }
        }
    }

    std::cout << to_string(P) << std::endl;
    
    return 0;
}

