#include "JSS.h"

#include "JLSIndices.h"
#include "MathUtils.h"

#include <algorithm>
#include <functional>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>

namespace relhel {

//-------------------------
const std::string to_string(unsigned omega, unsigned epsilon, unsigned phi)
{
    return (phi == 0) ? std::string(omega, 'o') + std::string(epsilon, 'e')
        : std::string(phi, '#') + "(" + to_string(omega, epsilon, 0) + ")";
}

//-------------------------
JSS::JSS(const QuantumNumbers& parent, const std::vector<QuantumNumbers>& daughters) :
    Parent_(parent),
    Daughters_(daughters)
{
    if (Daughters_.size() < 2)
        throw std::invalid_argument("At least two daughters must be specified.");

    if (Daughters_.size() > 2)
        throw std::invalid_argument("Only two daughters may be specified.");

    if (is_odd(Parent_.twoJ()) or std::any_of(Daughters_.begin(), Daughters_.end(), is_fermion))
        throw std::invalid_argument("Only valid for integer spin");
}

//-------------------------
void JSS::CalcAmpl()
{
    // get S allowed by daughters' J
    auto two_S = triangle(Daughters_[0].twoJ(), Daughters_[1].twoJ());
    
    // get intrinsic parity of parent * daughters
    auto intr_parity = Parent_.parity() * parity(Daughters_);

    // get set of allowed total orbital angular momentum
    std::set<unsigned> L;
    for (auto two_s : two_S) 
        for (auto two_l : triangle(Parent_.twoJ(), two_s))
            if ((two_l / 2) % 2 == (intr_parity == Parity::negative))
                L.insert(two_l / 2);

    auto daughter_two_j_sum = spin_sum(Daughters_);
    
    // whether an even number of contractions are needed
    const bool even_contraction = is_even((Parent_.twoJ() + daughter_two_j_sum) / 2 + *L.begin());

    // Phi = wave function of parent state
    // Psi = wave function of two-particle final state psi(j1, j2; S, m)
    // Chi = wave function for orbital angular momentum chi(l, 0)
    
    // loop over allowed orbital angular momenta
    for (auto l : L) {
        
        // get intrinsic spins allowed by J and L
        auto two_S_L = triangle(Parent_.twoJ(), 2 * l);
        // remove those not allowed by the daughter spins
        two_S_L.erase(std::remove_if(two_S_L.begin(), two_S_L.end(), [&](unsigned two_s){return two_s < two_S.front() or two_s > two_S.back();}), two_S_L.end());

        // loop over spin projections
        for (auto two_s : two_S_L) {

            // create a counter object
            // \todo: move up to outside the highest loop?
            JLSIndices K(Parent_.twoJ(), l, two_s);
            
            // unsigned max_contractionNumber = 0;

            for (K.psiInternal() = 0; K.psiInternal() <= std::min(Daughters_[0].twoJ(), Daughters_[1].twoJ()) / 2; ++K.psiInternal()) {
                for (K.chiPhi() = 0; K.chiPhi() <= std::min(Parent_.twoJ() / 2, l); ++K.chiPhi()) {
                    for (K.psiChi() = 0; K.psiChi() <= std::min(daughter_two_j_sum / 2, l); ++K.psiChi()) {
                        for (K.psiPhi() = 0; K.psiPhi() <= std::min(daughter_two_j_sum, Parent_.twoJ()) / 2; ++K.psiPhi()) {

                            /////////////////////////
                            // check selection rules
                            if (K.psiInternal() + K.psiChi() + K.chiPhi() + K.psiPhi() != ((K.twoJ() + daughter_two_j_sum) / 2 + K.L() - (even_contraction ? 0 : 3)) / 2)
                                continue;
                            
                            if (even_contraction) {
                                if (2 * K.psiInternal() + K.psiChi() + K.psiPhi() != daughter_two_j_sum / 2 )
                                    continue;
                                
                                if (K.psiChi() + K.chiPhi() != K.L())
                                    continue;
                                
                                if (K.chiPhi() + K.psiPhi() != K.twoJ() / 2)
                                    continue;
                            } else {
                                if (K.L() - K.psiChi() - K.chiPhi() > 1)
                                    continue;
                                
                                if (K.twoJ() / 2 - K.psiPhi() - K.chiPhi() > 1)
                                    continue;
                            }
                            /////////////////////////

                            /////////////////////////
                            // Get contraction ranks
                            std::vector<unsigned> ranks;
                            ranks.reserve(Daughters_.size());
                            std::transform(Daughters_.begin(), Daughters_.end(), std::back_inserter(ranks), [&K](const QuantumNumbers& q){return q.twoJ() / 2 - K.psiInternal();});
                            if (!even_contraction) {
                                switch (daughter_two_j_sum / 2 - 2 * K.psiInternal() - K.psiChi() - K.psiPhi()) {
                                case 0:
                                case 1:
                                    break;
                                case 2:
                                    if (std::any_of(ranks.begin(), ranks.end(), [](unsigned r){return r == 0;}))
                                        break;
                                    std::for_each(ranks.begin(), ranks.end(), [](unsigned& r){--r;});
                                    break;
                                default:
                                    throw std::runtime_error("PsiRest < 0 or > 2");
                                }
                            }
                            /////////////////////////
                            
                            for (K.chiOmega() = 0; K.chiOmega() <= ranks[0] and K.chiOmega() <= K.psiChi(); ++K.chiOmega()) {
                                for (K.phiOmega() = 0; K.phiOmega() <= ranks[0] - K.chiOmega() and K.phiOmega() <= K.psiPhi(); ++K.phiOmega()) {

                                    K.chiEps() = K.psiChi() - K.chiOmega();
                                    K.phiEps() = K.psiPhi() - K.phiOmega();

                                    if (K.phiEps() + K.chiEps() > ranks[1])
                                        continue;

                                    for (K.delta() = 0; K.delta() <= std::min(two_s, K.twoJ()) / 2; ++K.delta()) {

                                        std::cout << spin_to_string(K.twoJ()) << " "
                                                  << std::to_string(K.L()) << " "
                                                  << spin_to_string(K.twoS()) << " "
                                                  << std::to_string(K.psiInternal()) << " "
                                                  << std::to_string(K.chiPhi()) << " "
                                                  << std::to_string(K.psiChi()) << " "
                                                  << std::to_string(K.psiPhi()) << " "
                                                  << std::to_string(K.chiOmega()) << " "
                                                  << std::to_string(K.phiOmega()) << " "
                                                  << std::to_string(K.phiEps()) << " "
                                                  << std::to_string(K.chiEps()) << " "
                                                  << std::to_string(K.delta()) << " "
                                                  << std::endl;
                                        // auto it = std::find(LSAmplitudes_.begin(), LSAmplitudes_.end(), K);
                                        // K.contractionNumber() = (it != LSAmplitudes_.end()) ? it->contractionNumber() : MaxContractionNumber + 1;

                                        // LSAmplitudes_.emplace_back(K, Daughters_);
                                        // if (LSAmplitudes_.back().tensorSum().GetNterms() > 0)
                                        //     MaxContractionNumber = std::max(LSAmplitudes_.back().contractionNumber(), MaxContractionNumber);
                                        // else
                                        //     LSAmplitudes_.pop_back();
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // if (!FhhAmpl_.empty())
    //     throw std::runtime_error("FhhAmpl_ not empty.");

    // for (unsigned lambda = 0; lambda <= Daughters_[0].twoJ(); ++lambda)
    //     for (unsigned nu = (lambda == 0) ? 0 : -Daughters_[1].twoJ(); nu <= Daughters_[1].twoJ(); ++nu)
    //         FhhAmpl_.emplace_back(Parent_.twoJ(), Daughters_.twoJ(), lambda, nu, LSAmplitudes_, even_contraction);

    // // remove amps with 0 terms
    // FhhAmpl_.erase(std::remove_if(FhhAmpl_.begin(), FhhAmpl_.end(), [](const TFhh& f){return f.LSt().empty();}), FhhAmpl_.end());
    
    // WHAT IS THE POINT OF THIS?
    // if identical particles:
    // if (identical_particle) {
    //     FhhIdAmpl_.reserve(FhhAmpl_.size());
    //     for (const auto& f : FhhAmpl_)
    //         if (is_nu_nu(f))
    //             FhhIdAmpl_.emplace_back(f, TFhh::symmetry::nu_nu);
    //         else if (is_nu_minus_nu(f))
    //             FhhIdAmpl_.emplace_back(f, TFhh::symmetry::nu_minus_nu);
    //         else {
    //             auto it = std::find_if(FhhAmpl_.begin(), FhhAmpl_.end(), std::bind(nu_lambda_partners, f, std::placeholders::_1));
    //             if (it == FhhAmpl_.end())
    //                 throw std::runtime_error("Na partner for amplitude " + f.name());
    //             FhhIdAmpl_.emplace_back(f, *it);
    //         }
    // }
    
    // cout << FhhAmpl_.size() << " amplitudes: non-relativistic limit" << endl;
    // for (auto& f : FhhAmpl_)
    //     f.NonRelLimit();

    // cout << "Check non-relativistic G's" << endl;
    // for (const auto& f : FhhAmpl_)
    //     f.PrintNRG();
}

}

