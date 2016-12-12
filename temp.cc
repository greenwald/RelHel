std::vector<unsigned> b_vec(size_t s, size_t z)
{
    std::vector<unsigned> v(s, 1);
    v[z] = 0;
    return v;
}

void temp()
{
    std::vector<unsigned> r = {1, 1, 1, 0, 1};
    
    std::vector<std::vector<unsigned> > B;

    for (unsigned i = 0; i <= r.size() - 3; ++i)
        if (r[i] > 0)
            for (unsigned j = i + 1; j <= r.size() - 2; ++j)
                if (r[j] > 0)
                    for (unsigned k = j + 1; k <= r.size() - 1; ++k)
                        if (r[k] > 0)
                            std::cout << i << j << k << std::endl;
    
}
//     std::vector<unsigned> r(4, 0);
//     for (r[3] = 0; r[3] <= 2; ++r[3])
//         for (r[0] = 0; r[0] <= 2; ++r[0])
//             for(r[1] = r[0]; r[1] <= 2; ++r[1])
//                 for(r[2] = 0; r[2] <= 2; ++r[2]) {

//                     if (r != std::vector<unsigned>({2, 2, 1, 2}))
//                         continue;
                    
//                     if (abs(static_cast<int>(r[3]) - static_cast<int>(r[2])) > r[0] + r[1] or
//                         abs(static_cast<int>(r[0]) - static_cast<int>(r[1])) > r[2] + r[3])
//                         continue;

//                     // possible determine epsilon contraction
//                     std::vector<std::vector<unsigned> > B;
                    
//                     // if even rank
//                     if (std::accumulate(r.begin(), r.end(), 0u) % 2 == 0)
//                         B.push_back(std::vector<unsigned>(r.size(), 0));
//                     else {
//                         for (size_t i = 0; i < r.size(); ++i) {
//                             std::vector<unsigned> b(r.size(), 1);
//                             b[i] = 0;
//                             if (b < r)
//                                 B.push_back(b);
//                         }
//                     }
                    
//                     std::cout << r[3] << " -> " << r[0] << " + " << r[1] << " l = " << r[2] << std::endl;
//                     for (const auto& b : B) {
                    
//                         // store available ranks
//                         std::vector<unsigned> rb;
//                         rb.reserve(r.size());
//                         std::transform(r.begin(), r.end(), b.begin(), std::back_inserter(rb), [](unsigned rr, unsigned bb){return rr - bb;});

//                         std::vector<std::vector<unsigned> > C(r.size() - 1, std::vector<unsigned>(r.size(), 0));

//                         while (C.back().back() <= std::min(rb[rb.size() - 2], rb.back())) {

//                             std::vector<unsigned> c(rb.size(), 0);
//                             for (size_t i = 0; i < C.size(); ++i)
//                                 for (size_t j = i + 1; j < C[i].size(); ++j) {
//                                     c[i] += C[i][j];
//                                     c[j] += C[i][j];
//                                 }
//                             if (c == rb) {
//                                 for (size_t i = 0; i < C.size(); ++i)
//                                     for (size_t j = i + 1; j < C[i].size(); ++j)
//                                         std::cout << C[i][j] << " * (" << i << "," << j << "), ";
//                                 if (std::accumulate(b.begin(), b.end(), 0u) == 0)
//                                     std::cout << "\b\b";
//                                 else {
//                                     std::cout << "1 * (";
//                                     for (size_t i = 0; i < b.size(); ++i)
//                                         if (b[i] > 0)
//                                             std::cout << i << ", ";
//                                     std::cout << "\b\b)";
//                                 }
//                                 std::cout << std::endl;
//                             }

//                             ++C[0][1];
                            
//                             size_t i = 0;
//                             size_t j = 1;
//                             while(C[i][j] > std::min(rb[i], rb[j])) {
//                                 C[i][j] = 0;
                        
//                                 ++j;
//                                 if (j >= C[i].size()) {
//                                     ++i;
//                                     j = i + 1;
//                                 }

//                                 ++C[i][j];

//                                 if (i == C.size() - 1 and j == C.size())
//                                     break;
//                             }
//                         }
                    
//                         // for (const auto& b : B) {
//                         //     for (unsigned bi : b)
//                         //         std::cout << bi;
//                         //     std::cout << std::endl;
//                         // }
//                     }
//                 }
// }                    
/*                    
                    for (int b1 = 0; b1 <= r1 and b1 <= B; ++b1)
                        for (int b2 = 0; b2 <= r2 and b2 <= B; ++b2)
                            for (int b0 = 0; b0 <= r0 and b0 <= B; ++b0)
                                for (int bl = 0; bl <= l and bl <= B; ++bl) {
                                    
                                    if ((b1 + b2 + b0 + bl) != 3 * B)
                                        continue;
                                    
                                    for (int psi1_psi2 = 0; psi1_psi2 <= std::min(r1 - b1, r2 - b2); ++psi1_psi2)
                                        
                                        for (int psi1_chi = 0; psi1_chi <= std::min(r1 - b1 - psi1_psi2, l - bl); ++psi1_chi)
                                            for (int psi2_chi = 0; psi2_chi <= std::min(r2 - b2 - psi1_psi2, l - bl - psi1_chi); ++psi2_chi)
                                                
                                                for (int psi1_phi = 0; psi1_phi <= std::min(r1 - b1 - psi1_psi2 - psi1_chi, r0 - b0); ++psi1_phi)
                                                    for (int psi2_phi = 0; psi2_phi <= std::min(r2 - b2 - psi1_psi2 - psi2_chi, r0 - b0 - psi1_phi); ++psi2_phi)
                                                        
                                                        for (int chi_phi = 0; chi_phi <= std::min(l - bl - psi1_chi - psi2_chi, r0 - b0 - psi1_phi - psi2_phi); ++chi_phi) {
                                                            
                                                            if (2 * (psi1_psi2 + psi1_chi + psi2_chi + psi1_phi + psi2_phi + chi_phi) + 3 * B != r0 + r1 + r2 + l)
                                                                continue;
                                                            
                                                            std::cout << "\t"
                                                                      << " psi1_psi2 = " << psi1_psi2
                                                                      << " psi1_chi = " << psi1_chi
                                                                      << " psi2_chi = " << psi2_chi
                                                                      << " psi1_phi = " << psi1_phi
                                                                      << " psi2_phi = " << psi2_phi
                                                                      << " chi_phi = " << chi_phi
                                                                      << " B = " << b1 << b2 << bl << b0
                                                                      << std::endl;
                                                        }
                                    
                                }
                }
            }
}
*/
