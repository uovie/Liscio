/* MD_NVE_CPP_ */

// Liscio headers
#include "md/nve.h"

namespace uovie {
namespace md {

    /*** ================================================================== ***/
    /*** Molecular Dynamics (NVE) [customized]                              ***/
    /*** ================================================================== ***/

    void md_nve_ver::exec() {
        std::ofstream otime("sho.time", std::ios_base::binary);
        std::ofstream opos("sho.pos", std::ios_base::binary);
        std::ofstream omom("sho.mom", std::ios_base::binary);
        std::ofstream oetot("sho.etot", std::ios_base::binary);

        int nstep = 0;
        for (double t = 0; t < rt; t += ss) {
            p += F * ss / 2;
            q += p * ss / m;
            F = -m * pow(omega, 2.0) * q;
            p += F * ss / 2;

            nstep++;
            if (nstep == dcp) {
                otime.write(prime::as_bytes(t), sizeof(double));
                opos.write(prime::as_bytes(q), sizeof(double));
                omom.write(prime::as_bytes(p), sizeof(double));
                Etot = pow(p, 2) / (2 * m) + 0.5 * m * pow(omega * q, 2);
                oetot.write(prime::as_bytes(Etot), sizeof(double));
                nstep = 0;
            }
        }

        otime.close();
        opos.close();
        omom.close();
        oetot.close();
    }

} // !md
} // !uovie