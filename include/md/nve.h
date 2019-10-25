#ifndef MD_NVE_H_
#define MD_NVE_H_

// Liscio headers
#include "prime/base.h"

namespace uovie {
namespace md {

    /*** ================================================================== ***/
    /*** Molecular Dynamics (NVE) [customized]                              ***/
    /*** ================================================================== ***/

    class md_nve_ver {
    public:
        md_nve_ver() = default;

        void exec();

    private:
        double rt = 100;    // run time
        double ss = 0.0001; // time step size
        double dcp = 1000;  // data collection period

        double omega = 1;   // angular frequency

        double m = 1;       // mass
        double q = 1;       // position
        double p = 0;       // momentum
        double F = -m * pow(omega, 2.0) * q; // force

        double Etot = 0;    // total energy
    };

} // !md
} // !uovie
#endif // !MD_NVE_H_

