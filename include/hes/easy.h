#ifndef HES_EASY_H_
#define HES_EASY_H_

// Liscio headers
#include "prime/base.h"

namespace uovie {
namespace hes {

    constexpr auto pi = prime::pi;
    constexpr auto h_bar = prime::h_bar;

    /*** ================================================================== ***/
    /*** General Hamiltonian Eigen Solver                                   ***/
    /*** ================================================================== ***/

    class ghes {
    public:
        ghes() = default;
        ghes(const int& _nbasis, const int& _nstate) : nbasis(_nbasis), nstate(_nstate) { }

    protected:
        const int nbasis = 1;
        const int nstate = 1;

        Eigen::MatrixXd H; // Hamiltonian
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;

        virtual double mel(const int& i, const int& j) = 0; // matrix element
        virtual double bsf(const int& i, const double& x) = 0; // basis function

    public:
        void eigen_solve();
        void print_evals();
        void print_evecs(const std::pair<double, double>& xlim);
    };

    /*** ================================================================== ***/
    /*** GHES with ISW Basis Function                                       ***/
    /*** ================================================================== ***/

    class ghes_isw : public ghes {
    public:
        ghes_isw() = default;
        ghes_isw(const double& _m, const double& _a, const int& _nbasis, const int& _nstate)
            : m(_m), a(_a), ghes(_nbasis, _nstate) { }

    protected:
        const double m = 1;
        const double a = 1;  // isw: well width

        double bsf(const int& i, const double& x) override; // basis function
    };

    /*** ================================================================== ***/
    /*** GHES with SHO Basis Function                                       ***/
    /*** ================================================================== ***/

    class ghes_sho : public ghes {
    public:
        ghes_sho() = default;
        ghes_sho(const double& _m, const double& _omega_s, const int& _nbasis, const int& _nstate)
            : m(_m), omega_s(_omega_s), ghes(_nbasis, _nstate) { }

    protected:
        const double m = 1;
        const double omega_s = 1;  // sho: angular frequency

        double bsf(const int& i, const double& x) override; // basis function
    };

    /*** ================================================================== ***/
    /*** Simple Harmonic Oscillator on Infinite Square Well Basis           ***/
    /*** ================================================================== ***/

    class sho_on_isw : public ghes_isw {
    public:
        sho_on_isw() = default;
        sho_on_isw(const double& _m, const double& _omega_d, const double& _a,
            const int& _nbasis, const int& _nstate): omega_d(_omega_d),
            ghes_isw(_m, _a, _nbasis, _nstate) { }

    private:
        const double omega_d = 1;  // sho: angular frequency

        double mel(const int& i, const int& j) override;
    };

    /*** ================================================================== ***/
    /*** Infinite Square Well on Simple Harmonic Oscillator Basis (US)      ***/
    /*** ================================================================== ***/

    class isw_on_sho : public ghes_sho {
    public:
        isw_on_sho() = default;
        isw_on_sho(const double& _m, const double& _a, const double& _omega_s,
            const int& _nbasis, const int& _nstate) : a(_a),
            ghes_sho(_m, _omega_s, _nbasis, _nstate) { }

    private:
        const double a = 1; // isw: well width

        double mel(const int& i, const int& j) override;
    };

    /*** ================================================================== ***/
    /*** Double Parabola on Infinite Square Well Basis                      ***/
    /*** ================================================================== ***/

    class dpb_on_isw : public ghes_isw {
    public:
        dpb_on_isw() = default;
        dpb_on_isw(const double& _m, const double& _b, const double& _omega_d, const double& _a,
            const int& _nbasis, const int& _nstate) : b(_b), omega_d(_omega_d),
            ghes_isw(_m, _a, _nbasis, _nstate) { }

    private:
        const double b = 1; // dpb: minimum point
        const double omega_d = 1; // dpb: angular frequency

        double mel(const int& i, const int& j) override;
    };

    /*** ================================================================== ***/
    /*** Double Parabola on Simple Harmonic Oscillator Basis                ***/
    /*** ================================================================== ***/

    class dpb_on_sho : public ghes_sho {
    public:
        dpb_on_sho() = default;
        dpb_on_sho(const double& _m, const double& _b, const double& _omega_d, const double& _omega_s,
            const int& _nbasis, const int& _nstate) : b(_b), omega_d(_omega_d),
            ghes_sho(_m, _omega_s, _nbasis, _nstate) { }

    private:
        const double b = 1; // dpb: minimum point
        const double omega_d = 1; // dpb: angular frequency

        double mel(const int& i, const int& j) override;
    };

    /*** ================================================================== ***/
    /*** Quartic Potential on Infinite Square Well Basis                    ***/
    /*** ================================================================== ***/

    class qua_on_isw : public ghes_isw {
    public:
        qua_on_isw() = default;
        qua_on_isw(const double& _m, const double& _b, const double& _omega_d, const double& _a,
            const int& _nbasis, const int& _nstate) : b(_b), omega_d(_omega_d),
            ghes_isw(_m, _a, _nbasis, _nstate) { }

    private:
        const double b = 1; // dpb: minimum point
        const double omega_d = 1; // dpb: angular frequency

        double mel(const int& i, const int& j) override;
    };

    /*** ================================================================== ***/
    /*** Quartic Potential on Simple Harmonic Oscillator Basis              ***/
    /*** ================================================================== ***/

    class qua_on_sho : public ghes_sho {
    public:
        qua_on_sho() = default;
        qua_on_sho(const double& _m, const double& _b, const double& _omega_d, const double& _omega_s,
            const int& _nbasis, const int& _nstate) : b(_b), omega_d(_omega_d),
            ghes_sho(_m, _omega_s, _nbasis, _nstate) { }

    private:
        const double b = 1; // dpb: minimum point
        const double omega_d = 1; // dpb: angular frequency

        double mel(const int& i, const int& j) override;
    };

    /*** ================================================================== ***/
    /*** Simple Harmonic Oscillator on Simple Harmonic Oscillator Basis     ***/
    /*** ================================================================== ***/

    class sho_on_sho : public ghes_sho {
    public:
        sho_on_sho() = default;
        sho_on_sho(const double& _m, const double& _omega_d, const double& _omega_s,
            const int& _nbasis, const int& _nstate) : omega_d(_omega_d),
            ghes_sho(_m, _omega_s, _nbasis, _nstate) { }

    private:
        const double omega_d = 1; // sho: angular frequency

        double mel(const int& i, const int& j) override;
    };

    /*** ================================================================== ***/
    /*** Morse Potential on Infinite Square Well Basis                      ***/
    /*** ================================================================== ***/

    class mrs_on_isw : public ghes_isw {
    public:
        mrs_on_isw() = default;
        mrs_on_isw(const double& _m, const double& _De, const double& _alpha, const double& _xeq,
            const double& _a, const int& _nbasis, const int& _nstate) : De(_De), alpha(_alpha), xeq(_xeq),
            ghes_isw(_m, _a, _nbasis, _nstate) { }

    private:
        const double De = 1; // mrs: depth
        const double alpha = 1; // mrs: alpha
        const double xeq = 1; // mrs: equilibrium position

        double mel(const int& i, const int& j) override;
    };

    /*** ================================================================== ***/
    /*** Eckart Potential on Infinite Square Well Basis                     ***/
    /*** ================================================================== ***/

    class ekt_on_isw : public ghes_isw {
    public:
        ekt_on_isw() = default;
        ekt_on_isw(const double& _m, const double& _V0, const double& _c, const double& _a,
            const int& _nbasis, const int& _nstate) : V0(_V0), c(_c),
            ghes_isw(_m, _a, _nbasis, _nstate) { }

    private:
        const double V0 = 1; // ekt: V_0
        const double c = 1; // ekt: c

        double mel(const int& i, const int& j) override;
    };

} // !hes
} // !uovie
#endif // !HES_EASY_H_