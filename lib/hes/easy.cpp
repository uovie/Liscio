/* HES_EASY_CPP_ */

// Liscio headers
#include "hes/easy.h"

namespace uovie {
namespace hes {

    /*** ================================================================== ***/
    /*** General Hamiltonian Eigen Solver                                   ***/
    /*** ================================================================== ***/

    void ghes::eigen_solve()
    {
        H.resize(nbasis, nbasis);
        for (int i = 0; i < nbasis; i++)
            for (int j = 0; j < nbasis; j++)
                H(i, j) = mel(i, j);
        es.compute(H);
    }

    void ghes::print_evals()
    {
        std::ofstream oval("val.dat", std::ios_base::binary);
        for (int i = 0; i < nstate; i++) {
            double val = es.eigenvalues()(i);
            oval.write(prime::as_bytes(val), sizeof(double));
        }
        oval.close();
    }

    void ghes::print_evecs(const std::pair<double, double>& xlim)
    {
        std::ofstream ocrd("crd.dat", std::ios_base::binary);
        for (double x = xlim.first; x <= xlim.second; x += 0.01)
            ocrd.write(prime::as_bytes(x), sizeof(double));
        for (int i = 0; i < nstate; i++) {
            std::ofstream ovec("vec_" + std::to_string(i) + ".dat", std::ios_base::binary);
            for (double x = xlim.first; x <= xlim.second; x += 0.01) {
                double efc = 0;
                for (int j = 0; j < nstate; j++)
                    efc += es.eigenvectors()(j, i) * bsf(j, x);
                ovec.write(prime::as_bytes(efc), sizeof(double));
            }
            ovec.close();
        }
    }

    /*** ================================================================== ***/
    /*** GHES with ISW Basis Function                                       ***/
    /*** ================================================================== ***/

    double ghes_isw::bsf(const int& i, const double& x)
    {
        return sqrt(2 / a) * sin((i + 1) * pi / a * (x + a / 2));
    }

    /*** ================================================================== ***/
    /*** GHES with SHO Basis Function                                       ***/
    /*** ================================================================== ***/

    double ghes_sho::bsf(const int& i, const double& x)
    {
        double xi = sqrt(m * omega_s / h_bar) * x;
        return pow(m * omega_s / (pi * h_bar), 0.25) / sqrt(pow(2, i) * math::factorial(i))
            * boost::math::hermite(i, xi) * exp(-pow(xi, 2) / 2);
    }

    /*** ================================================================== ***/
    /*** Simple Harmonic Oscillator on Infinite Square Well Basis           ***/
    /*** ================================================================== ***/

    double sho_on_isw::mel(const int& _i, const int& _j)
    {
        double i = _i + 1, j = _j + 1;
        return i == j ? pow(j * pi * h_bar, 2) / (2 * m * pow(a, 2)) + m * pow(omega_d * a, 2) / 4 * (1.0 / 6 - 1 / pow(j * pi, 2))
            : m * pow(omega_d * a, 2) / (2 * pow(pi, 2)) * ((cos((i - j) * pi) + 1) / pow(i - j, 2) - (cos((i + j) * pi) + 1) / pow(i + j, 2));
    }

    /*** ================================================================== ***/
    /*** Infinite Square Well on Simple Harmonic Oscillator Basis           ***/
    /*** ================================================================== ***/

    double isw_on_sho::mel(const int& i, const int& j)
    {
        return 0.5 * h_bar * omega_s * (i == j ? j + 0.5 : (i == j - 2 ? -0.5 * sqrt((j - 1) * j)
            : (i == j + 2 ? -0.5 * sqrt((j + 1) * (j + 2)) : 0)));
    }

    /*** ================================================================== ***/
    /*** Double Parabola on Infinite Square Well Basis                      ***/
    /*** ================================================================== ***/

    double dpb_on_isw::mel(const int& _i, const int& _j)
    {
        double i = _i + 1, j = _j + 1, rst = 0;
        if (i == j) rst += pow(j * pi * h_bar, 2) / (2 * m * pow(a, 2));
        
        auto ligd = [this, &i, &j](double x) {
            return sin(i * pi / a * (x + a / 2)) * pow(x + b, 2) * sin(j * pi / a * (x + a / 2));
        }; // left integrand
        double lint = boost::math::quadrature::trapezoidal(ligd, -a / 2, 0.0, 1e-9); // left integral

        auto rigd = [this, &i, &j](double x) {
            return sin(i * pi / a * (x + a / 2)) * pow(x - b, 2) * sin(j * pi / a * (x + a / 2));
        }; // right integrand
        double rint = boost::math::quadrature::trapezoidal(rigd, 0.0, a / 2, 1e-9); // right integral

        rst += m * pow(omega_d, 2) / a * (lint + rint);
        return rst;
    }

    /*** ================================================================== ***/
    /*** Double Parabola on Simple Harmonic Oscillator Basis                ***/
    /*** ================================================================== ***/

    double dpb_on_sho::mel(const int& i, const int& j)
    {
        using namespace boost::math;
        using namespace boost::math::quadrature;

        double rst = 0;
        if (i == j) rst += 0.5 * h_bar * omega_s * (j + 0.5);

        auto ligd = [this, &i, &j](double xi) {
            return pow(sqrt(h_bar / (m * omega_s)) * xi + b, 2) * hermite(i, xi) * hermite(j, xi) * exp(-pow(xi, 2));
        }; // left integrand
        double lint = boost::math::quadrature::trapezoidal(ligd, -10.0, 0.0, 1e-9); // left integral

        auto rigd = [this, &i, &j](double xi) {
            return pow(sqrt(h_bar / (m * omega_s)) * xi - b, 2) * hermite(i, xi) * hermite(j, xi) * exp(-pow(xi, 2));
        }; // right integrand
        double rint = boost::math::quadrature::trapezoidal(rigd, 0.0, 10.0, 1e-9); // right integral

        rst += m * pow(omega_d, 2) / sqrt(pow(2, i + j + 2) * math::factorial(i) * math::factorial(j) * pi) * (lint + rint);
        return rst;
    }

    /*** ================================================================== ***/
    /*** Quartic Potential on Infinite Square Well Basis                    ***/
    /*** ================================================================== ***/

    double qua_on_isw::mel(const int& _i, const int& _j)
    {
        double i = _i + 1, j = _j + 1, rst = 0;
        if (i == j) rst += pow(j * pi * h_bar, 2) / (2 * m * pow(a, 2));

        auto cigd = [this, &i, &j](double x) {
            return sin(i * pi / a * (x + a / 2)) * pow(pow(x, 2) - pow(b, 2), 2) * sin(j * pi / a * (x + a / 2));
        }; // complete integrand
        double cint = boost::math::quadrature::trapezoidal(cigd, -a / 2, a / 2, 1e-9); // complete integral

        rst += m * pow(omega_d, 2) / (4 * a * pow(b, 2)) * cint;
        return rst;
    }

    /*** ================================================================== ***/
    /*** Quartic Potential on Simple Harmonic Oscillator Basis              ***/
    /*** ================================================================== ***/

    double qua_on_sho::mel(const int& i, const int& j)
    {
        using namespace boost::math;
        using namespace boost::math::quadrature;

        double rst = 0;
        if (i == j) rst += 0.5 * h_bar * omega_s * (j + 0.5);

        auto cigd = [this, &i, &j](double xi) {
            return pow(h_bar / (m * omega_s) * pow(xi, 2) - pow(b, 2), 2) * hermite(i, xi) * hermite(j, xi) * exp(-pow(xi, 2));
        }; // complete integrand
        double cint = boost::math::quadrature::trapezoidal(cigd, -10.0, 10.0, 1e-9); // complete integral

        rst += m * pow(omega_d, 2) / (sqrt(pow(2, i + j + 6) * math::factorial(i) * math::factorial(j) * pi) * pow(b, 2)) * cint;
        return rst;
    }

    /*** ================================================================== ***/
    /*** Simple Harmonic Oscillator on Simple Harmonic Oscillator Basis     ***/
    /*** ================================================================== ***/

    double sho_on_sho::mel(const int& i, const int& j)
    {
        using namespace boost::math;
        using namespace boost::math::quadrature;

        double rst = 0;
        if (i == j) rst += 0.5 * h_bar * omega_s * (j + 0.5);

        auto cigd = [this, &i, &j](double xi) {
            return pow(xi, 2) * hermite(i, xi) * hermite(j, xi) * exp(-pow(xi, 2));
        }; // complete integrand
        double cint = boost::math::quadrature::trapezoidal(cigd, -10.0, 10.0, 1e-9); // complete integral

        rst += h_bar * pow(omega_d, 2) / (sqrt(pow(2, i + j + 2) * math::factorial(i) * math::factorial(j) * pi) * omega_s) * cint;
        return rst;
    }

    /*** ================================================================== ***/
    /*** Morse Potential on Infinite Square Well Basis                      ***/
    /*** ================================================================== ***/

    double mrs_on_isw::mel(const int& _i, const int& _j)
    {
        double i = _i + 1, j = _j + 1, rst = 0;
        if (i == j) rst += pow(j * pi * h_bar, 2) / (2 * m * pow(a, 2));

        auto cigd = [this, &i, &j](double x) {
            return sin(i * pi / a * x) * pow(1 - exp(-alpha * (x - xeq)), 2) * sin(j * pi / a * x);
        }; // complete integrand
        double cint = boost::math::quadrature::trapezoidal(cigd, 0.0, a, 1e-9); // complete integral

        rst += 2 * De / a * cint;
        return rst;
    }

    /*** ================================================================== ***/
    /*** Eckart Potential on Infinite Square Well Basis                     ***/
    /*** ================================================================== ***/

    double ekt_on_isw::mel(const int& _i, const int& _j)
    {
        double i = _i + 1, j = _j + 1, rst = 0;
        if (i == j) rst += pow(j * pi * h_bar, 2) / (2 * m * pow(a, 2));

        auto cigd = [this, &i, &j](double x) {
            return sin(i * pi / a * (x + a / 2)) / pow(cosh(c * x), 2) * sin(j * pi / a * (x + a / 2));
        }; // complete integrand
        double cint = boost::math::quadrature::trapezoidal(cigd, -a / 2, a / 2, 1e-9); // complete integral

        rst -= 2 * V0 / a * cint;
        return rst;
    }

} // !hes
} // !uovie