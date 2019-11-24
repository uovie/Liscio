/* Liscio: a concise computational chemistry program */

// Liscio headers
#include "prime/base.h"
#include "hes/easy.h"
#include "md/nve.h"
#include "density/rho.h"

using namespace uovie;

int main(int argc, char* argv[])
{
    prime::fairy Mavis;
    prime::init(argv[1], Mavis);

    std::pair<double, double> xlim; // x-axis limit

    auto exec_hes = [&xlim](hes::ghes& gensol) { // execute hes
        gensol.eigen_solve();
        gensol.print_evals();
        gensol.print_evecs(xlim);
    };

    if (Mavis.job.type == "hes") {
        if (Mavis.job.dscp[0] == "sho_on_isw") {
            hes::sho_on_isw soi(Mavis.job.para["m"], Mavis.job.para["omega_d"], Mavis.job.para["a"],
                Mavis.job.para["nbasis"], Mavis.job.para["nstate"]);
            const double a = Mavis.job.para["a"];
            xlim = std::pair<double, double>(-a / 2, a / 2);
            exec_hes(soi);
        }
        else if (Mavis.job.dscp[0] == "isw_on_sho") {
            hes::isw_on_sho ios(Mavis.job.para["m"], Mavis.job.para["a"], Mavis.job.para["omega_s"],
                Mavis.job.para["nbasis"], Mavis.job.para["nstate"]);
            const double a = Mavis.job.para["a"];
            xlim = std::pair<double, double>(-a, a);
            exec_hes(ios);
        }
        else if (Mavis.job.dscp[0] == "dpb_on_isw") {
            hes::dpb_on_isw doi(Mavis.job.para["m"], Mavis.job.para["b"], Mavis.job.para["omega_d"],
                Mavis.job.para["a"], Mavis.job.para["nbasis"], Mavis.job.para["nstate"]);
            const double a = Mavis.job.para["a"];
            xlim = std::pair<double, double>(-a / 2, a / 2);
            exec_hes(doi);
        }
        else if (Mavis.job.dscp[0] == "dpb_on_sho") {
            hes::dpb_on_sho dos(Mavis.job.para["m"], Mavis.job.para["b"], Mavis.job.para["omega_d"],
                Mavis.job.para["omega_s"], Mavis.job.para["nbasis"], Mavis.job.para["nstate"]);
            const double b = Mavis.job.para["b"];
            xlim = std::pair<double, double>(-3 * b, 3 * b);
            exec_hes(dos);
        }
        else if (Mavis.job.dscp[0] == "qua_on_isw") {
            hes::qua_on_isw doi(Mavis.job.para["m"], Mavis.job.para["b"], Mavis.job.para["omega_d"],
                Mavis.job.para["a"], Mavis.job.para["nbasis"], Mavis.job.para["nstate"]);
            const double a = Mavis.job.para["a"];
            xlim = std::pair<double, double>(-a / 2, a / 2);
            exec_hes(doi);
        }
        else if (Mavis.job.dscp[0] == "qua_on_sho") {
            hes::qua_on_sho qos(Mavis.job.para["m"], Mavis.job.para["b"], Mavis.job.para["omega_d"],
                Mavis.job.para["omega_s"], Mavis.job.para["nbasis"], Mavis.job.para["nstate"]);
            const double b = Mavis.job.para["b"];
            xlim = std::pair<double, double>(-3 * b, 3 * b);
            exec_hes(qos);
        }
        else if (Mavis.job.dscp[0] == "sho_on_sho") {
            hes::sho_on_sho sos(Mavis.job.para["m"], Mavis.job.para["omega_d"], Mavis.job.para["omega_s"],
                Mavis.job.para["nbasis"], Mavis.job.para["nstate"]);
            const double ep = 10; // endpoint
            xlim = std::pair<double, double>(-ep, ep);
            exec_hes(sos);
        }
        else if (Mavis.job.dscp[0] == "mrs_on_isw") {
            hes::mrs_on_isw moi(Mavis.job.para["m"], Mavis.job.para["De"], Mavis.job.para["alpha"],
                Mavis.job.para["xeq"], Mavis.job.para["a"], Mavis.job.para["nbasis"], Mavis.job.para["nstate"]);
            xlim = std::pair<double, double>(0, 20);
            exec_hes(moi);
        }
        else if (Mavis.job.dscp[0] == "ekt_on_isw") {
            hes::ekt_on_isw eoi(Mavis.job.para["m"], Mavis.job.para["V0"], Mavis.job.para["c"],
                Mavis.job.para["a"], Mavis.job.para["nbasis"], Mavis.job.para["nstate"]);
            xlim = std::pair<double, double>(-5, 5);
            exec_hes(eoi);
        }
    }
    else if (Mavis.job.type == "md") {
        if (Mavis.job.dscp[0] == "nve") {
            if (Mavis.job.dscp[1] == "ver") {
                md::md_nve_ver smd;
                smd.exec();
            }
        }
    }
	else if (Mavis.job.type == "density") {
		if (Mavis.job.dscp[0] == "sho_rho") {
			double T = Mavis.job.para["T"]; T /= prime::a_u_energy;
			density::sho_rho spho(Mavis.job.para["m"], Mavis.job.para["omega"], T);
			spho.calc_range(Mavis.job.para["begin"], Mavis.job.para["end"], Mavis.job.para["step"]);
		}
	}
    
    std::cout << "Normal termination. Congratulations!" << std::endl;
    return 0;
}