#include "density/rho.h"

namespace uovie {
namespace density {

	void sho_rho::calc_point(const double q)
	{
		int n_max = 20;

		double xi = sqrt(m * omega / prime::h_bar) * q;

		double left_term = (1 - exp(-beta * prime::h_bar * omega))
			* sqrt(m * omega / (prime::pi * prime::h_bar)) * exp(-pow(xi, 2));

		double right_term = 0;
		for (int n = 0; n < n_max; n++)
			right_term += exp(-n * beta * prime::h_bar * omega) / (pow(2, n) * math::factorial(n))
				* pow(boost::math::hermite(n, xi), 2);

		std::cout << left_term << "    " << right_term << std::endl;
		
		rho = left_term * right_term;
	}

	void sho_rho::calc_range(const double begin, const double end, const double step)
	{
		std::ofstream ocrd("crd.dat", std::ios_base::binary);
		std::ofstream oval("val.dat", std::ios_base::binary);

		for (double q = begin; q < end; q += step) {
			calc_point(q);
			ocrd.write(prime::as_bytes(q), sizeof(double));
			oval.write(prime::as_bytes(rho), sizeof(double));
		}

		ocrd.close();
		oval.close();
	}


} // !density
} // !uovie