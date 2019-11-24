#ifndef DENSITY_PHO_H_
#define DENSITY_PHO_H_

#include "prime/base.h"

namespace uovie {
namespace density {

	class sho_rho {
	public:
		sho_rho() = default;
		sho_rho(const double& _m, const double& _omega, const double& _T) :
			m(_m), omega(_omega), T(_T) { }

		void calc_point(const double q);
		void calc_range(const double begin, const double end, const double step);

	private:
		const double m = 0;
		const double omega = 0;
		const double T = 0;
		const double beta = 1 / (prime::k * T);

		double rho = 0;
	};

} // !density
} // !uovie
#endif // !DENSITY_PHO_H_