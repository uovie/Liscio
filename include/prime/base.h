#ifndef PRIME_BASE_H_
#define PRIME_BASE_H_

// C++ Standard library
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <utility>
#include <chrono>

// Liscio headers
#include "math/func.h"

// JSON for Modern C++
#include <nlohmann/json.hpp>

// Eigen library
#include <Eigen/Dense>

// Boost library
#include <boost/math/special_functions/hermite.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>

namespace uovie {
namespace prime {

    // constants
    constexpr auto pi = M_PIl;
    constexpr double h_bar = 1;

    //------------------------------------------------------------------------//

    struct task {
        std::string type;
        nlohmann::json dscp;    // description
        nlohmann::json para;
    };

    struct fairy {
        task job;
    };

    //------------------------------------------------------------------------//

    template<typename T>
    char* as_bytes(T& i)
    {
        return reinterpret_cast<char*>(&i); // treat that memory as bytes
    }

    //------------------------------------------------------------------------//

    void init(const std::string& filename, fairy& Mavis);

} // !prime
} // !uovie
#endif // !PRIME_BASE_H_