/* PRIME_BASE_CPP_ */

// C++ Standard library
#include <cstdlib>

// Lisce headers
#include "prime/base.h"

namespace uovie {
namespace prime {

    void read(const std::string& filename, nlohmann::json& load)
    {
        // check the extension
        if (filename.rfind(".json") == std::string::npos)
            throw std::invalid_argument("Only .json file is a valid input file.");
        std::ifstream in(filename);

        if (in.fail()) {
            std::cout << "Can not open the file " << filename << '.' << std::endl;
            exit(EXIT_FAILURE);
        }
        std::cout << "Read Infomation from " << filename << '.' << std::endl;
        in >> load;
        in.close(); // close files
    }

    void init(const std::string& filename, fairy& Mavis)
    {
        nlohmann::json load;
        read(filename, load);

        Mavis.job.type = load["job"]["type"];
        Mavis.job.dscp = load["job"]["dscp"];
        Mavis.job.para = load["job"]["para"];
    }

} // !prime
} // !uovie