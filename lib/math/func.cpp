/* MATH_FUNC_CPP_ */

// Liscio headers
#include "math/func.h"

namespace uovie {
namespace math {

    // factorial
    int factorial(int val) {
        int rst = 1;
        while (val > 1)
            rst *= val--;
        return rst;
    }

} // !math
} // !uovie