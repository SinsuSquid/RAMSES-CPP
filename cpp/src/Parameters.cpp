#include "ramses/Parameters.hpp"

namespace ramses {
namespace params {

    bool verbose = false;
    bool hydro = false;
    bool pic = false;
    bool poisson = false;
    bool cosmo = false;

    int nx = 1, ny = 1, nz = 1;
    int levelmin = 1;
    int nlevelmax = 1;
    int ngridmax = 0;
    real_t boxlen = 1.0;

    int nstepmax = 1000000;
    real_t trestart = 0.0;

} // namespace params
} // namespace ramses
