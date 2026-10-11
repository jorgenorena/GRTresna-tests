/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#ifndef MATTERPARAMS_HPP_
#define MATTERPARAMS_HPP_

#include "GRParmParse.hpp"
#include "REAL.H"

namespace MatterParams
{

struct params_t
{
    Real v;
    Real lambda;
    Real thickness;
    Real Ro;
};

inline void read_params(GRParmParse &pp, params_t &matter_params)
{
    Real compactness_o_Ro;
    pp.get("compactness_o_Ro", compactness_o_Ro);
    pp.get("thickness", matter_params.thickness);
    pp.get("Ro", matter_params.Ro);

    matter_params.v = sqrt((3.0/32.0/M_PI)*compactness_o_Ro*matter_params.thickness);
    matter_params.lambda = 2.0/(matter_params.v*matter_params.v)/(matter_params.thickness*matter_params.thickness);
}

}; // namespace MatterParams

#endif
