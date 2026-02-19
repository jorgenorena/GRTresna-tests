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
    Real R_o;
};

inline void read_params(GRParmParse &pp, params_t &matter_params)
{
    Real compactness, thickness;
    pp.get("compactness", compactness);
    pp.get("thickness", matter_params.thickness);
    pp.get("R_o", matter_params.R_o);

    matter_params.v = sqrt((3.0/32.0/M_PI)*compactness*thickness);
    matter_params.lambda = 2.0/(matter_params.v*matter_params.v)/(matter_params.R_o*matter_params.R_o)/(matter_params.thickness*matter_params.thickness);
}

}; // namespace MatterParams

#endif
