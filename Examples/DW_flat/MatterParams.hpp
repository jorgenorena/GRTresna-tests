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
    Real R_o;
};

inline void read_params(GRParmParse &pp, params_t &matter_params)
{
    pp.get("v", matter_params.v);
    pp.get("lambda", matter_params.lambda);
    pp.get("R_o", matter_params.R_o);
}

}; // namespace MatterParams

#endif
