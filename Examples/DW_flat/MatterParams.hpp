/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#ifndef MATTERPARAMS_HPP_
#define MATTERPARAMS_HPP_

#include "GRParmParse.hpp"
#include "REAL.H"
#include "DomainWallPotential.hpp"

namespace MatterParams
{

struct params_t
{
    Real v;
    Real lambda;
    Real thickness;
    Real Ro;
    Real q = 0.0;
    Real phi_true;
    Real phi_false;
};

inline void read_params(GRParmParse &pp, params_t &matter_params)
{
    Real compactness_o_Ro;
    pp.get("compactness_o_Ro", compactness_o_Ro);
    pp.get("thickness", matter_params.thickness);
    pp.get("Ro", matter_params.Ro);
    pp.load("q", matter_params.q, 0.0);
    if (!(std::isfinite(compactness_o_Ro) && compactness_o_Ro > 0.0 &&
          std::isfinite(matter_params.thickness) && matter_params.thickness > 0.0 &&
          std::isfinite(matter_params.Ro) && matter_params.Ro > 0.0))
        MayDay::Error("DW requires positive finite thickness, Ro and compactness_o_Ro");

    matter_params.v = sqrt((3.0/32.0/M_PI)*compactness_o_Ro*matter_params.thickness);
    matter_params.lambda = 2.0/(matter_params.v*matter_params.v)/(matter_params.thickness*matter_params.thickness);
    matter_params.phi_true = DomainWallPotential::true_vacuum(
        matter_params.v, matter_params.lambda, matter_params.q);
    matter_params.phi_false = DomainWallPotential::false_vacuum(
        matter_params.v, matter_params.phi_true);
    pout() << "DW potential: q=" << matter_params.q
           << " phi_true=" << matter_params.phi_true
           << " phi_false=" << matter_params.phi_false << std::endl;
}

}; // namespace MatterParams

#endif
