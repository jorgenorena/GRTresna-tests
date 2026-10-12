/* GRTresna
 * Copyright 2024 The GRTL collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#include "ScalarField.hpp"

Real ScalarField::my_potential_function(const Real &phi_here) const
{
    return DomainWallPotential::value(phi_here, m_matter_params.v,
                                     m_matter_params.lambda, m_matter_params.q,
                                     m_matter_params.phi_true);
}

Real ScalarField::my_phi_function(const RealVect &loc) const
{
    /*
     * Spherical domain wall
     * $$
     * \phi = v\tanh\left(\sqrt{\frac{\lambda}{2}}v(r - R_o))\,.
     * $$
     * Note that R_o is in units of the domain length
    */
    Real rr = sqrt(loc[0] * loc[0] + loc[1] * loc[1] + loc[2] * loc[2]);
    const Real wall = tanh((1.0/m_matter_params.thickness) * (rr - m_matter_params.Ro));
    if (m_matter_params.q == 0.0)
        return m_matter_params.v * wall;
    // False-vacuum interior, zero-energy true-vacuum exterior. Negative q
    // reverses the field orientation, retaining an asymptotically flat exterior.
    return 0.5 * (m_matter_params.phi_true + m_matter_params.phi_false) +
           0.5 * (m_matter_params.phi_true - m_matter_params.phi_false) * wall;
}

Real ScalarField::my_Pi_function(const RealVect &loc) const
{
    /* We just take it to be zero */
    return 0.0;
}
