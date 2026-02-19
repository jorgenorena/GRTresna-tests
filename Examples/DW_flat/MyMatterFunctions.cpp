/* GRTresna
 * Copyright 2024 The GRTL collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#include "ScalarField.hpp"

Real ScalarField::my_potential_function(const Real &phi_here) const
{
    /*
     * Symmetric potential
     * $$
     * V(\phi) = \frac{\lambda}{4}\left(\phi^2 - v^2\right)^2\,,
     * $$
     * Note that lambda is in units of the domain length
    */
    Real phi2_minus_v2 = (phi_here*phi_here - m_matter_params.v*m_matter_params.v);
    Real L = domainLength[0];
    Real lambda = m_matter_params.lambda/L/L;

    return 0.25*lambda*phi2_minus_v2*phi2_minus_v2;
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
    Real L = domainLength[0];
    return m_matter_params.v * tanh((1.0/m_matter_params.thickness/L) * (rr/m_matter_params.R_o*L - 1.0));
}

Real ScalarField::my_Pi_function(const RealVect &loc) const
{
    /* We just take it to be zero */
    return 0.0;
}
