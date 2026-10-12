# Cubic bias in the spherical wall

The current example reads `q` directly (default `0`) and solves with

    V(phi) = lambda/4 * (phi^2-v^2)^2 - q*phi^3 + V0,

where the global minimum has exactly zero potential. `thickness`, `Ro`,
`compactness_o_Ro`, `G_Newton`, and the grid retain their existing meanings.
`Ro` is a coordinate radius, not a fraction of the box length.
The thickness/compactness conversion assumes `G_Newton=1`.
At nonzero bias, `compactness_o_Ro` retains the unbiased-wall calibration;
it is not a measurement of the total compactness including the vacuum energy.

The initial wall interpolates from the false minimum at small radius to the
true minimum outside. For positive `q`, the true minimum is positive. Negative
`q` reverses the field orientation, retaining a zero-energy exterior. The
`q=0` branch preserves the original arithmetic and scalar profile.

`DomainWallPotential.hpp` evaluates a nonnegative factored expression; it does
not subtract two nearly equal vacuum energies or clip negative values. The
same header is maintained in GRChombo's `DW_flat` example. All local objects
depend on the example headers, since matter-parameter changes alter layouts.

Use identical `q`, `thickness` and `compactness_o_Ro` in the evolution input.
A completed nonlinear loop or an output HDF5 file does not establish adequate
initial constraints. Check the fresh final residuals and independently inspect
the written data. Reduce `q` if the solve does not converge; failed large-bias
solves do not establish physical non-collapse.

The workspace study in
`gr-dw-notes/results/bias-constraints-convergence-2026-10-11/` contains clean-build
regressions, independent potential/constraint checks, exact parameters and the
new evolution and convergence results.

The tested small-bias initial-data recipe uses `q=5`, `Ro=2`, `thickness=0.5`,
`compactness_o_Ro=0.25`, `G_Newton=1`, CTTK-Hybrid with `sign_of_K=1`, and an
octant of length 16. Fresh 128/192/256-cell solves with 60 nonlinear passes
all have final Hamiltonian
and momentum relative residuals below `1e-5 %`, finite fields, and positive
exported chi including ghosts. The matched `q=20` attempt fails these checks;
a decreasing `q=10` pilot misses the fixed-pass residual gate. These are
numerical decisions, not physical collapse thresholds.

Initial-data acceptance does not certify the later evolution. The matched
biased evolution series completes `t=8` but retains weak central momentum
convergence and boundary residuals; see the study report and the GRChombo
example's `README_BIAS.md` before using these data for quantitative physics.
