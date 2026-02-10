TOY2DAC – Modified regularization routines for edge-preserving FWI
=================================================================

Purpose
-------
This repository provides modified/alternative regularization subroutines for use with
the 2D acoustic frequency-domain FWI code TOY2DAC (SEISCOPE Consortium), enabling
edge-preserving update conditioning inspired by rank-order/median regularization.

Third-party software notice (TOY2DAC)
-------------------------------------
TOY2DAC is developed and distributed by the SEISCOPE / SEISCOPEII consortium.
License: FreeBSD-like (with diffusion restrictions). Please refer to the official
TOY2DAC distribution and the included 0LEGAL_STATEMENT/manual for full terms:
https://seiscope2.osug.fr/TOY2DAC-273

How to use (TOY2DAC integration)
-------------------------------
1) Choose a 1D or 2D window implementation and replace the corresponding TOY2DAC
   regularization subroutine with the provided "sub_Tikhonov.f90" (or the selected
   variant in this repository).
2) Recompile TOY2DAC.

Parameter settings
------------------
- For 1D windows (directional regularization):
  Set directional weights separately in the TOY2DAC "fwi_input" file:
  - lambda_x : horizontal weight
  - lambda_z : vertical weight

- For 2D window types (isotropic; no directional separation):
  Use either lambda_x or lambda_z > 0 (one is sufficient).
  If both are set > 0, the penalty may be effectively applied twice.

Filter coefficients (optional)
------------------------------
To generate custom filter coefficients, use the script:
  "filtCoefficients.m"
and update the selected 1D/2D regularization subroutine accordingly.

IMPORTANT:
The same coefficients must be defined consistently in both:
  - sub_Tikhonov_fgrad
  - sub_Tikhonov_Hv

Citation / Attribution
----------------------
If you use these routines in academic work, please cite:

Karcıoğlu, G., Tekkeli, A. B., Üge, M. A., & Arslan, M. S. (2023).
Inversion of VLF-R data using a non-linear rank order smoothing constraint and its
implementation to image surface rupture of the North Anatolian Fault in Başiskele,
Kocaeli, NW Turkey. Journal of Applied Geophysics, 216, 105145.

Also cite TOY2DAC as requested by SEISCOPE (see the TOY2DAC documentation).
