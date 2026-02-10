TOY2DAC – Modified regularization routines for edge-preserving FWI
=================================================================

License notice (TOY2DAC)
------------------------
Copyright 2008–2011 SEISCOPE project, All rights reserved.
Copyright 2013–2021 SEISCOPEII project, All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the SEISCOPE project nor the names of its contributors
   may be used to endorse or promote products derived from this software without
   specific prior written permission.

Warranty Disclaimer:
THIS SOFTWARE IS PROVIDED BY THE SEISCOPE PROJECT AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE SEISCOPE PROJECT OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.

-------------------------------------------------------------------

Purpose
-------
This repository provides modified/alternative regularization subroutines for use with
the 2D acoustic frequency-domain FWI code TOY2DAC (SEISCOPE Consortium), enabling
edge-preserving update conditioning inspired by rank-order/median regularization.

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
To generate custom filter coefficients, use the MATLAB script:
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

In addition, please cite the TOY2DAC software reference as requested by SEISCOPE
(see the TOY2DAC documentation and distribution page).
