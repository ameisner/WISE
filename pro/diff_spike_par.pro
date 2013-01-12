;+
; NAME:
;   diff_spike_par
;
; PURPOSE:
;   repository for diffraction spike masking related parameters
;
; CALLING SEQUENCE:
;   par = diff_spike_par(w4=)
;
; KEYWORDS:
;   w4 - set for W4 (default W3)
;
; OUTPUTS:
;   par - structure containing diffraction spike masking parameters
;
; EXAMPLES:
;   see diff_spike.pro and diff_spike_mask.pro
;
; COMMENTS:
;   
;
; REVISION HISTORY:
;   2012-Nov-20 - Aaron Meisner
;----------------------------------------------------------------------
function diff_spike_par, w4=w4

  w4 = keyword_set(w4)
  
; ----- polynomial coefficients for diffraction spike half length as a function
;       of coverage-corrected source brightness, entry i has units arcmin/mag^i
  coeff = w4 ? [ 5.711499, -1.3256876, -0.0475158, -0.48735428] : $
               [16.257549, -8.8601018,  2.0040373, -0.15745750]
; ----- faintest coverage-corrected source brightness contributing to
;       fit of coeff
  meff_max = w4 ? 1.50 : 4.25
; ----- brightest coverage-corrected source brightness contributing to
;       fit of coeff
  meff_min = w4 ? -3.30 : -2.20
; ----- integer coverage to which all magnitudes are corrected
  covmed = 16
; ----- w?m above which diffraction spike rescaling with m untested/undesired
  mmax = 60
; ----- factor by which to multiply log(i100rms) to effectively
;       convert to magnitudes
  i100fac = 1.1
; ----- L1b diffraction spike mask dilation kernel sidelength (L1b pixels)
  dil = 3
; ----- fiducial angular extent of single-epoch diffraction spikes, based
;       on examination of W3 PSF model (in radians)
  dtheta0 = 0.015 ; didn't consider better tailoring for W4

  par = { coeff         : coeff,      $
          meff_max      : meff_max,   $
          meff_min      : meff_min,   $ 
          covmed        : covmed,     $ 
          mmax          : mmax,       $
          i100fac       : i100fac,    $
          dil           : dil,        $
          dtheta0       : dtheta0      }

  return, par

end
