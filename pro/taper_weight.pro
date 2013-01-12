;+
; NAME:
;   taper_weight
;
; PURPOSE:
;   create a circularly symmetric weight mask that ramps linearly between
;   0 and 1
;
; CALLING SEQUENCE:
;   
; INPUTS:
;   szx  - x sidelength (pix) of desired output weight map
;   szy  - y sidelength (pix) of desired output weight map
;   r0   - radial coordinate at which weight becomes 0, must have r1 < r0
;   r1   - outermost radial coordinate for which weight remains exactly 1
;
; KEYWORDS:
;   
; OUTPUTS:
;   weight - weight image with weight in interval [0,1]
;
; OPTIONAL OUTPUTS:
;   
; EXAMPLES:
;   see wrapper taper_cutout.pro
;
; COMMENTS:
;   
; REVISION HISTORY:
;   2012-Nov-26 - Aaron Meisner
;----------------------------------------------------------------------
function taper_weight, szx, szy, r0, r1

; ----- assume szx, szy odd
  dist_ellipse, dell, [szx, szy], szx/2, szy/2, float(szy)/szx, 0.
  weight = 1+((dell-r1) < (r0-r1))*(-1./(r0-r1))*(dell GT r1)

  return, weight

end
