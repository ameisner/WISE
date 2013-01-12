;+
; NAME:
;   wise_l1b_maskinterp
;
; PURPOSE:
;   use L1b mask files to interpolate over NaNs, CRs, and other bad pixels
;
; CALLING SEQUENCE:
;   intim = wise_l1b_maskinterp(raw, msk)
;
; INPUTS:
;   raw   - raw L1b image, i.e. image contained in *-w?-int-1b.fits
;   msk   - corresponding L1b mask file, i.e. *-w?-msk-1b.fits.gz
;   
; OUTPUTS:
;   intim - L1b intensity image with bad pixels replaced via interpolation
;
; OPTIONAL OUTPUTS:
;   badmask - mask indicating those pixels which have been
;             interpolated over based on raw L1b mask (1 = filled via
;             interpolation, 0 = unchanged)
;
; EXAMPLES:
;   see wise_l1b_clean
;
; COMMENTS:
;   
;
; REVISION HISTORY:
;   2012-Sep-22 - Code originally written by Doug in wise_l1b_clean
;                 moved here for the sake of modularity - Aaron Meisner
;
;----------------------------------------------------------------------
function wise_l1b_maskinterp, raw, msk, badmask=mask, nointerp=nointerp, $ 
                              satmask=satmask

; -------- dilate CR mask a bit
  crmask = (msk AND 2L^28) NE 0
  kpix   = 3     ; size of dilation kernel -- make this an odd number
  kern   = shift(dist(kpix), kpix/2, kpix/2) LT (kpix/2+0.5)
  crmask_dilate = dilate(crmask, kern)

; -------- interpolate over NaN mask, CR mask, and other bad pixels
  mask  = ((msk AND 255) NE 0) OR crmask_dilate OR (finite(raw) EQ 0)
  if arg_present(satmask) then begin
; -------- create mask for pixels saturated by bright signal, but not
;          part of static mask
      satmask = ((msk AND 523264) NE 0) AND ((msk AND 255) EQ 0)
      satmask = dilate(satmask NE 0, kern)
  endif
  if keyword_set(nointerp) then return, -1
  intx  = djs_maskinterp(raw, mask, iaxis=0, /const)
  inty  = djs_maskinterp(raw, mask, iaxis=1, /const)
  intim = (intx+inty)/2

  return, intim

end
