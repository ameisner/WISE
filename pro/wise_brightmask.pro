;+
; NAME:
;  wise_brightmask
;   
; PURPOSE:
;  make a bit mask for bright stars identifying saturated pixels,
;  bright pixels, and ghost-contaminated pixels
;
; CALLING SEQUENCE:  
;  outmsk = wise_brightmask(psf)
;
; INPUTS:
;  psf    - PSF image scaled to relevant star's flux
;
; OUTPUTS:
;  outmsk - bit mask for bright stars identifying saturated pixels, 
;           bright pixels, and ghost-contaminated pixels. See header
;           of wise_l1b_clean for bit details.
;
; REVISION HISTORY:
;   2012-Feb-19 - Written by Aaron Meisner
;----------------------------------------------------------------------
function wise_brightmask, psf, multi_epoch=multi_epoch, w4=w4, faint=faint

  if ~keyword_set(multi_epoch) then multi_epoch = 0
  par = psf_par_struc(w4=w4, /everything)
  PSFPIX = keyword_set(faint) ? par.pfaint : par.psfpix
  SAT    = 30000 ; DN

  BTHRESH = 250;threshold for pixel to be bright, DN

  outmsk = intarr(PSFPIX, PSFPIX)
  
  ; ---- pixel saturated -> 2^0
  outmsk OR= (psf GT SAT)

  ; ---- bright -> 2^2
  outmsk OR=  4*(psf GT BTHRESH)

  if ~keyword_set(faint) then begin
      BGHOST_THRESH = 20 ; DN
      GHOST_THRESH = 0.0001 ;flux relative to central PSF pixel
      indpsf = lindgen(PSFPIX, PSFPIX)
  ; ---- ghost -> 2^1
      outmsk OR= 2*(((indpsf / PSFPIX) LT par.ygpix) AND $ 
                    ((indpsf MOD PSFPIX) GE (PSFPIX/2-par.xgpix/2)) AND $ 
                    ((indpsf MOD PSFPIX) LE (PSFPIX/2+par.xgpix/2)) AND $
                     (psf/max(psf) GT GHOST_THRESH))
  ; ---- bright ghost
      outmsk OR= 32*(((outmsk AND 2) EQ 2) AND (psf GT BGHOST_THRESH))
  ; ---- mark ghost-affected pixels that can eventually be interpolated over,
  ;      as they inhabit locations on sky covered by other exposures in which
  ;      ghost has moved elsewhere
      if multi_epoch then outmsk OR= 128*((outmsk AND 2) NE 0)
  endif

  return, outmsk
end
