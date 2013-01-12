;+
; NAME:
;   wise_exposure_warp
;
; PURPOSE:
;   warp a WISE exposure in attempt to correct for time-varying foregrounds
;
; CALLING SEQUENCE:
;   wise_exposure_warp, im, grad
;
; INPUTS:
;   im   - 500x500 cleaned L1b image
;   grad - [xgrad, ygrad] in units of WISE DN per (2.754*2) asec pixel
;   
; EXAMPLES:
;   see wise_mosaic1b.pro, wise_healpix_index.pro
;
; COMMENTS:
;   the initial version of this will be pretty trivial but I can
;   imagine the applied warps becoming more complex than simple
;   gradients and wise_mosaic1b is already very cluttered as is
;
;   
; REVISION HISTORY:
;   2011-Oct-12 - Aaron Meisner
;----------------------------------------------------------------------
pro wise_exposure_warp, im, grad, w4=w4

  par = psf_par_struc(w4=w4)
  cpix = par.pclean
  crpix = cpix/2-0.5
  sz = size(im, /DIM)

  if (sz[0] NE cpix) OR (sz[1] NE cpix) then begin
      print, 'unexpected image size'
      return
  endif

  xbox = lindgen(cpix,cpix) MOD cpix
  ybox = lindgen(cpix,cpix) / cpix

  warp = (xbox-crpix)*grad[0]+(ybox-crpix)*grad[1]
  im -= warp

end
