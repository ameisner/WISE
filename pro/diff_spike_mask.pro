;+
; NAME:
;   diff_spike_mask
;
; PURPOSE:
;   return a single star diffraction spike mask
;
; CALLING SEQUENCE:
;   msk = diff_spike_mask(mag, w4=w4)
;
; INPUTS:
;   mag - star magnitude
;
; KEYWORDS:
;   w4  - set for W4 (default W3)
;
; OUTPUTS:
;   msk - diffraction spike mask for star of interest (1=spike, 0=no spike)
;
; EXAMPLES:
;   
; COMMENTS:
;   w4 not yet implemented
;
; REVISION HISTORY:
;   2012-Nov-14 - Aaron Meisner
;----------------------------------------------------------------------
function diff_spike_mask, mag, m, w4=w4, sz=sz, dil=dil, i100rms=i100rms, $ 
                          radec=radec

  par = psf_par_struc(w4=w4, /everything)
  len = spike_half_length(mag, m, w4=w4, i100rms=i100rms, $ 
                          radec=radec)*(60/par.pscl) ; pix
  sz = 2*round(len/sqrt(2))+1
  sz = sz[0] ; make sure sz is not an array ...
  msk = bytarr(sz,sz)
  xbox = lindgen(sz,sz) MOD sz
  ybox = lindgen(sz,sz) / sz
  
  msk = (xbox EQ ybox) OR ((sz-1-xbox) EQ ybox)

; ----- at some point add check on whether dil is odd
  if keyword_set(dil) then begin
      kern = shift(dist(dil), dil/2, dil/2) LT (dil/2+0.5)
      msk = dilate(msk, kern)
  endif

  return, msk
end
