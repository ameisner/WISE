;+
; NAME:
;   sfd_i100_rms
;
; PURPOSE:
;   generate healpix map containing RMS of SFD i100 within each pixel
;
; CALLING SEQUENCE:
;   rms = sfd_i100_rms(nside_out=, nside_samp=)
;
; OPTIONAL INPUTS:
;   nside_out  - nside of output SFD i100 RMS map
;   nside_samp - nside dictating the coordinates/number of samples to calculate
;                RMS (must have nside_samp > nside_out)
; OUTPUTS:
;   rms        - output RMS map, 1d array with length
;                12*(nside_out^2), units of MJy/Sr
;
; EXAMPLES:
;   
; COMMENTS:
;   there are a lot of combinations of nside_out, nside_samp
;   that are not sensible, but this code will not do any checks
;
; REVISION HISTORY:
;   2012-Nov-16 - Aaron Meisner
;----------------------------------------------------------------------
function sfd_i100_rms, nside_out=nside_out, nside_samp=nside_samp

  if ~keyword_set(nside_out) then nside_out = 64
  if ~keyword_set(nside_samp) then nside_samp = 512

  np_out = 12L*nside_out*nside_out
  np_samp = 12L*nside_samp*nside_samp
  rms = fltarr(np_out)

  healgen_lb, nside_samp, lsamp, bsamp
  ang2pix_ring, nside_out, (90.-bsamp)/!radeg, lsamp/!radeg, pix

  sind = sort(pix)
  pix = pix[sind]
  lsamp = lsamp[sind]
  bsamp = bsamp[sind]

  u = uniq(pix)

  i100 = dust_getval(lsamp, bsamp, map='i100', /interp, /noloop)
  for i=0L, np_out-1 do begin
      indl = (i EQ 0) ? 0 : (u[i-1]+1)
      indu = u[i]
      rms[i] = $ 
          sqrt(total((i100[indl:indu]-mean(i100[indl:indu]))^2)/(indu-indl+1))
  endfor

  return, rms
end
