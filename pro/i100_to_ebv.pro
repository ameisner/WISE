;+
; NAME:
;   i100_to_ebv
;
; PURPOSE:
;   convert a WISE-enhanced SFD i100 tile into E(B-V) using tile bit masks
;   to interpolate where necessary and applying temperature correction
;
; CALLING SEQUENCE:
;   ebv = i100_to_ebv(tnum)
;
; INPUTS:
;   tnum     - integer (1-430) ISSA tile number
;
; KEYWORDS:
;   
;
; OUTPUTS:
;   ebv      - tile in units of E(B-V) and with artifacts interpolated over
;
; COMMENTS:
;   
; REVISION HISTORY:
;   2011-Jul-20 - Written by Aaron Meisner
;----------------------------------------------------------------------
function get_sfd, astr, wsfd, _extra=extra

  nx = astr.naxis[0]
  x = wsfd MOD nx
  y = wsfd / nx
  xy2ad, x, y, astr, a, d
  euler, a, d, l, b, 1
  sfd_i100 = dust_getval(l, b, /interp, /noloop, _extra=extra)
  return, sfd_i100

end

; _extra enables passing "enhancement" parameter to sfd_match
function i100_to_ebv, tnum, sfdfill=sfdfill, noneg=noneg, astr=astr, $ 
                      i100=i100, bitmask=bitmask, _extra=extra

  if (tnum LT 1) OR (tnum GT 430) then return, -1
  i100 = sfd_match(tnum, omask=omask, h=h, _extra=extra)
  extast, h, astr

  latmask = omask AND 8 ; latents
  ghostmask = omask AND 32 ; ghosts
  coremask = omask AND 4 ; bright stellar cores
  kpix =21
  kern = shift(dist(kpix), kpix/2, kpix/2) LT (kpix/2+0.5)

  bitmask = (8*dilate(latmask, kern) OR 32*dilate(ghostmask, kern) OR $ 
             4*dilate(coremask, kern))
  badmask = bitmask NE 0
  print, total(badmask)
  if ~keyword_set(sfdfill) then begin
      intx = djs_maskinterp(i100, badmask, iaxis=0)
      inty = djs_maskinterp(i100, badmask, iaxis=1)
      i100 = (intx+inty)/2
  endif else begin
; ----- fill in mask with SFD i100, how terrible is this idea???
      wbad = where(badmask, nwbad)
      if (nwbad GT 0) then i100[wbad] = get_sfd(astr, wbad)
  endelse

  if keyword_set(noneg) then begin
      wneg = where(i100 LT 0, nwneg)
      if (nwneg GT 0) then i100[wneg] = get_sfd(astr, wneg)
  endif

; ----- convert to E(b-v) according to SFD 98 prescription
  p = 0.0184
  X = get_sfd(astr, lindgen(astr.naxis[0], astr.naxis[1]), map='X')
  ebv = i100*X*p
  return, ebv

end
