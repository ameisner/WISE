;+
; NAME:
;   wise_subtract_other
;
; PURPOSE:
;   PSF subtract all but one particular source in an L1b image
;
; CALLING SEQUENCE:
;   subtract_other, im, h, x, y, allsky=allsky, starmask=starmask
;
; INPUTS:
;   im       - 1016x1016 WISE L1b image
;   h        - L1b header corresponding to im
;   x        - x location of source to exclude from subtraction process
;   y        - y location of source to exclude from subtraction process
;
; KEYWORDS:
;   allsky   - set for allsky release
;
; OPTIONAL OUTPUTS:
;   starmask - mask for compact source related artifacts
;
; COMMENTS:
;   the right way to do this is to eventually keep cntr identifier
;   from WISE catalog, but that causes the catalog to occupy more memory
;
; REVISION HISTORY:
;   2012-Oct-11 - Aaron Meisner
;----------------------------------------------------------------------
pro wise_subtract_other, im, h, x, y, allsky=allsky, starmask=starmask

  wise_starlist, h, xlist, ylist, maglist, allsky=allsky, mjdlim=mjdlim, wm=m
  nobj = n_elements(xlist)

  if (nobj LE 1) then begin 
      print, 'no other sources to subtract??'
      return
  endif

  dobj = sqrt((xlist-x)^2+(ylist-y)^2)
  dmin = min(dobj, indmin)

  if (dmin GT 0.) then begin
      print, dmin
      print, 'specific source of interest not found, giving up??'
      return
  endif

  keep = bytarr(nobj)+1
  keep[indmin] = 0
  wkeep = where(keep)
  srcstr = replicate({xlist:0.d, ylist:0.d, maglist:0.d, mjdlim:[0.d,0.d], $ 
                      m:0}, $ 
      nobj-1)
  srcstr.xlist = xlist[wkeep] 
  srcstr.ylist = ylist[wkeep]
  srcstr.maglist = maglist[wkeep]
  srcstr.mjdlim = transpose(mjdlim[wkeep,*])
  srcstr.m = m[wkeep]
  
  im = wise_l1b_substar(im, h, starmask=starmask, allsky=allsky, $ 
                        srcstr=srcstr)

end
