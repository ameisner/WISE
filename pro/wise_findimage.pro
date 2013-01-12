;+
; NAME:
;   wise_findimage
;
; PURPOSE:
;   Find WISE Atlas images containing (RA, dec) positions
;
; CALLING SEQUENCE:
;   wise_findimage, ra, dec, wcs, ind1, ind2
;
; INPUTS:
;   (ra, dec)   - positions to match
;   wcs         - array of WISE WCS information
;  
; OUTPUTS:
;   ind1        - index list of objects
;   ind2        - index into wcs of WISE pointings
;
; EXAMPLES:
;   see wisephot.pro
;
; COMMENTS:
;   For each (ra, dec), find the WISE atlas image with the nearest
;     center, if there is one within 1.12 degrees.
;
; REVISION HISTORY:
;   2011-May-15 - Written by Douglas Finkbeiner, CfA
;
;----------------------------------------------------------------------
pro wise_findimage, ra, dec, wcs, ind1, ind2

; -------- match object list with WCS info from WISE
  spherematch, ra, dec, wcs.crval[0], wcs.crval[1], 1.12,  $
    m1, m2, d12, maxmatch=0

  sind = sort(m1)
  sm1 = m1[sind]
  nearest = bytarr(n_elements(sm1)) 
  u1 = uniq(sm1)
  j0 = 0L
  for i=0L, n_elements(u1)-1 do begin 
     j1 = u1[i]

     if (j0 EQ j1) then begin 
        nearest[j0] = 1B
     endif else begin 
        mindist = min(d12[sind[j0:j1]], thisind)
        nearest[j0+thisind] = 1B
     endelse

     j0 = j1+1
  endfor
  
  best = bytarr(n_elements(sm1))
  best[sind] = nearest
  ind = where(best, nbest)  ; index into m1, m2, d12

  ind1 = m1[ind]
  ind2 = m2[ind]

  return
end
