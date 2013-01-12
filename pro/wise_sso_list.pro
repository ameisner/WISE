;+
; NAME:
;   wise_sso_list
;
; PURPOSE:
;   identify SSOs in L1b exposure based on WISE SSO catalog
;
; CALLING SEQUENCE:
;   wise_sso_list, h, x, y, mag, w4=w4
;
; INPUTS:
;   h - header of L1b exposure of interest
;   
; KEYWORDS:
;   w4 - set for W4 (default W3)
;
; OUTPUTS:
;   x  - list of SSO x coordinates
;   y  - list of SSO y coordinates
;   mag - list of SSO magnitudes
;
; EXAMPLES:
;   see wise_remove_sso.pro
;
; COMMENTS:
;   meant to be analagous to wise_starlist.pro (although this situation
;   is different and more clear cut since each SSO catalog entry
;   corresponds identically to one exposure)
;
; REVISION HISTORY:
;   2012-Nov-14 - Written by Aaron Meisner
;----------------------------------------------------------------------
pro wise_sso_list, h, x, y, mag, w4=w4

  init_sso_catalog, w4=w4
  COMMON SSO, ra, dec, w3mag, w4mag, id, scan_id, frame_num

  this_scan_id = strtrim(sxpar(h, 'SCAN'))
  this_frame_num = fix(sxpar(h, 'FRNUM'))

  suffuniq = ['a', 'b', 'c', 'j']
  suff2int = ['0', '1', '2', '9']
  
  this_id = long(strmid(this_scan_id, 0, 5) + $
            suff2int[where(strmid(this_scan_id,5,1) EQ suffuniq)] + $ 
            string(this_frame_num, format='(I03)'))

  if (this_id GT max(id)) then return
  binary_search, id, this_id, indu
  binary_search, id, this_id-1, indl
  indl += 1
  nsso = indu-indl+1

  if (nsso EQ 0) then return
  w = indl + lindgen(nsso)
  thisra = ra[w]
  thisdec = dec[w]
  mag = keyword_set(w4) ? w4mag[w] : w3mag[w]
  extast, h, astr
  ad2xy, thisra, thisdec, astr, x, y

end
