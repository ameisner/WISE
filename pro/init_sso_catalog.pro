;+
; NAME:
;   init_sso_catalog
;
; PURPOSE:
;   handle reading/filtering/caching of WISE SSO catalog
;
; CALLING SEQUENCE:
;   init_sso_catalog, w4=
;
; KEYWORDS:
;   w4   - set for W4 SSO catalog (default is W3)
;
; EXAMPLES:
;   see wise_remove_sso.pro
;
; COMMENTS:
;   
; REVISION HISTORY:
;   2012-Nov-11 - Aaron Meisner
;----------------------------------------------------------------------
pro init_sso_catalog, w4=w4

  COMMON SSO, ra, dec, w3mag, w4mag, id, scan_id, frame_num
  if n_elements(ra) EQ 0 then begin
    cat  = mrdfits('$WISE_DATA/sso_catalog.fits', 1)
    cat = keyword_set(w4) ? cat[where(cat.w4mag NE -9999)] : $ 
                            cat[where(cat.w3mag NE -9999)]
    ra = cat.ra
    dec = cat.dec
    w3mag = cat.w3mag
    w4mag = cat.w4mag
    id = cat.id
    scan_id = cat.scan_id
    frame_num = cat.frame_num
  endif

end
