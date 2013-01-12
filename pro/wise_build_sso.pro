;+
; NAME:
;   wise_build_sso
;
; PURPOSE:
;   ingest WISE All-Sky Known Solar System Object Possible Association List
;   for use in W3, W4 image processing
;
; CALLING SEQUENCE:
;   wise_build_sso, cat
;   
; OUTPUTS:
;   cat - IDL structure containing relevant
;
; EXAMPLES:
;   wise_build_sso, cat
;   mwrfits, cat, '$WISE_DATA/sso_catalog.fits'
;
; COMMENTS:
;   long id added to make binary searching possible based on scan_id/frame_num
;
; REVISION HISTORY:
;   2012-Sep-11 - Written by Aaron Meisner
;
;----------------------------------------------------------------------
pro wise_build_sso, cat

  catfile = $ 
      '/n/panlfs/ameisner/wise_allsky.wise_allsky_4band_p1ba_mch15169.tbl'
  readcol, catfile, ra, dec, scan_id, frame_num, dra, ddec, w3mag, w4mag, $ 
      F='X, F, F, X, X, X, X, X, X, X, A, I, X, X, A, A, X, X, X, X, A, X, A',$
       skipline=67, /silent
; dra, ddec, w3mag, w4mag have nulls

  dra[where(dra EQ 'null')] = '-9999'
  ddec[where(ddec EQ 'null')] = '-9999'
  w3mag[where(w3mag EQ 'null')] = '-9999'
  w4mag[where(w4mag EQ 'null')] = '-9999'
  dra = float(dra)
  ddec = float(ddec)
  w3mag = float(w3mag)
  w4mag = float(w4mag)

  rawise = ra + (dra NE -9999)*dra/3600. ; this is the right sign convention
  decwise = dec + (ddec NE -9999)*ddec/3600.

  ncat = n_elements(ra)

; -----------  a    b    c    j
  suffuniq = ['a', 'b', 'c', 'j']
  suff2int = ['0', '1', '2', '9']
  suff = strmid(scan_id, 5, 1)
  matchlist, suffuniq, suff, mu, ms
  suff = strarr(ncat)
  suff[ms] = suff2int[mu]

  id = $ 
  long(strmid(scan_id, 0, 5)+suff+string(frame_num, format='(I03)'))
  cat = replicate({id:0L, ra_pred:0.,dec_pred:0., ra:0., dec:0., $ 
                   scan_id:'',frame_num:0, dra:0., ddec:0., w3mag:0., $ 
                   w4mag:0.}, ncat)

  cat.id = id
  cat.ra_pred = ra ; orbit-predicted ra
  cat.dec_pred = dec ; orbit-predicted dec
  cat.ra = rawise ; ra of WISE-detected counterpart
  cat.dec = decwise ; dec of WISE-detected counterpart
  cat.scan_id = scan_id
  cat.frame_num = frame_num
  cat.dra = dra
  cat.ddec = ddec
  cat.w3mag = w3mag
  cat.w4mag = w4mag

  sind = sort(id)
  cat = cat[sind]

end
