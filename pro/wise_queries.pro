;+
; NAME:
;   wise_queries
;
; PURPOSE:
;   Make query files for IRSA database query
;   
; CALLING SEQUENCE:
;   wise_queries
;
; KEYWORDS:
;   rarange    - range of ra values for plate centers to include
;   decrange   - range of dec values
; 
; OUTPUTS:
;   query files
;   
; EXAMPLES:
;   wise_queries,dec=[-3,3]
;
; COMMENTS:
;   Database queries can be uploaded here:
;     http://irsa.ipac.caltech.edu/applications/wise/
;
; REVISION HISTORY:
;   2011-Apr-16 - Written by Douglas Finkbeiner, CfA
;
;----------------------------------------------------------------------
pro wise_queries, nrow, rarange=rarange, decrange=decrange

  query_dir = '$WISE_DATA/queries'

  start = 0
  if ~ keyword_set(nrow) then nrow  = 500  ; split into chunks of 500 rows

; -------- read list of plate centers downloaded from web
  tilefile = '$WISE_DIR/etc/wise_tiles.dat'
  readcol, tilefile, ii, aa, ra, dec

; -------- select desired fields
  if keyword_set(rarange) then begin 
     wra = where(ra GE rarange[0] AND ra LE rarange[1], nra)
     if nra LE 0 then begin 
        message, 'OOPS -- no files satisfy RA range requirement'
     endif
     ra = ra[wra]
     dec = dec[wra]
  endif 

  if keyword_set(decrange) then begin 
     wdec = where(dec GE decrange[0] AND dec LE decrange[1], ndec)
     if ndec LE 0 then begin 
        message, 'OOPS -- no files satisfy dec range requirement'
     endif
     ra = ra[wdec]
     dec = dec[wdec]
  endif 

; -------- header for file
  head1 = '|    ra      |   dec       |'
  head2 = '|     r      |     r       |'

  nfield = n_elements(ra) 
  fieldno = start

  k = 1
  while fieldno LT nfield do begin 
     outfile = string('wisequery', k, '.dat', format='(A,I2.2,A)')
     outname = concat_dir(query_dir, outfile)
     print, 'Writing ', outname
     openw, wlun, outname, /get_lun
     printf, wlun, head1
     printf, wlun, head2
     for i = fieldno, ((fieldno+nrow) < nfield)-1 do $
       printf, wlun, ra[i], dec[i], format='(2F13.8)'
     free_lun, wlun
     fieldno += nrow
     k++
  endwhile
  

  return
end
