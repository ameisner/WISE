;+
; NAME:
;   wise_clean_batch
;
; PURPOSE:
;   Write many cleaned in some MJD range
;
; CALLING SEQUENCE:
;   wise_clean_batch, mjdrange, outpath=outpath
;
; INPUTS:
;   mjdrange  - range of MJDs to process
;
; OUTPUTS:
;   a bunch of files in outpath
;
; OPTIONAL OUTPUTS:
;   
; EXAMPLES:
;   see wise_make_mosaic.pro (sort of)
;
; COMMENTS:
;   Only works on files that are not there yet. 
;   Can run multiple threads and they will play nice. 
;
; REVISION HISTORY:
;   2011-Dec-01 - Written by Douglas Finkbeiner, CfA
;   2012-Jan-25 - Added outpath keyword, Aaron Meisner
;   2012-Feb-21 - Added extension to cleaned outputs containing star mask - AMM
;   2012-Feb-28 - Taken from wise_write_clean1b - DPF & AMM
;   2012-Sep-16 - Call wise_clean_loop.pro for actual processing - AMM
;
;----------------------------------------------------------------------
pro wise_clean_batch, mjdrange, outpath=outpath, allsky=allsky, $ 
                      indstart=indstart, nproc=nproc, w4=w4, _extra=extra

  if ~ keyword_set(mjdrange) then message, 'set mjdrange'
; ----- indstart is starting index within subset of exposures satisfying
;       mjdrange constraint
  if ~keyword_set(indstart) then indstart = 0

; -------- resolve the environment variable here, because spawn,
;          /noshell below requires it. 
  if ~ keyword_set(outpath) then outpath = getenv('WISE_DATA') + '/clean'
  file_mkdir, outpath
  
; -------- read metadata table
  ind = wise_index_metadata([0., 0.], /NOSORT, allsky=allsky, w4=w4)
  binary_search, ind.mjd, mjdrange[0], i0
  binary_search, ind.mjd, mjdrange[1], i1

  i0 >= 0
  i1 <= (n_elements(ind)-1) 
  if i1 EQ -1 then i1 = n_elements(ind)-1 
  ind = ind[i0:i1]
  nimage = n_elements(ind)
  indend = keyword_set(nproc) ? $ 
      ((long(indstart)+nproc-1) < (nimage-1)) : (nimage-1)

  wise_clean_loop, ind[indstart:indend], outpath=outpath, allsky=allsky, $ 
      w4=w4, _extra=extra

  return
end
