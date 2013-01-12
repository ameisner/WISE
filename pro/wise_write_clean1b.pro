;+
; NAME:
;   wise_write_clean1b
;
; PURPOSE:
;   Write many cleaned images near Galactic coordinates lb
;
; CALLING SEQUENCE:
;   wise_write_clean1b, lb, nimage=nimage
;
; INPUTS:
;   lb  - 2-element array, Galacitc coords [deg]
;
; KEYWORDS:
;   nimage - Number of images to process (default 2000)
;
; OUTPUTS:
;   a bunch of files is $WISE_DATA/clean
;
; OPTIONAL OUTPUTS:
;   
; EXAMPLES:
;   see wise_make_mosaic.pro
;
; COMMENTS:
;   Only works on files that are not there yet. 
;   Can run multiple threads and they will play nice. 
;
; REVISION HISTORY:
;   2011-Dec-01 - Written by Douglas Finkbeiner, CfA
;   2012-Jan-25 - Added outpath keyword, Aaron Meisner
;   2012-Feb-21 - Added extension to cleaned outputs containing star mask - AM
;   2012-Sep-16 - Call wise_clean_loop.pro for actual processing - AMM
;   2012-Sep-18 - Added keywords to enable parallel processing of exposures
;                 near specified pointing - Aaron Meisner
;
;----------------------------------------------------------------------
pro wise_write_clean1b, lb, nimage=nimage, outpath=outpath, allsky=allsky, $ 
                        indstart=indstart, nproc=nproc, w4=w4, _extra=extra

  if ~keyword_set(nimage) then nimage = 2000
  if ~keyword_set(indstart) then indstart = 0

  indend = keyword_set(nproc) ? $ 
      ((long(indstart)+nproc-1) < (nimage-1)) : (nimage-1)

; -------- resolve the environment variable here, because spawn,
;          /noshell below requires it. 
  if ~ keyword_set(outpath) then outpath = getenv('WISE_DATA') + '/clean'
  file_mkdir, outpath
  
; -------- read metadata table
  indstr = wise_index_metadata(lb, nimage=nimage, allsky=allsky, w4=w4)

  wise_clean_loop, indstr[indstart:indend], outpath=outpath, allsky=allsky, $ 
      w4=w4, _extra=extra

  return
end
