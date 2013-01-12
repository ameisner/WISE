;+
; NAME:
;   wise_check_clean
;
; PURPOSE:
;   check for missing or corrupt cleaned L1b images
;
; CALLING SEQUENCE:
;   wise_check_clean, indstart, nproc, allsky=, w4=, cleanpath=
;
; INPUTS:
;   indstart  - starting index within MJD-sorted index structure
;   nproc     - number of files to examine
;
; KEYWORDS:
;   allsky    - set for allsky release
;   w4        - set for W4 (default W3)
;   cleanpath - directory containing cleaned L1b images
;
; OUTPUTS:
;   you should pipe to idl to save the log file and later read/grep that
;   
; COMMENTS:
;   
; REVISION HISTORY:
;   2012-Dec-6 - Written by Aaron Meisner
;----------------------------------------------------------------------
pro wise_check_clean, indstart, nproc, allsky=allsky, w4=w4, $ 
                      cleanpath=cleanpath

  indstr = wise_index_metadata([0.,0.], /nosort, allsky=allsky, w4=w4)
  nfile = n_elements(indstr)

  fname = wise_name('clean', indstr.fname, cleanpath=cleanpath)
  indstart = long(indstart)
  indend = (indstart+nproc-1) < (nfile-1)

  t0 = systime(1)
  for i=indstart, indend do begin
      if (i MOD 100) EQ 0 then print, i
      if ~file_test(fname[i]) then begin
          print, 'missing: ', fname[i]
          continue
      endif
      im  = readfits(fname[i], /silent)
      dirt = readfits(fname[i], ex=1, /silent)
      msk = readfits(fname[i], ex=2, /silent)
      if (size(im, /TYPE) EQ 3) || (size(dirt, /TYPE) EQ 3) || $ 
         (size(msk, /TYPE) EQ 3) then print, 'corrupt file: ', fname[i]
  endfor

  print, 'dt = ', systime(1)-t0

end
