;+
; NAME:
;   akari_find_nan
;
; PURPOSE:
;   Find NaNs in the 430 Akari tiles
;
; CALLING SEQUENCE:
;   akari_find_nan

; OUTPUTS:
;   list of files
;
; EXAMPLES:
;   just run it. 
;
; COMMENTS:
;   
; REVISION HISTORY:
;   2012-Feb-20 - Written by Douglas Finkbeiner, CfA
;
;----------------------------------------------------------------------
pro akari_find_nan

  t0 = systime(1)
  bands = ['WideS', 'WideL', 'N60', 'N160']

  for i=0, n_elements(bands)-1 do begin 
     band = bands[i]

     for tilenum=1, 430 do begin 
        
        tilestr = string(tilenum, format='(I3.3)')
        fname = '$AKARI_DATA/FITS/Release1.0/'+band+'/'+tilestr+'_'+band+'.fits'
        im = readfits(fname, h, /silent)
        
        winf = where(finite(im) EQ 0, ninf)
        if ninf GT 0 then begin 
           splog, 'NaNs found in file ', fname
           im[winf] = 0
        endif
        
     endfor
  endfor
  print, 'Elapsed time: ', systime(1)-t0
  
  return
end
