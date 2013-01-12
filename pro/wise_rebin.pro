;+
; NAME:
;   wise_rebin
;
; PURPOSE:
;   Pad and rebin wise images
;   
; CALLING SEQUENCE:
;   wise_rebin, image, h, binned, hout
;
; INPUTS:
;   image   - image to bin (assume atlas image, 4095x4095)
;   h       - header from atlas image file
;  
; OUTPUTS:
;   binned  - binned image (2048x2048)
;   hout    - output header, with WCS info updated
; 
; OPTIONAL OUTPUTS:
;   
; EXAMPLES:
;   
; COMMENTS:
;   
; REVISION HISTORY:
;   2011-May-17 - Written by Douglas Finkbeiner, CfA
;
;----------------------------------------------------------------------
pro wise_rebin, image, h, binned, hout, noextrap=noextrap

  nx = 2048   ; bin down to this size

; -------- determine image size, subtract median
  sz = size(image, /dimen)  ; will be [4095, 4095]

  npix = sz[0]+1  ; always 4096
; -------- assert that image is too small by 1. 
  if sz[0] NE (npix-1) or sz[1] NE (npix-1) then stop
  full = fltarr(npix, npix)
  full[0:npix-2, 0:npix-2] = image   

; -------- pad out to 4096 with linear interpolation
  if ~ keyword_set(noextrap) then begin 
     full[npix-1, *] = full[npix-2, *]*2 - full[npix-3, *]
     full[*, npix-1] = full[*, npix-2]*2 - full[*, npix-3]
  endif 

; -------- bin down 
  binned = rebin(temporary(full), nx, nx)

; -------- update astrometry after rebinning
  extast, h, astr

; -------- first deal with padding
  xy2ad, npix/2-0.5, npix/2-0.5, astr, racen, decen

  astr.crval = [racen, decen]
  astr.crpix += 0.5

; -------- now adjust for binning
;  astr.naxis = [1, 1]*nx
  astr.cdelt *= (4096/nx)
  astr.crpix = (astr.crpix-0.5)/(4096/nx)+0.5

; -------- overwrite WCS header and return
  hout = h
  putast, hout, astr
  sxaddpar, hout, 'NAXIS1', nx
  sxaddpar, hout, 'NAXIS2', nx
  sxaddpar, hout, 'PXSCAL1', 1.375d * (4096/nx)
  sxaddpar, hout, 'PXSCAL2', 1.375d * (4096/nx)
  sxaddpar, hout, 'HISTORY', 'WISE_REBIN.pro:  Rebinned by D. Finkbeiner, CfA '
  sxdelpar, hout, 'ELON'
  sxdelpar, hout, 'ELAT'
  sxdelpar, hout, 'GLON'
  sxdelpar, hout, 'GLAT'

  magzp0 = sxpar(hout, 'MAGZP')
  magzp = magzp0-2.5*alog10(4)
  print, magzp0, magzp
  
  sxaddpar, hout, 'MAGZP0', magzp0, ' [mag] original (pre-rebin) zero point', before='MAGZP'
  sxaddpar, hout, 'MAGZP', magzp, ' [mag] relative photometric zero point', after='MAGZP0'
  
  return
end


pro doit
  
  t0 = systime(1)
  input_dir  = 'L3a'
  output_dir = 'L3a_2048'
  flist = file_search(input_dir+'/*/*w?-int-3.fits', count=nfile)
  unclist = repstr(flist, '-int-', '-unc-')

  print, nfile, ' files' 
  for i=0L, nfile-1 do begin 
; -------- output file name and path
     file = fileandpath(flist[i], path=path)
     outpath = repstr(path, input_dir+'/', output_dir+'/')
     file_mkdir, outpath

     outfile = outpath+file
     uncfile = fileandpath(unclist[i], path=path)
     outunc  = outpath+uncfile
     
     print, i, '  Writing: ', outfile, format='(I5,A,A)'

; -------- read image and rebin
     if ~ file_test(outfile) then begin 
        if ~ file_test(flist[i]) then message, 'file not found: '+ flist[i], /info
        im = readfits(flist[i], h, /silent)
        
        if n_elements(im) GT 1 then begin
           wise_rebin, im, h, imout, hout
           writefits, outfile, imout, hout
        endif else begin 
           print, 'WARNING!!!  Short file!'
        endelse
     endif

; -------- read uncertainty and rebin
     if ~ file_test(outunc) then begin 
        if ~ file_test(unclist[i]) then message, 'file not found: '+ unclist[i], /info
        im = readfits(unclist[i], h, /silent)
        if n_elements(im) GT 1 then begin
           wise_rebin, im, h, imout, hout
           writefits, outunc, imout, hout
        endif else begin 
           print, 'WARNING!!!  Short file!'
        endelse
     endif
  endfor

  print, systime(1)-t0, ' seconds'

  return
end
