;+
; NAME:
;   wise_write_tile
;
; PURPOSE:
;   Write a mosaic on ISSA tile number tilenum, 1..430
;
; CALLING SEQUENCE:
;   wise_write_tile, tilenum
;
; INPUTS:
;   tilenum  - ISSA tile number (1..430)
;
; OPTIONAL INPUTS:
;   
; KEYWORDS:
;   
; OUTPUTS:
;   
; OPTIONAL OUTPUTS:
;   
; EXAMPLES:
;   
; COMMENTS:
;   tilenum is 1-indexed, 1..430
;   
; REVISION HISTORY:
;   2011-Feb-10 - Written by Douglas Finkbeiner, CfA
;
;----------------------------------------------------------------------
pro wise_write_tile, tilenum, nimage=nimage, skip=skip, cleanpath=cleanpath, $ 
                     outpath=outpath, allsky=allsky, warp=warp, w4=w4

  if ~ keyword_set(cleanpath) then $ 
    cleanpath = '/n/panlfs/ameisner/clean.prelim.v0'
  print
  splog, 'ISSA tile', tilenum
  splog, 'Using cleanpath= ', cleanpath
  wait, 2

  index = '$AKARI_DATA/akari-index-WideS.fits'
  indstr = mrdfits(index, 1)

  fname = fileandpath(indstr.fname)
  tnum = strmid(fname, 0, 3)
  w = where(long(tilenum) EQ long(tnum), nw)
  if nw NE 1 then message, 'tile number is not in file '+index

  h = headfits('$AKARI_DATA/'+indstr[w].fname)
  extast, h, astr

  euler, indstr[w].ra, indstr[w].dec, lgal, bgal, 1
  lb = [lgal, bgal]

; -------- figure out how many images we need
  ind = wise_index_metadata(lb, angle=9.0, allsky=allsky, w4=w4)
  if ~keyword_set(ind) then begin
     splog, 'No images found for tile', tilenum, ' -- bailing!!!'
     return
  endif
  if ~keyword_set(nimage) then nimage = n_elements(ind) 

  if ~ keyword_set(skip) then begin 

; -------- write the cleaned files
     wise_write_clean1b, lb, nimage=nimage, outpath=cleanpath, allsky=allsky, $
         w4=w4

; -------- read files and map to healpix index grids
     wise_healpix_index, indh, flux, nside=2048, lb=lb, nimage=nimage, $ 
         cleanpath=cleanpath, indstr=ind, warp=warp, w4=w4
     
; -------- compute variance of pairwise image differences
     pair = wise_pairwise_compare(lb, indh, flux, allsky=allsky, w4=w4)
     
; -------- remove the worst outliers
     wise_pairwise_reject, pair
     
; -------- generate a mosaic
     goodlist = pair.ndiff GT 0
     
  endif else begin 
     goodlist = bytarr(nimage)+1B
  endelse
  wise_mosaic1b, im, raw, wt, minim, maxim, amask, omask, $
    dust=dust, astr=astr, hdr=hdr, goodlist=goodlist, cleanpath=cleanpath, $ 
    allsky=allsky, warp=warp, w4=w4

; -------- add more info to header
  band = keyword_set(w4) ? 4 : 3
  wav = keyword_set(w4) ? 22 : 12
  par = psf_par_struc(allsky=allsky, w4=w4, /everything)

  sxaddpar, hdr, 'EXTNAME', 'clean image'
  sxaddpar, hdr, 'TELESCOP', 'WISE'
  sxaddpar, hdr, 'BUNIT', 'DN'
  sxaddpar, hdr, 'WAVELEN', wav,     ' [micron] Passband center wavelength'
  sxaddpar, hdr, 'BAND', band,       ' Band number'
  sxaddpar, hdr, 'MAGZP', par.magzp, ' Magnitude zero point, Vega system'
  sxaddpar, hdr, 'NIMAGE',   nimage, ' Number of images used in mosaic'
  sxaddpar, hdr, 'COMMENT', 'Mosaic built from WISE level 1b images obtained from', after='MAGZP'
  sxaddpar, hdr, 'COMMENT', '  http://irsa.ipac.caltech.edu/Missions/wise.html'
  sxaddpar, hdr, 'COMMENT', 'Mosaic code by Aaron Meisner & Doug Finkbeiner'

; -------- subtract min, max images and reweight
  imclean = ((wt GT 2)*(wt*im-(minim+maxim)))/((wt-2) > 1)
  art = (wt GT 2)*(im-imclean)

  outname = 'wise_'+string(tilenum, format='(I3.3)')+'.fits'
  if keyword_set(outpath) then outname = concat_dir(outpath, outname)
  writefits, outname, imclean, hdr, /checksum

; -------- modify header for extensions
  sxaddpar, hdr, 'XTENSION', 'IMAGE', ' IMAGE extension',before='SIMPLE'
  sxdelpar, hdr, 'SIMPLE'
  sxdelpar, hdr, 'EXTEND'

  sxaddpar, hdr, 'EXTNAME', 'dirty image'
  writefits, outname, raw, hdr, /checksum, /append

  sxaddpar, hdr, 'EXTNAME', 'pixel weights'
; -------- writing weight of appropriate data type saves tens of MB/tile
  if (max(wt) LT 256) then begin
      wt = byte(wt)
      sxaddpar, hdr, 'BITPIX', 8
  endif else begin
      wt = fix(wt)
      sxaddpar, hdr, 'BITPIX', 16
  endelse
  writefits, outname, wt, hdr, /checksum, /append
  sxaddpar, hdr, 'BITPIX', -32

  hrebin, imclean, hdr, _, dhdr, OUTSIZE=size(imclean,/DIM)/2
  delvarx, _
  sxaddpar, dhdr, 'EXTNAME', 'SFD dust'
  writefits, outname, dust, dhdr, /checksum, /append

  sxaddpar, hdr, 'EXTNAME', 'minimum pixel value'
  writefits, outname, minim, hdr, /checksum, /append

  sxaddpar, hdr, 'EXTNAME', 'maximum pixel value'
  writefits, outname, maxim, hdr, /checksum, /append

; -------- header for uint masks
  bhdr = hdr
  sxaddpar, bhdr, 'BITPIX', 16
  sxaddpar, bhdr, 'BZERO', 32768, after='NAXIS2'

  sxaddpar, bhdr, 'EXTNAME', 'AND mask'
  writefits, outname, uint(amask), bhdr, /checksum, /append

  sxaddpar, bhdr, 'EXTNAME', 'OR mask'
  writefits, outname, uint(omask), bhdr, /checksum, /append

  sxaddpar, hdr, 'EXTNAME', 'artifact image'
  writefits, outname, art, hdr, /checksum, /append

  return
end
