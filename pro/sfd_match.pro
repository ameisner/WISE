;+
; NAME:
;   sfd_match
;
; PURPOSE:
;   recalibrate WISE 12um ISSA tiles to SFD i100, assuming a constant
;   conversion factor between WISE DN and i100 intensity
;
; CALLING SEQUENCE:
;   
; INPUTS:
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
;   full sky SFD i100/WISE 12um hybrid results in
;   /n/wise/ameisner as wise-ebv-recalib.4096.fits, wise-i100-recalib.4096.fits
; REVISION HISTORY:
;   2012-Jun-19 - Written by Aaron Meisner
;
;----------------------------------------------------------------------
pro write_i100_tiles, outpath=outpath

;----- square, subsampled relative to wise tiles but still well-sampled
  IMPIX = 600 
  if ~keyword_set(outpath) then outpath = '/n/panlfs/ameisner/tile-sfd100'

  tstr = mrdfits('~ameisner/wise/pro/wisetile-index-allsky.fits', 1)
  tilepath = '/n/panlfs/ameisner/tile-allsky'
  flist = concat_dir(tilepath, tstr.fname)

  xbox = lindgen(IMPIX, IMPIX) MOD IMPIX
  ybox = lindgen(IMPIX, IMPIX) / IMPIX
  
  ntile = n_elements(flist)
  t0 = systime(1)
  for i=0,ntile-1 do begin
    print, i
    im = readfits(flist[i], h)
    hrebin, im, h, _, hsmall, OUTSIZE=[IMPIX, IMPIX]
    delvarx, im
    extast, hsmall, astr
    xy2ad, xbox, ybox, astr, abox, dbox
    euler, abox, dbox, lbox, bbox, 1
    i100 = dust_getval(lbox, bbox, map='i100', /interp, /noloop)
    outname = repstr(tstr[i].fname, 'wise', 'i100')
    writefits, concat_dir(outpath, outname), i100, hsmall
  endfor
  print, 'dt = ', systime(1) - t0, ' seconds'
end

; match to SFD i100 when smoothed to 6.1 arcmin. Assume SFD i100 has
; properly accounted for unresolved sources so that it is true diffuse emission
; that can be compared to star-subtracted WISE data.
; enhancement is factor by which to increase small angular scale WISE 
; fluctuations relative to their default scaling
function sfd_match, tnum, wisetile_path=wisetile_path, $ 
                     i100tile_path=i100tile_path, enhancement=enhancement, $ 
                     omask=omask, h=h

  if ~keyword_set(wisetile_path) then $ 
      wisetile_path = '/n/panlfs/ameisner/tile-allsky-flat'
  if ~keyword_set(i100tile_path) then $ 
      i100tile_path = '/n/panlfs/ameisner/tile-sfd100'
  if ~keyword_set(enhancement) then enhancement=1.

  wisefile = $ 
      concat_dir(wisetile_path, 'wise_'+string(tnum, format='(I03)')+'.fits')
  i100file = $ 
      concat_dir(i100tile_path, 'i100_'+string(tnum, format='(I03)')+'.fits')

  PIXSMALL = 600

;----- read in wise clean image
  wise = readfits(wisefile, h, /silent)
;----- put WISE DN into units of i100 MJy/sr
;      i have found in general that at high latitudes (and low
;      latitudes away from galactic center, wise ~ 3.5 x i100
  calfac = (1/3.5)*enhancement
  wise = wise*calfac
  wiseclean = wise ; save the full sized version

;----- use mask to interpolate over ghosts and latents
  omask = readfits(wisefile, ex=7, /silent)
  lpix = 41
  latentkern = shift(dist(lpix), lpix/2, lpix/2) LT (lpix/2+0.5)
  spix = 23
  satkern = shift(dist(spix), spix/2, spix/2) LT (spix/2+0.5)
  gpix = 11
  ghostkern = shift(dist(gpix), gpix/2, gpix/2) LT (gpix/2+0.5)
  mask_dilate = bytarr(size(wise, /DIM))
  mask_dilate OR= dilate((omask AND 8) NE 0, latentkern) ; latent = 2^3
  mask_dilate OR= dilate(omask AND 1, satkern) ; sat = 2^0
  mask_dilate OR= dilate((omask AND 2) NE 0, ghostkern) ; ghost = 2^1
  intx = djs_maskinterp(wise, mask_dilate, iaxis=0)
  inty = djs_maskinterp(wise, mask_dilate, iaxis=1)
  wise = (intx+inty)/2

;----- rebin to avoid doing unnecessary work during convolution
  wise = rebin(wise, PIXSMALL, PIXSMALL)
;----- create 6.1 arcmin gaussian smoothing kernel
  kpix = 21 ; 4+ fwhm
  kern = psf_gaussian(NPIXEL=kpix, FWHM=[4.88, 4.88], /NORM)
;----- smooth image
  smth = convol(wise, kern, /edge_zero)
  wt = convol(float(wise NE 0), kern, /edge_zero)
  smth = smth/(wt + (wt EQ 0))
;----- read in corresponding sfd i100 image (assume properly sized, 600x600)
  i100 = readfits(i100file, /silent)
  offset = i100-smth
  offset = rebin(offset, 3000, 3000)
  wisewt = readfits(wisefile, ex=2, /silent)
  wise_recalib = (wiseclean+ offset)*(wisewt NE 0) ; this really preferable???

  return, wise_recalib

end

pro match_i100_batch, indstart, nproc, outpath=outpath, _extra=extra

  if ~keyword_set(outpath) then outpath='/n/panlfs/ameisner/tile-recalib-i100'

  tstr = mrdfits('~ameisner/wise/pro/wisetile-index-allsky.fits', 1)
  ntile = n_elements(tstr)
  tnum = fix(strmid(tstr.fname, 5, 3))

  indend = (indstart + nproc - 1) < (ntile - 1)

  for i=indstart, indend do begin
      print, i
      recalib = sfd_match(tnum[i], omask=omask, h=h, _extra=extra)
      outname = concat_dir(outpath, tstr[i].fname)
      writefits, outname, recalib, h
; -------- modify header for extension
      sxaddpar, h, 'XTENSION', 'Image', ' IMAGE extension',before='SIMPLE'
      sxdelpar, h, 'SIMPLE'
      sxdelpar, h, 'EXTEND'
      writefits, outname, omask, h, /append
      delvarx, omask, h
  endfor
end
