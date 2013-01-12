;+
; NAME:
;   dirbe
;
; PURPOSE:
;   recalibrate WISE 12um ISSA tiles to DIRBE 12um (band 5)
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
;   
; REVISION HISTORY:
;   2012-Jun-17 - Written by Aaron Meisner
;
;----------------------------------------------------------------------
pro dirbe_starlist, h, xlist, ylist

;----- mask size; 37 pix for 300x300 tiles, 361 pix for 3000 x 3000 tiles  
  MSKDIAM = 1.5 ; deg
  MSKPIX = round(MSKDIAM/abs(sxpar(h, 'CD1_1'))) + 1

  COMMON DIRBECATALOG, catra, catdec
  if n_elements(catra) EQ 0 then begin
    cat = mrdfits('$WISE_DATA/dirbe-psc.fits', 1)
    wbright = where(cat.flux5/cat.unc5 GT 2)
    catra = cat[wbright].ra
    catdec = cat[wbright].dec
  endif

  racen = sxpar(h, 'CRVAL1')
  deccen = sxpar(h, 'CRVAL2')
  XPIX = sxpar(h, 'NAXIS1')
  YPIX = sxpar(h, 'NAXIS2')

  maxsep = 12.5/sqrt(2) ; assume tiles 12.5 deg x 12.5 deg

  dangle = djs_diff_angle(catra, catdec, racen, deccen)
  wclose = where(dangle LE maxsep, nclose)
  if (nclose EQ 0) then return

  extast, h, astr
  ad2xy, catra[wclose], catdec[wclose], astr, x, y
  ix = round(x)
  iy = round(y)
  wtile = where((ix GE -MSKPIX/2) AND (ix LT XPIX+MSKPIX/2) AND $ 
                (iy GE -MSKPIX/2) AND (iy LT YPIX+MSKPIX/2), nwtile)
  if nwtile EQ 0 then return
  xlist = x[wtile]
  ylist = y[wtile]

end

pro test_dirbe_starlist, nstar, xmin, xmax, ymin, ymax, tilepath=tilepath, $ 
                         index=index

  if ~keyword_set(index) then $ 
    index = '~ameisner/wise/pro/wisetile-index-allsky.fits'
  tstr = mrdfits(index, 1)
  if ~keyword_set(tilepath) then tilepath = '/n/panlfs/ameisner/tile-allsky'
  
  flist = concat_dir(tilepath, tstr.fname)
  
  nstar = fltarr(n_elements(flist))
  xmax = -1e9
  ymax = -1e9
  xmin = 1e9
  ymin = 1e9
  for i=0, n_elements(flist)-1 do begin
    h = headfits(flist[i])
    dirbe_starlist, h, xlist, ylist
    nstar[i] = n_elements(xlist)
    if n_elements(xlist) EQ 0 then continue
    xmax >= max(xlist)
    ymax >= max(ylist)
    xmin <= min(xlist)
    ymin <= min(ylist)
    test_starmask = 1
    if test_starmask then mask = dirbe_star_mask(h, xlist, ylist)
    if (i MOD 10) EQ 0 then print, i
  endfor
end

pro starlist_single_tile, tnum, im, msk=msk

  if ((tnum LT 1) OR (tnum GT 430)) then return
  tilepath = '/n/panlfs/ameisner/tile-dirbe'
  tname = concat_dir(tilepath, 'dirbe_'+string(tnum, format='(I03)')+'.fits')
  
  im = readfits(tname, h)
  dirbe_starlist, h, xlist, ylist
  if n_elements(xlist) EQ 0 then begin
    print, 'no stars on this tile'
    return
  endif
  atv, tname ; want wcs
  if arg_present(msk) then msk = dirbe_star_mask(h, xlist, ylist)
  atvplot, xlist, ylist,ps=6,color='red',syms=2
end

function dirbe_star_mask, h, xlist, ylist

   IMPIX = sxpar(h, 'NAXIS1')
   msk = bytarr(sxpar(h, 'NAXIS1'), sxpar(h, 'NAXIS2'))
   STARSZ = round((9./6.)/abs(sxpar(h, 'CD1_1'))) + 1 ;37 pix for 300x300 tiles
                                                   ;361 pix for 3000x3000 tiles

;----- try circular rather than square mask
   kern = shift(dist(STARSZ), STARSZ/2, STARSZ/2) LT (STARSZ/2+0.5)

   nstar = n_elements(xlist)
   if (nstar EQ 0) then return, msk

   for i=0, nstar-1 do begin
     ix = round(xlist[i])
     iy = round(ylist[i])
     msk[((ix-STARSZ/2) > 0):((ix+STARSZ/2) < (IMPIX-1)), $ 
         ((iy-STARSZ/2) > 0):((iy+STARSZ/2) < (IMPIX-1))] OR= $ 
     kern[((STARSZ/2-ix) > 0):((IMPIX-ix+STARSZ/2-1) < (STARSZ-1)), $ 
          ((STARSZ/2-iy) > 0):((IMPIX-iy+STARSZ/2-1) < (STARSZ-1))]
   endfor
   return, msk

end

function dirbe_star_remove, im, h

  dirbe_starlist, h, xlist, ylist
  if n_elements(xlist) EQ 0 then return, im
  starmask = dirbe_star_mask(h, xlist, ylist)
  if total(starmask) EQ 0 then return, im
  intx = djs_maskinterp(im, starmask, iaxis=0)
  inty = djs_maskinterp(im, starmask, iaxis=1) ; use /const flag??
  clean = (intx+inty)/2

;----- remove residual pattern noise
  pix_kern = 7
  pix_fwhm = 2. ;precise nyquist value for being "well-sampled"??

  kern = psf_gaussian(NPIXEL=pix_kern, FWHM=[pix_fwhm, pix_fwhm], /NORM)
  
;----- pad the image
  IMPIX = sxpar(h, 'NAXIS1') ; assume square for now
  impad = fltarr(IMPIX+2*(pix_kern/2), IMPIX+2*(pix_kern/2))
  impad[pix_kern/2:pix_kern/2+IMPIX-1, pix_kern/2:pix_kern/2+IMPIX-1] = $ 
    clean*starmask
  smth = convol(impad, kern)
  wt = convol(float(impad NE 0), kern)
  filler = smth/(wt + (wt EQ 0))
  filler = filler[pix_kern/2:pix_kern/2+IMPIX-1, pix_kern/2:pix_kern/2+IMPIX-1]
  clean[where(starmask)] = filler[where(starmask)] ;prior check prevents failur
  return, clean

end

pro write_dirbe_tiles, outpath=outpath

  dirbe_heal = readfits('/n/panlfs/ameisner/dirbe5_1024.fits')
  IMPIX = 300 ; square
  if ~keyword_set(outpath) then outpath = '/n/panlfs/ameisner/tile-dirbe'

  xbox = lindgen(IMPIX, IMPIX) MOD IMPIX
  ybox = lindgen(IMPIX, IMPIX) / IMPIX

  tstr = mrdfits('~ameisner/wise/pro/wisetile-index-allsky.fits', 1)
  tilepath = '/n/panlfs/ameisner/tile-allsky'
  
  flist = concat_dir(tilepath, tstr.fname)

  ntile = n_elements(flist)
  t0 = systime(1)
  for i=0,ntile-1 do begin
    print, i
    im = readfits(flist[i], h)
    hrebin, im, h, _, hsmall, OUTSIZE=[IMPIX, IMPIX]
    extast, hsmall, astr
    xy2ad, xbox, ybox, astr, abox, dbox
    euler, abox, dbox, lbox, bbox, 1 ; dirbe healpix map in galactic projection
    pix_interp = heal_interp(dirbe_heal, lbox, bbox)
    im_interp = reform(pix_interp, IMPIX, IMPIX)
    dirbe_starlist, hsmall, xlist, ylist
    starmask = dirbe_star_mask(hsmall, xlist, ylist)
    clean = dirbe_star_remove(im_interp, hsmall)
    outname = repstr(tstr[i].fname, 'wise', 'dirbe')
    writefits, concat_dir(outpath, outname), im_interp, hsmall
; -------- modify header for extensions
    sxaddpar, hsmall, 'XTENSION', 'Image', ' IMAGE extension', before='SIMPLE'
    sxdelpar, hsmall, 'SIMPLE'
    sxdelpar, hsmall, 'EXTEND'
    writefits, concat_dir(outpath, outname), clean, hsmall, /append
    writefits, concat_dir(outpath, outname), starmask, hsmall, /append
  endfor
  print, 'dt = ', systime(1) - t0, ' seconds'

end

pro wise_match, indstart, nproc

  wisetile_path = '/n/panlfs/ameisner/tile-allsky-flat'
  dirbetile_path = '/n/panlfs/ameisner/tile-dirbe'
  outpath = '/n/panlfs/ameisner/tile-recalib'

  tstr = mrdfits('~ameisner/wise/pro/wisetile-index-allsky.fits', 1)
  flist_wise = tstr.fname
  flist_dirbe = repstr(flist_wise, 'wise', 'dirbe')
  ntile = n_elements(tstr)

  indend = (indstart + nproc - 1) < (ntile - 1)

  t0 = systime(1)
  for i=indstart, indend do begin
    print, i
;----- read in wise dirty extension, with stars + artifacts
    wise = readfits(concat_dir(wisetile_path, flist_wise[i]), ex=1 , /silent)
;----- make wise gain equal to that of DIRBE (pretending filters are identical)
    calfac = 0.0163402
    wise = wise*calfac

;----- subtract min/max artifacts
    artwise = readfits(concat_dir(wisetile_path, flist_wise[i]), ex=8, /silent)
    wise = wise-calfac*artwise
    
;----- use mask to interpolate over ghosts and latents
    omask = readfits(concat_dir(wisetile_path, flist_wise[i]), ex=7, /silent)
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

;----- rebin by 10 in each dimension to 300 x 300 
    wise = rebin(wise, 300, 300)

;----- create 0.7 deg gaussian PSF
    kpix = 51
    kern = psf_gaussian(NPIXEL=51, FWHM=[16.8, 16.8], /NORM)
;----- pad rebinned image
    IMPIX = 300
    impad = fltarr(IMPIX+2*(kpix/2), IMPIX+2*(kpix/2))
    impad[kpix/2:kpix/2+IMPIX-1, kpix/2:kpix/2+IMPIX-1] = wise
;----- convolve by gaussian to match dirbe resolution
    smth = convol(impad, kern)
    wt = convol(float(impad NE 0), kern)
    smth = smth/(wt + (wt EQ 0))
    smth = smth[kpix/2:kpix/2+IMPIX-1, kpix/2:kpix/2+IMPIX-1]
;----- find difference of dirbe and wise smoothed to dirbe resolution
    dirbe5 = readfits(concat_dir(dirbetile_path, flist_dirbe[i]), ex=1,/silent)
    offset = dirbe5-smth
    offset = rebin(offset, 3000, 3000)
;-----
    wiseclean = readfits(concat_dir(wisetile_path, flist_wise[i]), h, /silent)
    wisewt = readfits(concat_dir(wisetile_path, flist_wise[i]), ex=2, /silent)
    wise_recalib = ((wiseclean*calfac)+ offset)*(wisewt NE 0)
    outfile = concat_dir(outpath, flist_wise[i])
    writefits, outfile, wise_recalib, h
  endfor
  print, 'dt = ', systime(1) - t0, ' seconds'
end

pro append_dirbe_extensions

  tiledir =  '/n/panlfs/ameisner/tile-recalib'
  dirbedir = '/n/panlfs/ameisner/tile-dirbe'
  pushd, tiledir
  spawn, 'ls', flist
  popd

  flist_dirbe = repstr(flist, 'wise', 'dirbe')
  flist = concat_dir(tiledir, flist)  
  flist_dirbe = concat_dir(dirbedir, flist_dirbe)
  for i = 0, n_elements(flist)-1 do begin
    if (i MOD 10) EQ 0 then print, i
    dirbe_raw = readfits(flist_dirbe[i], hdirbe)
    dirbe_clean = readfits(flist_dirbe[i], ex=1)
    dirbe_mask = readfits(flist_dirbe[i], ex=2)
    writefits, flist[i], dirbe_raw, h, /append
    writefits, flist[i], dirbe_clean, h, /append
    writefits, flist[i], dirbe_mask, h, /append
  endfor
end
