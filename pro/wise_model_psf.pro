;+
; NAME:
;   wise_model_psf
; PURPOSE:
;   gather cutouts and other information necessary to construct PSF
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
;   NB: the x, y, ra, dec coordinates in cutout summary structures
;   ARE THOSE OF THE STAR to which the feature being modeled pertains, 
;   not necessarily the coordinates of the feature itself, for example
;   in the cases of latents/ghosts. however, fname does rever to
;   the image from which cutout was extracted, not necessarily the
;   image on which relevant bright star was detected.
;
;
; REVISION HISTORY:
;   2012-Oct-15 - Aaron Meisner
;----------------------------------------------------------------------
function sshift_bitmask, im, xyshift, wt

; sinc shift a rectangular subsection of an image based on a coverage bitmask 
; that indicates where data is present (1) or absent (0)

  wshift = where(wt, nwshift)
  if nwshift EQ 0 then begin
     print, 'entire image missing??'
     return, -1
  endif

  xsz = (size(im, /DIM))[0]
  xwshift = wshift MOD xsz
  ywshift = wshift / xsz

  shifted = im
; ----- don't try to shift something with sidelength=1 along one or
;       more dimensions
  if ((max(xwshift)-min(xwshift)) EQ 0) OR $ 
     ((max(ywshift)-min(ywshift)) EQ 0) then return, shifted
  shifted[min(xwshift):max(xwshift),min(ywshift):max(ywshift)] = $
      sshift2d(im[min(xwshift):max(xwshift),min(ywshift):max(ywshift)], $ 
               xyshift)

  return, shifted

end

function determine_background, im, x, y, w4=w4, feat=feat

  if ~keyword_set(feat) then feat='wings'
  par = psf_par_struc(w4=w4, feat=feat)
  bg = djs_photsky(x, y, par.skyrad, im)
  return, bg

end

function check_mask, msk, ix, iy, w4=w4, feat=feat

; sometimes bit masks can flag real PSF features as spurious "spikes" e.g.
; in W4 PSF "core", with 2^28 bit set
; ix, iy must be integer x,y coords of the feature being modeled
; return 1 if everything is OK, return 0 if not OK

  if ~keyword_set(w4) OR (feat NE 'core') then return, 1

; ---- examine bitmask within a 3x3 rectangle about ix,iy
  sz = 3
  incl = wise_l1b_cutout(_, ix, iy, sz, sz, w4=w4, /bool, /silent)
  if ~incl then return, 1
  mcut = long(wise_l1b_cutout(msk, ix, iy, sz, sz, w4=w4, /silent))
  good = (total((mcut AND 2L^28) NE 0) EQ 0)
  return, good

end

function prepare_single_cutout, im, msk, unc, h, x, y, bg=bg, feat=feat, $ 
                                subtract=subtract, wt=wt, w4=w4

; given file name and compact source coordinates, retrieve relevant cutout
; and corresponding coverage bitmask

  if ~keyword_set(feat) then feat='wings'
 ; full = (feat EQ 'latent') ; demand full cutout ?
  full = 0
  par = psf_par_struc(w4=w4, feat=feat)

  mgood = check_mask(msk, round(x), round(y)-par.yoffs, w4=w4, feat=feat)
  if ~mgood then return, -1

  intim = wise_l1b_maskinterp(im, msk)
  intunc = wise_l1b_maskinterp(unc, msk)
  if keyword_set(subtract) then $ 
      wise_subtract_other, intim, h, x, y, /allsky, starmask=starmask

  xpix = par.szx
  ypix = par.szy

  bg = determine_background(intim, x, y, w4=w4, feat=feat)
  cutout = wise_l1b_cutout(intim, round(x), round(y)-par.yoffs, xpix, ypix, $ 
                           wt=wt, w4=w4, full=full)
  if max(wt) EQ 0 then return, -1
  unc = wise_l1b_cutout(intunc, round(x), round(y)-par.yoffs, xpix, ypix, $ 
                           wt=wt, w4=w4)
  cutout -= bg*(wt NE 0)
  cutout = sshift_bitmask(cutout, [round(x)-x, round(y)-y], wt)
  unc = sshift_bitmask(unc, [round(x)-x, round(y)-y], wt)

  return, cutout

end

function package_cutout, ra, dec, mag, x, y, fname, cutout, wt, unc, bg

  outstr = { ra        :  ra,     $ 
             dec       :  dec,    $
             mag       :  mag,    $
             x         :  x,      $
             y         :  y,      $
             fname     :  fname,  $ 
             cutout    :  cutout, $ 
             wt        :  wt,     $
             unc       :  unc,    $ 
             bg        :  bg        }

  return, outstr

end

function prepare_relevant_cutouts, ra, dec, mag, nmax=nmax, feat=feat, $ 
                                   subtract=subtract, w4=w4

; given ra,dec determine images from which cutouts can be drawn
; and call prepare_single_cutout to retrieve

  if ~keyword_set(feat) then feat='wings'
  if ~keyword_set(nmax) then nmax = 10000L

  par = psf_par_struc(w4=w4, feat=feat)
  xpix = par.szx
  ypix = par.szy

  maxsep = $
    sqrt((par.crpix+xpix/2)^2+(par.crpix+ypix/2+par.yoffs)^2)
; ----- convert from pixels to deg
  maxsep = maxsep*(par.pscl)/3600.

  euler, ra, dec, lgal, bgal, 1
  lb = [lgal, bgal]
  indstr = wise_index_metadata(lb, angle=maxsep, /allsky, w4=w4)
  nfile = n_elements(indstr)
  print, 'number of candidate L1b exposures: ', nfile
  if nfile EQ 0 then return, -1
  
  indend = keyword_set(all) ? (nfile-1) : (nfile-1) < (nmax-1)
  for i=0, indend do begin
      print, 'investigating cutout ', i, ' of ', (indend+1)
      thisfile = indstr[i].fname
      h = headfits(thisfile, /silent)
      extast, h, astr
      ad2xy, ra, dec, astr, x, y 
; ----- is the cutout "included" in this image or not?
      incl = $ 
          wise_l1b_cutout(_, round(x), round(y)-par.yoffs, xpix, ypix, $ 
                          wt=wt, w4=w4, /bool)
      if ~incl then continue
      if feat EQ 'latent' then begin
          prev = previous_exposure(indstr[i].mjd, /allsky, w4=w4, /forward)
          if (size(prev, /TYPE) EQ 2) then continue
          h = headfits(prev[0].fname, /silent)
          thisfile = prev[0].fname
      endif
      if (feat EQ 'latent2') then begin
          prev = previous_exposure(indstr[i].mjd,/allsky,w4=w4,/forward,N=2)
          if (size(prev, /TYPE) EQ 2) then continue
          w2 = where(prev.n EQ 2)
          if (w2 EQ -1) then continue
          h = headfits(prev[w2].fname, /silent)
          thisfile = prev[w2].fname
      endif
      im = readfits(thisfile, /silent)
      mname = repstr(thisfile, 'int', 'msk') + '.gz'
      msk = readfits(mname, /silent)
      uname = repstr(thisfile, 'int', 'unc') + '.gz'
      unc = readfits(uname, /silent)
; ----- don't crash on missing uncertainty image
      if n_elements(unc) EQ 1 then continue
      cutout = $
          prepare_single_cutout(im, msk, unc, h, x, y, subtract=subtract, $ 
                                wt=wt, w4=w4, bg=bg, feat=feat)
      if max(wt) EQ 0 then continue
      cutstr = package_cutout(ra, dec, mag, x, y, thisfile, cutout, $ 
                              wt, unc, bg)
      if n_elements(cubestr) EQ 0 then cubestr = cutstr else $
          cubestr = struct_append(cubestr, cutstr)
  endfor

  if n_elements(cubestr) EQ 0 then return, -1
  return, cubestr

end

pro get_star_sample, ra, dec, mag, w4=w4, feat=feat

  if ~keyword_set(feat) then feat = 'wings'
  if ~keyword_set(w4) then begin
; ----- read in bright star custom fitted mags, coords, and associated
;       flags
      if (feat EQ 'core') OR (feat EQ 'latent') then begin
           init_source_catalog, w4=w4, /allsky
           COMMON CATALOG, catra, catdec, catmag, mjdmin, mjdmax, catm
           good = bytarr(n_elements(catra))+1
      endif else begin
           _ = $ 
             analyze_fit_metrics(catra, catdec, catmag, fpath=fpath, good=good)
      endelse 
           
      euler, catra, catdec, lgal, bgal, 1
      mag_upper = (feat EQ 'core') ? 5. :  1.
      mag_lower = (feat EQ 'core') ? 4. : -2.
      bmin = (feat EQ 'core') ? 40 : 15
      if (feat EQ 'latent') then begin
          mag_upper = 4 & mag_lower = -3
          bmin = (15 > (15+(catmag-1.)*15.)) < 40
      endif
      if (feat EQ 'latent2') then begin
          mag_upper = -0.5
          mag_lower = -2.0
      endif
; ----- eventually absorb these mag/coord cut parameters into psf_par_struc.pro
      wpsf = where((catmag LT mag_upper) AND (catmag GT mag_lower) AND $ 
                   (abs(bgal) GT bmin) AND (good))
      ra = catra[wpsf]
      dec = catdec[wpsf]
      mag = catmag[wpsf]
  endif else begin
      cat = mrdfits('$WISE_DATA/w4_catalog-allsky.fits', 1)
      euler, cat.ra, cat.dec, lgal, bgal, 1
; ----- eventually absorb these into psf_par_struc
      mag_upper = (feat EQ 'core') ? 0.5  :   0.0
      mag_lower = (feat EQ 'core') ? 0.0  :  -5.0
      if (feat EQ 'latent2') then begin
         mag_upper = -1.
         mag_lower = -5.
      endif
      wpsf = where((cat.w4mpro LT mag_upper) AND $ 
                   (cat.w4mpro GT mag_lower) AND (abs(bgal) GT 15))
      if (feat EQ 'latent') then begin
        ra_lmc = 79.7 & dec_lmc = -68.7
        dlmc = djs_diff_angle(cat[wpsf].ra,cat[wpsf].dec,ra_lmc,dec_lmc)
        wgood = where(dlmc GT 4.5)
        wpsf = wpsf[wgood]
      endif
      ra = cat[wpsf].ra
      dec = cat[wpsf].dec
      mag = cat[wpsf].w4mpro
  endelse

end

function process_many_stars, indstart, nproc, w4=w4, feat=feat, $ 
                             subtract=subtract, nmax=nmax

  if ~keyword_set(feat) then feat='wings'
  get_star_sample, ra, dec, mag, w4=w4, feat=feat

  npsf = n_elements(ra)
  indend = (long(indstart)+nproc-1) < (npsf-1)

  for i=long(indstart), indend do begin
      print, 'processing source with i = ', i
      t0 = systime(1)
      cubestr = prepare_relevant_cutouts(ra[i], dec[i], mag[i], $ 
                     nmax=nmax, feat=feat, subtract=subtract, w4=w4)
      dt = systime(1)-t0
      print, 'processing source took: ', dt, ' seconds total'
; ----- handle case with no cutouts found corresponding to star
      if size(cubestr, /TYPE) EQ 2 then continue
      if n_elements(outstr) EQ 0 then outstr = cubestr else $
          outstr = struct_append(outstr, cubestr)
  endfor
  
  return, outstr

end

pro write_many_stars, indstart, nproc, outpath=outpath, w4=w4, feat=feat, $ 
                      subtract=subtract, nmax=nmax

  if ~keyword_set(feat) then feat='wings'
  if ~keyword_set(outpath) then outpath = '/n/panlfs/ameisner/psf'

  outstr = process_many_stars(indstart, nproc, w4=w4, feat=feat, $ 
                              subtract=subtract, nmax=nmax)

  outname = 'psfcube_'+string(indstart,format='(I04)')+'.fits'
  outname = concat_dir(outpath, outname)
  mwrfits, outstr, outname

end

