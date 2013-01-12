;+
; NAME:
;  wise_correct_latent
;   
; PURPOSE:
;  Correct for latents in an exposure due to bright stars in preceding
;  exposure ~11s prior
;
; CALLING SEQUENCE:  
;  wise_correct_latent, im, h, latentmask=
;
; INPUTS:
;  im     - image with latent artifacts
;  h      - FITS header containing WCS astrometry
;
; OPTIONAL OUTPUTS:
;  latentmask - mask marking central pixels of any latents
;
; DEPENDENCIES:
;  identifying the exposure preceding im relies on binary_search.pro
;  downloaded from:
;  astro.washington.edu/docs/idl/cgi-bin/getpro/library36.html?BINARY_SEARCH
;
; REVISION HISTORY:
;   2012-Feb-18 - Written by Aaron Meisner
;----------------------------------------------------------------------
pro wise_correct_latent, im, h, latentmask=latentmask, allsky=allsky, w4=w4

  par = psf_par_struc(w4=w4, allsky=allsky, /everything)
  LPIX   = par.szlat ;sidelength of latent cutout
  LP2    = par.szlat2 ; sidelength of 2nd latent
  IMPIX  = par.impix
  MASKSIZE = 9
  kern   = shift(dist(MASKSIZE), MASKSIZE/2, MASKSIZE/2) LT (MASKSIZE/2+0.5)
  latentmask = intarr(IMPIX, IMPIX)

;----- identify exposures preceding this one by j*11.1s, for j=1,2,3,4 
  thismjd = sxpar(h, 'MJD_OBS')
  prev = previous_exposure(thismjd, allsky=allsky, w4=w4, N=4)
  if size(prev, /TYPE) EQ 2 then return
  
  COMMON LATENT, latent_cutout, latent_cutout2
  if (n_elements(latent_cutout) EQ 0) then begin
      latent_cutout_full = readfits(par.latim)
      latent_cutout = taper_cutout(latent_cutout_full, feat='latent', $ 
                                   allsky=allsky, w4=w4)
      latent_cutout2_full = readfits(par.latim2)
      latent_cutout2 = taper_cutout(latent_cutout2_full, feat='latent2', $ 
                                    allsky=allsky, w4=w4)   
  endif
 
  sz = [LPIX, LP2, MASKSIZE, MASKSIZE]
  nprev = n_elements(prev)
  for j=0, nprev-1 do begin
    hprev = lfs_fits_access(prev[j].fname, /HEADERONLY)
    wise_starlist, hprev, xstar, ystar, maglist, allsky=allsky, w4=w4
    if n_elements(maglist) EQ 0 then continue
    n = prev[j].n
    len = sz[n-1]
    wlatent = where((maglist LT (par.lmskmax)[n-1]) AND $ 
      (round(xstar) GE -len/2) AND (round(xstar) LT IMPIX+len/2) AND $ 
      (round(ystar) GE -len/2) AND (round(ystar) LT IMPIX+len/2), nlatent)

    if (n EQ 1) then print, '# of latents removed: ', nlatent
    if (nlatent EQ 0) then continue

    xstar  = xstar[wlatent]
    ystar  = ystar[wlatent]
    maglist = maglist[wlatent]

    if (n EQ 1) then $ 
        nlin_corr = latent_nonlin_corr(maglist, allsky=allsky, w4=w4)

    for i = 0, nlatent - 1 do begin    
      thisx = xstar[i]
      thisy = ystar[i]
      thismag = maglist[i]
    
      ix    = round(thisx)
      iy    = round(thisy)
     
; ----- latent model stored in .fits has been scaled to w?mag of 0.0
      if (n EQ 1) then begin
          this_cutout = sshift2d(latent_cutout*(10^(-thismag/2.5)), $ 
            [thisx-ix,thisy-iy])*nlin_corr[i]
    
          im[((ix-LPIX/2) > 0):((ix+LPIX/2) < (IMPIX-1)), $ 
            ((iy-LPIX/2) > 0):((iy+LPIX/2) < (IMPIX-1))] -= $ 
          this_cutout[((LPIX/2-ix) > 0):((IMPIX-ix+LPIX/2-1) < (LPIX-1)), $ 
            ((LPIX/2-iy) > 0):((IMPIX-iy+LPIX/2-1) < (LPIX-1))]
      endif
; ----- subtract second latent
      if (n EQ 2) && (thismag LE par.l2max) then begin $
          this_cutout = sshift2d(latent_cutout2*(10^(-thismag/2.5)), $ 
            [thisx-ix,thisy-iy])
          im[((ix-LP2/2) > 0):((ix+LP2/2) < (IMPIX-1)), $ 
            ((iy-LP2/2) > 0):((iy+LP2/2) < (IMPIX-1))] -= $ 
          this_cutout[((LP2/2-ix) > 0):((IMPIX-ix+LP2/2-1) < (LP2-1)), $ 
            ((LP2/2-iy) > 0):((IMPIX-iy+LP2/2-1) < (LP2-1))]
      endif
      if (ix GE -MASKSIZE/2) AND (ix LT IMPIX+MASKSIZE/2) AND $ 
        (iy GE -MASKSIZE/2) AND (iy LT IMPIX+MASKSIZE/2) then begin
        latentmask[((ix-MASKSIZE/2) > 0):((ix+MASKSIZE/2) < (IMPIX-1)), $ 
            ((iy-MASKSIZE/2) > 0):((iy+MASKSIZE/2) < (IMPIX-1))] OR= $ 
            kern*((par.latflag)[n-1])
      endif
    endfor
  endfor

end
