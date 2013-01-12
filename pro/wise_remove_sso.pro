;+
; NAME:
;   wise_remove_sso
;
; PURPOSE:
;   interpolate over moving objects based on their detection in individual
;   exposures, using WISE All-Sky Known Solar System Object Possible 
;   Association List
;
; CALLING SEQUENCE:
;   wise_remove_sso, im, h, w4=w4, smask=smask
;
; INPUTS:
;   im     - L1b intensity image
;   h      - corresponding fits header
;
; KEYWORDS:
;   w4     - set if L1b image is band 4
;
; OPTIONAL OUTPUTS:
;   smask  - SSO interpolation mask (0 = unchanged,1 = interpolated over)
;
; REVISION HISTORY:
;   2012-Sep-11 - Written by Aaron Meisner
;
;----------------------------------------------------------------------
pro sso_mark_latent, h, smask, w4=w4, allsky=allsky

  prev = previous_exposure(sxpar(h,'MJD_OBS'), allsky=allsky, w4=w4, N=1)
  if size(prev, /TYPE) EQ 2 then return
  hprev = lfs_fits_access(prev[0].fname, /HEADERONLY)
  wise_sso_list, hprev, x, y, mag, w4=w4
  if n_elements(mag) EQ 0 then return
  par = psf_par_struc(w4=w4, allsky=allsky, /everything)

  magmax = keyword_set(w4) ? -1. : 1.
  wlat = where(mag LT magmax, nlat)
  if (nlat EQ 0) then return

  ix = round(x[wlat])
  iy = round(y[wlat])
  
  kpix = 9
  kern = shift(dist(kpix), kpix/2, kpix/2) LT (kpix/2+0.5)

  for i=0, nlat-1 do begin
      smask[((ix[i]-kpix/2) > 0):((ix[i]+kpix/2) < (par.impix-1)), $ 
            ((iy[i]-kpix/2) > 0):((iy[i]+kpix/2) < (par.impix-1))] OR= $ 
     4*kern[((kpix/2-ix[i]) > 0):((par.impix-ix[i]+kpix/2-1) < (kpix-1)), $ 
            ((kpix/2-iy[i]) > 0):((par.impix-iy[i]+kpix/2-1) < (kpix-1))]
  endfor

end

pro wise_remove_sso, im, h, w4=w4, smask=smask, allsky=allsky, x=x, y=y

  par = psf_par_struc(w4=w4, /everything)
  IMPIX = par.impix
  w4 = keyword_set(w4)
; ----- faintest magnitude of SSO for which to mask diffraction spikes 
  spkmax = w4 ? 0.2 : 3.0
; ----- faintest magnitude of SSO for which to interpolate over ghost
  gmax = w4 ? -0.3 :  2.5
; ----- mag beyond which ghost mask will not expand for brighter SSOs
  gmin = w4 ? -4.8 : -2.0
  ghx = par.xgpix/2 ; ghost half x
  ghy = par.ygpix/2 ; ghost half y
  smask = bytarr(IMPIX, IMPIX)
  sso_mark_latent, h, smask, w4=w4, allsky=allsky

  wise_sso_list, h, x, y, mag, w4=w4
  if n_elements(mag) EQ 0 then return else nsso = n_elements(mag)

; ----- mask size determined from visual inspection
  if ~w4 then $ 
      kpix = (2*(3+ceil(((10-mag) > 0)*4.3))+1) < 83 $ 
  else $ 
      kpix = (2*(1+ceil(((7-mag) > 0)*4.1))+1) < 83
  ix = round(x)
  iy = round(y)
  ginterp = (mag LT gmax)
  for i = 0L, nsso-1 do begin
      ksize = mag[i]
      khalf = kpix[i]/2
      kern = shift(dist(kpix[i]), kpix[i]/2, kpix[i]/2) LT (kpix[i]/2+0.5)
      smask[((ix[i]-khalf) > 0):((ix[i]+khalf) < (IMPIX-1)), $ 
            ((iy[i]-khalf) > 0):((iy[i]+khalf) < (IMPIX-1))] OR= $ 
      kern[((khalf-ix[i]) > 0):((IMPIX-ix[i]+khalf-1) < (kpix[i]-1)), $ 
           ((khalf-iy[i]) > 0):((IMPIX-iy[i]+khalf-1) < (kpix[i]-1))]
      if mag[i] LE spkmax then begin
          spmsk = diff_spike_mask(mag[i], 1, w4=w4, sz=sz, dil=11)
          smask[((ix[i]-sz/2) > 0):((ix[i]+sz/2) < (IMPIX-1)), $ 
                ((iy[i]-sz/2) > 0):((iy[i]+sz/2) < (IMPIX-1))] OR= $ 
          spmsk[((sz/2-ix[i]) > 0):((IMPIX-ix[i]+sz/2-1) < (sz-1)), $ 
                ((sz/2-iy[i]) > 0):((IMPIX-iy[i]+sz/2-1) < (sz-1))]
      endif
      if ginterp[i] then begin
; ----- determine if ghost in this image using wise_l1b_cutout
          incl = wise_l1b_cutout(_, ix[i], iy[i]-par.ygoffs, $ 
                                 par.xgpix, par.ygpix, w4=w4, /bool, /silent)
          if ~incl then continue
          dist_ellipse, dell, [par.xgpix, par.ygpix], par.xgpix/2, $ 
              par.ygpix/2+2, float(par.ygpix)/par.xgpix,0.
; ----- these sizes are based on inspection of W3 images, W4 not yet considered
          radg = (15.+(10./(gmax-gmin))*(gmax-mag[i])) < 25
          gmsk = 2*(dell LT radg)
          smask[((ix[i]-ghx) > 0):((ix[i]+ghx) < (IMPIX-1)), $ 
      ((iy[i]-par.ygoffs-ghy) > 0):((iy[i]-par.ygoffs+ghy) < (IMPIX-1))] OR= $ 
          gmsk[((ghx-ix[i]) > 0):((IMPIX-ix[i]+ghx-1) < (par.xgpix-1)), $ 
 ((ghy+par.ygoffs-iy[i]) > 0):((IMPIX+par.ygoffs-iy[i]+ghy-1) < (par.ygpix-1))]
      endif
  endfor

  intx = djs_maskinterp(im, (smask NE 0), iaxis=0, /const)
  inty = djs_maskinterp(im, (smask NE 0), iaxis=1, /const)
  im = (intx+inty)/2

end

pro test_xy_coords, indstart, nproc, w4=w4

; curious to know whether entries in SSO catalog are always at coordinates
; strictly within boundaries of relevant L1b exposure
; conclusions for W3, 4.3777182 < x < 1010.2629, 3.6281532  < y < 1009.1289
; conclusions for W4, 0.9905761 < x < 504.24374, 3.8662809 < y < 502.27559

  init_sso_catalog, w4=w4
  COMMON SSO, ra, dec, w3mag, w4mag, id, scan_id, frame_num
  ind_unique = uniq(id)

  scan_unique = scan_id[ind_unique]
  frame_unique = frame_num[ind_unique]

  id_unique = scan_unique + string(frame_unique, format='(I03)')

  indstr = wise_index_metadata([0.,0.], /allsky, /nosort, w4=w4)
  par = psf_par_struc(/allsky, w4=w4, /everything)
  id_all = indstr.scan_id + string(indstr.frame_num, format='(I03)')
  matchlist, id_all, id_unique, ind_all, _
  indstr = indstr[ind_all]

  indend = (long(indstart)+nproc-1) < (n_elements(indstr)-1)
  xmax = -1e6
  xmin = 1e6
  ymax = -1e6
  ymin = 1e6
  for i=long(indstart), indend do begin
      print, i
      h=headfits(indstr[i].fname, /silent)
      wise_remove_sso, fltarr(par.impix,par.impix), h, w4=w4, x=x, y=y
      if n_elements(x) EQ 0 then continue
      xmax >= max(x)
      xmin <= min(x)
      ymax >= max(y)
      ymin <= min(y)
      delvarx, x, y
  endfor

  print, ymax, ymin
  print, xmax, xmin
  outstr = {xmin: xmin, xmax: xmax, ymin: ymin, ymax: ymax}
  mwrfits, outstr, 's'+strtrim(string(indstart),1)+'.fits'
end
