;+
; NAME:
;   fit polynomial PSF model based on cutouts
;
; PURPOSE:
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
;   lots of other experimental PSF model fitting code included as well
;
; REVISION HISTORY:
;   2011-Oct-16 - Aaron Meisner
;----------------------------------------------------------------------
function gather_cutouts, fpath=fpath, mag_u=mag_u

  if ~keyword_set(fpath) then fpath='/n/panlfs/ameisner/psf/logs/w4'
  spawn, 'ls '+fpath+'/*.fits', flist

  nfile = n_elements(flist)
  if (nfile EQ 0) then begin
      print, 'no cutout files found!!!'
      return, -1
  endif

  outstr = mrdfits(flist[0], 1)
  for i=1,nfile-1 do begin
      print, i
      addstr = mrdfits(flist[i],1)
      if keyword_set(mag_u) then begin
          wkeep = where(addstr.mag LT mag_u)
          addstr = addstr[wkeep]
      endif
      outstr = struct_append(outstr, addstr)
  endfor

  return, outstr

end

function linfit_one_pixel, xvals, yvals, imvals, wtvals, w4=w4, order=order, $
                           iwt=iwt

  if ~keyword_set(order) then order = 1

  npar_list = [3, 6, 10]

  wgood = where(wtvals NE 0, nwgood)
  print, 'number of cutouts which include this pixel location : ', nwgood
  if nwgood EQ 0 then begin
      print, 'no valid cutouts for this pixel??'
      return, -1
  endif

  par = psf_par_struc(w4=w4)
  crpix = par.crpix

; linear gradients in x, y
  nparam = npar_list[order-1]
  dx = xvals[wgood]-crpix
  dy = yvals[wgood]-crpix
  AA = fltarr(nparam, nwgood)
  AA[0,*] = 1
  AA[1,*] = dx
  AA[2,*] = dy
  if (order GT 1) then begin
      AA[3,*] = dx*dy
      AA[4,*] = dx^2
      AA[5,*] = dy^2
  endif
  if (order EQ 3) then begin
      AA[6,*] = (dx^2)*dy
      AA[7,*] = (dy^2)*dx
      AA[8,*] = dx^3
      AA[9,*] = dy^3
  endif
  YY = imvals[wgood]
  WW = wtvals[wgood]

  hogg_iter_linfit, AA, YY, WW, coeff

  iwt = bytarr(n_elements(xvals))
  if (total(WW) GT 0) then iwt[wgood[where(WW NE 0)]] = 1
  return, coeff

end

function cutouts_per_star, dec

; determine, for each cutout, the number of times the relevant star
; appears in the data cube

  dec_unique = dec[uniq(dec, sort(dec))]
  nunique = n_elements(dec_unique)

; there is doubtless a more efficient/faster/more clever way to do this
; with uniq and matchlist, but i don't care at this point

  nexp = lonarr(n_elements(dec))
  for i=0L, nunique-1 do begin
      w = where(dec EQ dec_unique[i], nw)
      nexp[w] = nw
  endfor

  return, nexp

end

pro remove_nan, cutstr

  w = where(~finite(cutstr.cutout), nw)
  if (nw EQ 0) then return
  sz = size(cutstr[0].cutout, /DIM)
  szx = long(sz[0])
  szy = sz[1]

  indbad = w / (szx*szy)
  indbad = indbad[uniq(indbad, sort(indbad))]
  good = bytarr(n_elements(cutstr))+1
  good[indbad] = 0
  wgood = where(good)
  cutstr = cutstr[wgood]

end

pro psf_model_fit, coeff, w4=w4, fpath=fpath, feat=feat, order=order, $ 
                   iwtcube=iwtcube

  if ~keyword_set(feat) then feat = 'wings'
  ;if ~keyword_set(order) then order = 1 ; linear fit
  par = psf_par_struc(w4=w4, feat=feat)
  if ~keyword_set(order) then order = par.order
  print, 'order = ', order

  if (order LT 1) or (order GT 3) then begin
      print, 'this code is only intended for polynomials of order 1, 2, 3'
  endif
  npar_list = [3, 6, 10]
  nparam = npar_list[order-1]

  if (~keyword_set(w4) AND (feat EQ 'latent')) then mag_u = 1.
  cutstr = gather_cutouts(fpath=fpath, mag_u=mag_u)
; try this temporarily !!!!
  wgood = where(cutstr.bg NE 0) ; bg = 0 means djs_photsky failed
  cutstr = cutstr[wgood]
  remove_nan, cutstr
  ncut = n_elements(cutstr)

  par = psf_par_struc(w4=w4, feat=feat)
  szx = par.szx
  szy = par.szy 
  npix = long(szx)*szy

; --- W3 second latent is very close to linear, so don't make correction for 
;     feat = 'latent2' case
  corr = (feat EQ 'latent') ? latent_nonlin_corr(cutstr.mag, w4=w4) : $ 
      fltarr(ncut)+1.
  corr = 1./corr ; matter of semantics/defintions here...
  nexp = cutouts_per_star(cutstr.dec)
  cutcube = cutstr.cutout
  unccube = cutstr.unc
  mskcube = cutstr.wt
; ---- rescale arbitrarily to mag = 0
;  cutcube /= 10^(-cutstr.mag[lindgen(npix*ncut) / npix]/2.5)
  for i=0L, ncut-1 do begin
      cutcube[*,*,i] *= 10^(cutstr[i].mag/2.5)*corr[i]
      unccube[*,*,i] *= 10^(cutstr[i].mag/2.5)*corr[i] ; propagation of errors
  endfor
  wtcube = (mskcube NE 0)*(1/(unccube^2 + (unccube EQ 0)))

  for i=0L, ncut-1 do begin
       wtcube[*,*,i] /= nexp[i]
  endfor

; ---- generalize dimensions of coeff later depending on order of
;      polynomial fit for different features
  coeff = fltarr(szx, szy, nparam)
  iwtcube = bytarr(szx, szy, ncut)
  for x=0, szx-1 do begin
      for y=0, szy-1 do begin
        print, x, y 
        c = linfit_one_pixel(cutstr.x, cutstr.y, cutcube[x,y,*], $ 
                               wtcube[x,y,*], w4=w4, order=order, iwt=iwt)
        coeff[x,y,*] = c
        iwtcube[x,y,*] = iwt
      endfor
  endfor

end

; --- remaining code is relatively experimental ...

function single_latent_comparison, l0, l1

; do one pairwise comparison between two different latents, with a
; linear fit of latent1 = m*latent0+b

  sz = (size(l0,/DIM))[0]
  xbox = (lindgen(sz,sz) MOD sz)-sz/2
  ybox = (lindgen(sz,sz) / sz)-sz/2
  dist = sqrt(xbox^2+ybox^2)

  wann = where((dist GT 8) AND (dist LT 23), nwann) ; w3 first latent
;  wann = where((dist GT 5) AND (dist LT 20), nwann) ; w4 both, w3 2nd latent
  
  nparam = 2
  AA = fltarr(nparam, nwann)
  AA[0,*] = 1
  AA[1,*] = l0[wann]
  YY = l1[wann]
  WW = replicate(1., nwann)
  
  hogg_iter_linfit, AA, YY, WW, coeff

  return, coeff[1]

end

function latent_median_cube, fpath=fpath, w4=w4, latstr=latstr

  if ~keyword_set(fpath) then $ 
      fpath = '/n/panlfs/ameisner/psf/logs/w4/latent'
  if ~keyword_set(latstr) then latstr = gather_cutouts(fpath=fpath)

  dec = latstr.dec

  dec_unique = dec[uniq(dec, sort(dec))]
  nunique = n_elements(dec_unique)

  full = bytarr(n_elements(latstr))
  for i=0L,n_elements(latstr)-1 do begin
      full[i] = total(latstr[i].wt EQ 0) EQ 0
  endfor

  latcube = latstr.cutout
  for i=0L, nunique-1 do begin

      w = where((latstr.dec EQ dec_unique[i]) AND (full), nw)
      if nw EQ 0 then begin
          print, 'no full cutouts for this source ???'
          continue
      endif
      thismed  = (nw GT 1) ? median(latcube[*,*,w],dimension=3) : $ 
                             latcube[*,*,w]
      addstr = {ra: latstr[w[0]].ra, dec: latstr[w[0]].dec, $
                mag: latstr[w[0]].mag, cutout: thismed, ncut: nw}
      if ~keyword_set(medstr) then begin
          medstr = addstr
      endif else begin
          medstr = struct_append(medstr, addstr)
      endelse

  endfor

  return, medstr

end

function medcube_conserve_mem, fpath=fpath

  spawn, 'ls '+fpath+'/*.fits', flist
  
  for i=0, n_elements(flist)-1 do begin
      print, i
      latstr = mrdfits(flist[i], 1)
      addstr = latent_median_cube(w4=w4, latstr=latstr)
      if ~keyword_set(medstr) then begin
          medstr = addstr
      endif else begin
          medstr = struct_append(medstr, addstr)
      endelse
  endfor

  return, medstr
end

pro latent_nonlin_bin, medcube, scalefac, bincenters, fpath=fpath, w4=w4, $ 
                       feat=feat, medstr=medstr

; fit latent nonlinearity by binning in magnitudes, taking median
; and then fitting pairwise rescaling factors between bins

  if ~keyword_set(feat) then feat = 'latent'
  if ~keyword_set(medstr) then medstr = latent_median_cube(w4=w4, fpath=fpath)
  par = psf_par_struc(w4=w4, /allsky, feat=feat)

  if keyword_set(w4) then begin
; ----- brightest star in sample for w4 : -4.77
      if (feat EQ 'latent') then begin
          nbin = 30
          bincenters = 1.1-0.2*lindgen(nbin)
          bhalf = 0.1
      endif else begin
; ----- 2nd latent case
          nbin = 19
          bincenters = -1.1-0.2*lindgen(nbin)
          bhalf = 0.1
      endelse
  endif else begin
      if (feat EQ 'latent') then begin
          nbin = 35
          bincenters = -2.9 + 0.2*lindgen(nbin)
          bhalf = 0.1
      endif else begin
; ----- if you're here then you're fitting 2nd latent
        nbin = 14
        bincenters = -0.55 - 0.1*lindgen(nbin)
        bhalf = 0.05
      endelse
  endelse

  
  latcube = medstr.cutout
  medcube = fltarr(par.szx, par.szy, nbin)
  for i=0,nbin-1 do begin
      print, i, bincenters[i]-bhalf, bincenters[i]+bhalf
      wmag = where((medstr.mag LE (bincenters[i]+bhalf)) AND $ 
                   (medstr.mag GT (bincenters[i]-bhalf)))
      medcube[*,*,i] = median(latcube[*,*,wmag],dimension=3)
  endfor

; compare everything to brightest bin
  scalefac = fltarr(nbin-1)
  l0 = medcube[*,*,nbin-1]
  for i=0,nbin-2 do begin
      l1 = medcube[*,*,i]
      fac = single_latent_comparison(l0, l1)
      print, bincenters[i], fac
      scalefac[i] = fac
  endfor
  scalefac = [scalefac,1.]

end

pro nonlin_quad, magvals, correction, coeff

; fit quadratic to nonlinearity correction for w4 latent, as a
; function of MAGNITUDE w4mpro (not flux!!)

;  latent_nonlin_bin, medcube, scalefac, magvals
  scalefac = $ 
readfits('/n/panlfs/ameisner/psf/logs/w4/latent/faint/results/scalefac.fits')
  magvals = $ 
readfits('/n/panlfs/ameisner/psf/logs/w4/latent/faint/results/bincenters.fits')
  flux=10^(-(magvals+4.7)/2.5)

  correction = scalefac/flux
  nmag = n_elements(magvals)

  nparam = 3
  AA = fltarr(nparam, nmag)
  AA[0,*] = 1
  AA[1,*] = magvals
  AA[2,*] = magvals^2
  YY = correction
  WW = replicate(1., nmag)

  hogg_iter_linfit, AA, YY, WW, coeff

end

pro rescale_mag0, latcube, medstr, med

; use information from nonlin_quad fit to rescale each latent to w4=0 based on
; its w4 magnitude (including nonlinearity), then median filter
; to arrive at w4 = 0 latent image

  medstr = latent_median_cube(/w4)
  par = psf_par_struc(/w4, feat='latent')

  coeff = [2.1096788,-0.016110173,-0.054974536]

  magmin = -4.8
  magmax = 0.

; effective magnitude for computing nonlin correction
  meff = (medstr.mag > magmin) < magmax
  nlin_corr = (coeff[0]+meff*coeff[1]+(meff^2)*coeff[2])/coeff[0]
  rescale = 1./nlin_corr
  nlat = n_elements(medstr)
  latcube = medstr.cutout
  for i=0, nlat-1 do begin
      latcube[*,*,i] = latcube[*,*,i]*rescale[i]*(10^(medstr[i].mag/2.5))
  endfor

  wtcube = bytarr(par.szx,par.szy,nlat)
  for x=0,szx-1 do begin
      for y=0,szy-1 do begin
         djs_iterstat, latcube[x,y,*],mask=pixmask
         wtcube[x,y,*] = pixmask
      endfor
  endfor

  nrej = total(total(~wtcube,1),1)
  wgood = where(nrej LT 0.2*par.szx*par.szy) ; 0.2 could be better tuned???

  med = median(latcube[*,*,wgood], dimension=3)

end

function visualize_ghost, xstar, ystar, coeff, w4=w4

  par = psf_par_struc(w4=w4)
  dx = xstar-par.crpix
  dy = ystar-par.crpix
  
  cutout = coeff[*, *, 0] + $ 
           coeff[*, *, 1]*dx + $
           coeff[*, *, 2]*dy + $
           coeff[*, *, 3]*dx*dy + $
           coeff[*, *, 4]*dx^2 + $ 
           coeff[*, *, 5]*dy^2 + $ 
           coeff[*, *, 6]*(dx^2)*dy + $
           coeff[*, *, 7]*(dy^2)*dx + $
           coeff[*, *, 8]*(dx^3) + $ 
           coeff[*, *, 9]*(dy^3)
   return, cutout

end

function shifted_ghost, wfit, xyshift, template=template

  thiscutout = template ; don't modify keyword parameter
  ;print,xyshift
  thiscutout = sshift2d(thiscutout, xyshift)
  
  return, thiscutout[wfit]

end

function fit_ghost_centroid, cutout, template


  ;print,size(cutout,/DIM)
  ;print, n_elements(wann)

  xyshift = [0.1d, 0.1d]; starting value for shift
 

  wfit = lindgen(n_elements(cutout)) ; just a placeholder, does nothing
  functargs = {template : template}

  err = fltarr(n_elements(cutout))+1
  result = $ 
      MPFITFUN('shifted_ghost', wfit, cutout[wfit], err, $ 
               xyshift, functargs=functargs, /quiet)

  return, result

end

function scaled_ghost, wfit, scalefac, template=template

  thiscutout = template ; don't modify keyword parameter
  ;print,xyshift
  thiscutout = thiscutout*scalefac
  
  return, thiscutout[wfit]

end

function fit_ghost_amplitude, cutout, template

  scalefac = 1.1 ; starting value for scalefac
 
  wfit = lindgen(n_elements(cutout)) ; just a placeholder, does nothing
  functargs = {template : template}

  err = fltarr(n_elements(cutout))+1
  result = $ 
      MPFITFUN('scaled_ghost', wfit, cutout[wfit], err, $ 
               scalefac, functargs=functargs, /quiet)

  return, result

end

pro measure_ghost_centroid, coeff, w4=w4, xshift, yshift, scalefac

  binfac = 10
  par = psf_par_struc(w4=w4, /everything)
  xmin = 0
  xmax = par.impix-1
  ymin = par.ygoffs
  ymax = par.impix+par.ygoffs-1
  xntrial = floor((xmax-par.crpix)/binfac)+floor((par.crpix-xmin)/binfac)+1
  yntrial = floor((xmax-par.crpix)/binfac)+floor((par.crpix-xmin)/binfac)+1
  xtrial = binfac*(lindgen(xntrial)-xntrial/2)+par.crpix
  ytrial = binfac*(lindgen(yntrial)-yntrial/2)+par.crpix+par.ygoffs

; perform fits for 0 <= x <= IMPIX-1
; and ygoffs <= y <= IMPIX-1+ygoffs

  xshift = fltarr(xntrial, yntrial)
  yshift = fltarr(xntrial, yntrial)
  scalefac = fltarr(xntrial, yntrial)
  template = visualize_ghost(par.crpix, par.crpix, coeff, w4=w4)
  for ix=0, xntrial-1 do begin
      for iy=0, yntrial-1 do begin
          print, ix, iy
          ghost = visualize_ghost(xtrial[ix], ytrial[iy], coeff, w4=w4)
          xyshift = fit_ghost_centroid(ghost,template)
          shifted_template = sshift2d(template, xyshift)
          this_scalefac = fit_ghost_amplitude(ghost, shifted_template)
          print, this_scalefac
          xshift[ix, iy] = xyshift[0]
          yshift[ix, iy] = xyshift[1]
          scalefac[ix, iy] = this_scalefac
      endfor
  endfor

end

pro ghost_centroid_poly, shiftmap, xtrial, ytrial, order=order,coeff, $ 
                         model=model

; fit polynomials to output of measure_ghost_centroid
; fit linear function to x centroid, second order to y centroid

  if ~keyword_set(order) then order = 1
  npar_list = [3, 6, 10]
  nparam = npar_list[order-1]

  sz = (size(shiftmap, /DIM))[0] ; square

  ind_x = lindgen(sz,sz) MOD sz
  ind_y = lindgen(sz,sz) / sz

  xvals = xtrial[ind_x]
  yvals = ytrial[ind_y]

  dx = xvals - 507.5
  dy = yvals - 507.5
  ;dx = (lindgen(sz,sz) MOD sz)-sz/2
  ;dy = (lindgen(sz,sz) / sz)-sz/2
  
  dx = reform(dx,sz*sz)
  dy = reform(dy,sz*sz)

  AA = fltarr(nparam, sz*sz)
  AA[0,*] = 1
  AA[1,*] = dx
  AA[2,*] = dy
  if (order GT 1) then begin
      AA[3,*] = dx*dy
      AA[4,*] = dx^2
      AA[5,*] = dy^2
  endif
  if (order EQ 3) then begin
      AA[6,*] = (dx^2)*dy
      AA[7,*] = (dy^2)*dx
      AA[8,*] = dx^3
      AA[9,*] = dy^3
  endif
  YY = reform(shiftmap, sz*sz)
  WW = fltarr(sz*sz)+1.

  hogg_iter_linfit, AA, YY, WW, coeff

  model = coeff[0]+coeff[1]*(xvals-507.5)+coeff[2]*(yvals-507.5)

end
