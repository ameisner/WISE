; This files is a development sandbox.

;  Runs of possible interest: 
;  2578        2861        3438        5871        6476

; do sweep files for one run
pro wise_sweep, new, type=type

  if ~ keyword_set(type) then type = 'star'

; -------- WISE band (1..4)
  t0 = systime(1)
  bandname = ['1', '2', '3']
  runstr = '3438'

; -------- Aperture radii (these are 88, 95, 98% containment at 3.5 mu)
  aper = [7., 10., 12.5,  15., 20.] ; arcsec 
  skyrad = [20.0, 28.0] ; arcsec

; -------- Read SDSS sweep file
;  obj = sweep_readobj(2334, 5, rerun=137)

  for i=1, 6 do begin 
     camcol = string(i, format='(I1)')
     infile = '~/sweeps/calibObj-00'+runstr+'-'+camcol+'-'+type+'.fits.gz'
     outfile = '~/sweeps/calibObj-00'+runstr+'-'+camcol+'-'+type+'-wise.fits'

; -------- read sweep file
     obj = mrdfits(infile, 1)
     wise_photometry, obj.ra, obj.dec, bandname, str

; -------- add new fields to obj structure
     new = struct_addtags(obj, str)
     mwrfits, new, outfile, /create
     spawn, 'gzip -v '+outfile
     print
     print
  endfor
  print, 'Elapsed time: ', systime(1)-t0

  return
end



function objmags, obj, lambda=lambda
  apnum = 1
  

; -------- 2MASS offsets from http://coolwiki.ipac.caltech.edu/index.php/Units
  ABzp = 3631.0
  dJ = -2.5*alog10(1594./ABzp)
  dH = -2.5*alog10(1024./ABzp)
  dK = -2.5*alog10(666.7/ABzp)

; -------- WISE offsets from Table 5 of
; http://wise2.ipac.caltech.edu/docs/release/prelim/expsup/sec4_3g.html#WISEZMA
;  dW = [2.683, 3.319, 5.242, 6.604]  this is no longer needed

  wisemags = 22.5-2.5*alog10(reform(obj.wise_flux[apnum, 0:2] > 0.001))

; -------- put all mags on AB system
  mags  = [22.5-2.5*alog10(obj.psfflux > .001), $
           transpose(obj.tmass_j)+dJ, $
           transpose(obj.tmass_h)+dH, $
           transpose(obj.tmass_k)+dK, $
           wisemags[0, *]+dW[0], $
           wisemags[1, *]+dW[1], $
           wisemags[2, *]+dW[2]]

  
  lambda = [.35868, .47167, .61651, .74759, .89229, $ ; SDSS from Schlafly+
            1.24829, 1.65884, 2.18977, $ ; These are from UKIRT
            3.35, 4.60, 11.56, 22.08] ; WISE Exp. Supp.
  
  return, mags
end




pro play


; -------- read file
  obj = mrdfits('calibObj-003438-1-star-wise.fits.gz',1)

  Ebv = obj.extinction[2]/2.285

; -------- SDSS mags
  g = 22.5-2.5*alog10(obj.psfflux[1]>.001)
  r = 22.5-2.5*alog10(obj.psfflux[2]>.001)
  i = 22.5-2.5*alog10(obj.psfflux[3]>.001)
  z = 22.5-2.5*alog10(obj.psfflux[4]>.001)
  gr = g-r
  ri = r-i

; -------- WISE mags
  apnum = 0
  a = 22.5-2.5*alog10(obj.wise_flux[apnum, 0]>.001)
  b = 22.5-2.5*alog10(obj.wise_flux[apnum, 1]>.001)

; -------- select stars
  w = where(a lt 17.5 and a gt 10)


  djs_plot, gr[w], ri[w], ps=3, xr=[0, 2], yr=[0, 2]

  fstar = where(gr[w] gt 0.4 and gr[w] lt 0.6 and ri[w] gt 0.1 and ri[w] lt 0.3)
  mstar = where(gr[w] gt 1.35 and gr[w] lt 1.6 and ri[w] gt 0.5 and ri[w] lt 0.9)

  wf = w[fstar]
  wm = w[mstar]
  wq = w[where((a-b)[w] gt 0.3 and obj[w].tmass_k ne 0)]

  cols = ['red', 'green', 'blue', 'cyan', 'magenta', 'yellow', 'white']

  djs_oplot, gr[wf], ri[wf], ps=3, color='blue'
  djs_oplot, gr[wm], ri[wm], ps=3, color='red'

  plot, g[w]-z[w], z[w]-a[w]+off, ps=3, xr=[-2, 5], yr=[0, 4]
  djs_oplot, (g-z)[wf], (z-a)[wf]+off, ps=3, color='blue'
  djs_oplot, (g-z)[wm], (z-a)[wm]+off, ps=3, color='red'


  plot, [3], xr=[0.3, 15], yr=[20, 10], /nodata, xtit='lambda [micron]', $
    /xlog, /xst
  mmag = objmags(obj[wm], lambda=lambda)
  for ct=0, 19 do djs_oplot, lambda, mmag[*, ct], col=cols[ct mod 7]


  plot, [3], xr=[0.3, 15], yr=[20, 10], /nodata, xtit='lambda [micron]', $
    /xlog, /xst
  fmag = objmags(obj[wf], lambda=lambda)
  for ct=0, 69 do djs_oplot, lambda, fmag[*, ct], col=cols[ct mod 7]

  oplot, [2, 10], 11+2.5*alog10([2, 10]^2), line=2
  
  plot, [3], xr=[0.3, 15], yr=[20, 10], /nodata, xtit='lambda [micron]', $
    /xlog, /xst
  qmag = objmags(obj[wq], lambda=lambda)
  for ct=0, n_elements(wq)-1 do djs_oplot, lambda, qmag[*, ct], col=cols[ct mod 7]
  



  aivar = obj.wise_flux_ivar[apnum, 0]
  bivar = obj.wise_flux_ivar[apnum, 1]

  asn = obj.wise_flux[apnum, 0]*sqrt(aivar)
  bsn = obj.wise_flux[apnum, 1]*sqrt(bivar)

  plot, a[wm], (a-b)[wm]+0.08, ps=3, yr=[-1, 1]
  aberr = sqrt(1./asn^2 + 1./bsn^2 + 0.05^2) *0.7
  djs_oplot, a[wm], aberr[wm], ps=3, color='red'
  djs_oplot, a[wm], -aberr[wm], ps=3, color='red'

  dg = 3.303*Ebv
  dz = 1.263*Ebv
  da = 0.15 *Ebv


  return
end
