;----- this routine computes residuals of SFD vs. SSPP based on
;      Eddie's file
pro sspp_residuals, resid, specstr=specstr
;----- read file with stellar color predictions and SFD E(B-V) predictions
  specstr = mrdfits('/n/panlfs/ameisner/sspp-reddening.fits', 1)
  specstr = specstr[where(specstr.goodflag)] ; add a cut on r-i???

;----- SDSS observed colors, NOT dereddened
  catmag = specstr.catmag
  gr = catmag[1, *] - catmag[2, *]

;----- http://www.sdss.org/dr7/algorithms/sdssUBVRITransform.html
;All stars with Rc-Ic < 1.15
;        Transformation                RMS residual
;    g-r    =    1.02*(B-V)   - 0.22      0.04

  bv_obs = gr/1.02 ; + constant
;----- E(B-V) field identical to dust_getval values on specstr.ra, specstr.dec
  ebv_sfd = specstr.ebv

;----- predict B-V based on SSPP spectrum
  ssppmag = specstr.ssppmag
  bv_pred = (ssppmag[1, *] - ssppmag[2, *])/1.02 ; + same constant

;----- now find "true" E(B-V) based on observation and synthetic spectrum
  ebv_true = bv_obs - bv_pred
  resid = ebv_true - ebv_sfd ; +'tive when more reddened than SFD predicts

end

pro tile_residuals, fluxlocal, fluxsmth, tilepath=tilepath, $ 
                    akari=akari
;----- read file with stellar color predictions and SFD E(B-V) predictions
  specstr = mrdfits('/n/panlfs/ameisner/sspp-reddening.fits', 1)
  specstr = specstr[where(specstr.goodflag)] ; add a cut on r-i???

;----- use tile 15 asec tile pixels as "local" WISE/AKARI flux values
  tstr = keyword_set(akari) ? $
    mrdfits('/n/home09/ameisner/wise/pro/akaritile-index.fits', 1) : $ 
    mrdfits('/n/home08/dfink/runwise/tile/wisetile-index.fits', 1)
  spherematch, specstr.ra, specstr.dec, tstr.ra, tstr.dec, 8.766826, $ 
    mspec, mtile, maxmatch=7486

  tlist=mtile[uniq(mtile, sort(mtile))]
  if ~keyword_set(tilepath) then $ 
    tilepath = keyword_set(akari) ? '/n/itc1/ameisner/WideS' : $ 
      '/n/itc1/ameisner/tile150'

  fluxlocal = fltarr(n_elements(specstr))
  wtlocal = fltarr(n_elements(specstr))
  for i = 0, n_elements(tlist) - 1 do begin
    print, '.', format='($,A)'
    wtile = where(mtile EQ tlist[i]) ; guaranteed to be nonzero # of matches
    
    specra = specstr[mspec[wtile]].ra
    specdec = specstr[mspec[wtile]].dec

    tile = $ 
      readfits(concat_dir(tilepath, tstr[tlist[i]].fname), h, /silent)
    if ~keyword_set(akari) then tilewt = $ 
      readfits(concat_dir(tilepath, tstr[tlist[i]].fname), exten_no=2, /silent)

    extast, h, astr
    ad2xy, specra, specdec, astr, x, y

;----- good = within image boundary, eventually add cut on ORMASK also
    wgood = where((x GT -0.5) AND (x LT 2999.5) AND (y GT -0.5) $ 
                  AND (y LT 2999.5), ngood)
    if (ngood GT 0) then begin
      intflux = interpolate(tile, x[wgood], y[wgood])
      if ~keyword_set(akari) then begin
        intwt = interpolate(tilewt, x[wgood], y[wgood])
        wcov = where(intwt GT 4, nwcov) ;coverage cut
      endif else begin
        nwcov = n_elements(x[wgood])
        wcov = lindgen(nwcov)
      endelse
      if nwcov GT 0 then begin
        fluxlocal[mspec[wtile[wgood[wcov]]]] += intflux
        wtlocal[mspec[wtile[wgood[wcov]]]] += 1
      endif
    endif
  end
  print

  fluxlocal = fluxlocal/(wtlocal + (wtlocal EQ 0))

  print, 'total stars passing Eddies SSPP cuts: ', n_elements(specstr)
  print, 'total number of WISE image star positions used: ', $ 
    long(total(fluxlocal NE 0))
;----- now get 6.1' smoothed wise map
  smthmap = keyword_set(akari) ? $
readfits('/n/panlfs/ameisner/akariWideS_2048.ecl.smooth.fits', exten_no=2, /silent) : $ 
readfits('/n/panlfs/ameisner/wise12_2048.ecl.smooth.fits', exten_no=2, /silent)

;----- smooth map is in ecliptic projection...
  euler, specstr.ra, specstr.dec, lambda, beta, 3
  ang2pix_ring, 2048, (90.0-beta)/!radeg, lambda/!radeg, pix

  fluxsmth = smthmap[pix]
end

pro trendplot, akari=akari, residstr=residstr

  if ~keyword_set(residstr) then begin
    sspp_residuals, resid
    tile_residuals, fluxlocal, fluxsmth, akari=akari
  endif else begin
    print, 'READING RESIDUALS FROM FILE...'
    resid = residstr.ebvresid
    if keyword_set(akari) then begin
      fluxlocal = residstr.akariloc
      fluxsmth = residstr.akarismth
    endif else begin
      fluxlocal = residstr.wiseloc
      fluxsmth = residstr.wisesmth
    endelse
  endelse

  wcomp = where((fluxlocal NE 0) AND (fluxsmth NE 0))
  
  diff = fluxlocal[wcomp] - fluxsmth[wcomp]
  resid = resid[wcomp]

;----- try to visualize the trend in some way

  if ~keyword_set(akari) then begin
    bincenters = -10 + lindgen(201)*.1 ; in units of WISE DN
    bhalf = 0.05
  endif else begin
    bincenters = -1.2 + lindgen(241)*.01
    bhalf = 0.005
  endelse
  meddiff = fltarr(n_elements(bincenters))
  for i=0 , n_elements(bincenters) - 1 do begin
    meddiff[i] = median(resid[where((diff LT bincenters[i]+bhalf) $ 
                                   AND (diff GT bincenters[i]-bhalf))])
  endfor

  if keyword_set(akari) then begin
    djs_plot, bincenters/30., meddiff, psym=1, xtitle='LOCAL-SMOOTH (AKARI ISSA TILE UNITS)', $ 
    ytitle='SSPP-SFD (mags EBV)', xrange=[-0.04,0.04], $ 
    /xst, yrange=[-0.005, 0.01], /yst
  endif else begin
    djs_plot, bincenters, meddiff, psym=1, xtitle='LOCAL-SMOOTH (WISE DN)', $ 
      ytitle='SSPP-SFD (mags EBV)', xrange=[-10,10], /xst, /yst
  endelse

end

pro resid_file

  sspp_residuals, resid, specstr=specstr
  tile_residuals, wiseloc, wisesmth
  tile_residuals, akariloc, akarismth, /akari
  rstr = $ 
    replicate({ebvresid:0., wiseloc:0., wisesmth: 0., $ 
      akariloc: 0., akarismth: 0.}, n_elements(resid))
  rstr.ebvresid = transpose(resid)
  rstr.wiseloc = wiseloc
  rstr.wisesmth = wisesmth
  rstr.akariloc = akariloc
  rstr.akarismth = akarismth
  outstr = struct_addtags(specstr, rstr)
  outstr = struct_trimtags(outstr, except_tags='goodflag')
  mwrfits, outstr, '/n/panlfs/ameisner/ebv_summary.fits'

end

pro plotfootprint, outname

  !P.MULTI = [0, 1, 2]
  dfpsplot, outname, bits=8, /color

  indstr = mrdfits('$WISE_DATA/index-metadata-L1b.fits', 1)
  specstr = mrdfits('/n/panlfs/ameisner/sspp-reddening.fits', 1)
  spectsr = specstr[where(specstr.goodflag EQ 1)]
  euler, indstr.ra, indstr.dec, lambda, beta, 3
  euler, specstr.ra, specstr.dec, slambda, sbeta, 3
  euler, indstr.ra, indstr.dec, lgal, bgal, 1
  euler, specstr.ra, specstr.dec, slgal, sbgal, 1
  
  plot, lambda, beta, psym=3, xrange=[360,0], /xst, yrange=[-90,90], /yst, $
       xtitle=textoidl('\lambda')+' (deg)', $ 
      ytitle=textoidl('\beta')+' (deg)', $ 
      title='WISE preliminary release vs. SSPP spectroscopy footprint'
  oplot,slambda,sbeta,psym=3, color=djs_icolor('red')
  plot, lgal, bgal, psym=3, xrange=[360,0], /xst, yrange=[-90,90], /yst, $
      xtitle= 'l (deg)', $ 
      ytitle='b (deg)', $ 
      title='WISE preliminary release vs. SSPP spectroscopy footprint'
  oplot,slgal,sbgal,psym=3, color=djs_icolor('red')
  dfpsclose
     
end
