;+
; NAME:
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
;   
; REVISION HISTORY:
;   2011-Nov-25 - Written by Douglas Finkbeiner, CfA
;   2012-Feb-13 - Added cleanpath keyword to wise_name
;                 call and to wise_mosaic1b definition  - Aaron Meisner
;----------------------------------------------------------------------


function binary_search_defunct, list, val
  
  n = n_elements(list) 
  step = n/2
  
  ind = n/2

  for bar=1, (alog(n)/alog(2))+1 do begin 
     step = (step/2) > 1
     ind += (list[ind] LT val ? step : -step)
     ind = (ind < (n-1)) > 0
;     print, step, ind, list[ind]
  endfor

  return, ind
end

pro wise_mosaic1b, im, raw, wt, minim, maxim, amask, omask, nside=nside, nx=nx, sz=sz, lgal=lgal, bgal=bgal, dust=dust, hdr=hdr, astr=astr, imref=imref, wtref=wtref, sigthresh=sigthresh, goodlist=goodlist, cleanpath=cleanpath, wisezp=wisezp, allsky=allsky, warp=warp, w4=w4

  t0 = systime(1)

  if ~ keyword_set(cleanpath) then cleanpath = '$WISE_DATA/clean'
  if ~ keyword_set(sigthresh) then sigthresh = 25
  wisezp = keyword_set(wisezp)

  thishost = getenv('HOSTNAME')   ; get $HOSTNAME from shell
  runlocal = (thishost EQ 'pan') OR (thishost EQ 'nebel.rc.fas.harvard.edu') $ 
                                 OR (thishost EQ 'wise.rc.fas.harvard.edu')
  
  maxrad = 0.6  ; L1b

; -------- initialize mosaic projection WCS header (if not healpix)
  if keyword_set(nside) then begin 
     healpix = 1
     npix = 12LL * nside * nside
     print, 'HEALPIX nside:', nside, '   npix:', npix
  endif else begin 
     healpix = 0
     if keyword_set(astr) then begin 
        print, 'using the astrometry structure you passed'
        sz = astr.naxis[0]
        galcoords = strmid(astr.ctype[0], 0, 4) EQ 'GLON'
        if galcoords then begin 
           lgal = astr.crval[0]
           bgal = astr.crval[1]
        endif else begin
           euler, astr.crval[0], astr.crval[1], lgal, bgal, 1
        endelse

; -------- Make WCS header corresponding to astr        
        mkhdr, hdr, fltarr(sz, sz, /nozero), /extend
        putast, hdr, astr

     endif else begin 
        cd = 2.75 * (1000/nx) / 3600.d
        dummy = fltarr(sz, sz, /nozero)
        hdr = wise_mosaic_gal_header(dummy, cd, lgal=lgal, bgal=bgal)
        dummy = 0
        extast, hdr, astr
     endelse
     galcoords = strmid(astr.ctype[0], 0, 4) EQ 'GLON'
     eclcoords = strmid(astr.ctype[0], 0, 4) EQ 'ELON'
  endelse

; -------- compute (RA, dec) for each pixel in mosaic
  print, 'Preparing mosaic coordinates...'
  if healpix then begin 
     healgen_lb, nside, lbox, bbox
     abox = fltarr(npix)   ; try this with float -- should really use double!!
     dbox = fltarr(npix)
     nslice = 16
     for k=0L, nslice-1 do begin     ; loop for memory efficiency
        print, '.', format='($,A)'
        nrow = npix/nslice  ; rows per slice
        euler, lbox[k*nrow:(k+1)*nrow-1], bbox[k*nrow:(k+1)*nrow-1], atemp, dtemp, 2
        abox[k*nrow:(k+1)*nrow-1] = temporary(atemp)
        dbox[k*nrow:(k+1)*nrow-1] = temporary(dtemp)
     endfor
     delvarx, lbox, bbox
  endif else begin 
     nslice = 16
     nrow = ceil(sz/float(nslice)) ; rows per slice
     xbox = (lindgen(sz, nrow) mod sz)+0.d  ; these are pixel *centers*
     ybox = (lindgen(sz, sz)  /  sz)+0.d
     abox = dblarr(sz, sz)
     dbox = dblarr(sz, sz)
     for k=0L, nslice-1 do begin     ; loop for memory efficiency
        print, '.', format='($,A)'
        y0 = k*nrow
        y1 = ((k+1)*nrow-1) < (sz-1)
        xy2ad, xbox[*, 0:y1-y0], ybox[*, y0:y1], astr, atemp, dtemp
        if galcoords then euler, atemp, dtemp, select=2 ; atemp,dtemp were Galactic l,b
        if eclcoords then euler, atemp, dtemp, select=4 ; atemp,dtemp were ecliptic lam,beta
        abox[*, y0:y1] = atemp
        dbox[*, y0:y1] = dtemp
     endfor
     print
     delvarx, xbox, ybox, atemp, dtemp
     asmall = rebin(abox, sz/2, sz/2, /sample)  ; this fails at l=0!!!
     dsmall = rebin(dbox, sz/2, sz/2, /sample)
     euler, asmall, dsmall, lsmall, bsmall, 1
     dust = dust_getval(lsmall, bsmall, /noloop, /interp)
     delvarx, asmall, dsmall, lsmall, bsmall
  endelse

  print, 'Sorting declinations...'
  sind = sort(dbox)

  alist = abox[sind] 
  dlist = dbox[sind]

  delvarx, abox, dbox

  print, 'Computing unit vectors'
  uv = fltarr(n_elements(alist), 3)
  nslice = 16
  for k=0L, nslice-1 do begin 
     print, '.', format='($,A)'
     nind = n_elements(alist)/nslice
     ind0 = k*nind
     ind1 = (k+1)*nind-1
     uv[ind0:ind1, *] = float(ll2uv([[alist[ind0:ind1]], [dlist[ind0:ind1]]]))
  endfor
  print 

; -------- allocate some big arrays (this hurts!)
  fluxlist = fltarr(n_elements(dlist))
  rawlist  = fltarr(n_elements(dlist))
  wtlist   = fltarr(n_elements(dlist))
  minlist  = fltarr(n_elements(dlist))+1E20
  maxlist  = fltarr(n_elements(dlist))-1E20
  andlist  = intarr(n_elements(dlist))+fix(32511)
  orlist   = intarr(n_elements(dlist))

; -------- read metadata table
  indstr = wise_index_metadata([lgal, bgal], nimage=nimage, allsky=allsky, $ 
                               w4=w4)
  wgood  = where(goodlist)
  indstr = indstr[wgood]
  flist  = indstr.fname

; -------- loop over files and add them to mosaic. 
  npad = 0  ; for the 512x512 clean files!

; -------- make (x,y) coords for WISE image
  par = psf_par_struc(w4=w4, allsky=allsky)
  nc = par.pclean
  xbox = (lindgen(nc, nc) mod nc)+0.d ; these are pixel *centers*
  ybox = (lindgen(nc, nc)  /  nc)+0.d

  for i=0L, n_elements(flist) -1 do begin 

     ;if n_elements(clean) EQ 1 then print, 'missing extension, skipped...'
;    nc = (size(clean, /dimen))[0]
     crval = [indstr[i].ra, indstr[i].dec]
; -------- check if crval is even in the target projection
     if healpix then goodcen = 1 else $ 
     begin 
        if galcoords then begin ; if GAL projection, convert to (l,b) first
           euler, crval[0], crval[1], thisl, thisb, 1
           ad2xy, thisl, thisb, astr, xcen, ycen
        endif else begin 
           if eclcoords then begin 
              euler, crval[0], crval[1], thisl, thisb, 3 ; ecliptic
              ad2xy, thisl, thisb, astr, xcen, ycen
           endif else begin 
              ad2xy, crval[0], crval[1], astr, xcen, ycen
           endelse
        endelse
        goodcen = (xcen GE 0) AND (xcen LT sz) AND (ycen GE 0) AND (ycen LT sz)
     endelse
        
; -------- Get indices for relevant dec range     
     if goodcen then begin 
        print, i, flist[i], format='(I8,"  ",A)'
        binary_search, dlist, crval[1]-0.90, ind0
        binary_search, dlist, crval[1]+0.90, ind1
        if (ind0 EQ -1) then ind0 = 0L
        if (ind1 EQ -1) then ind1 = n_elements(dlist)-1

        uvcen = float(ll2uv(transpose(crval)))
        near = where(transpose(uvcen)##uv[ind0:ind1, *] GT cos(maxrad*!dpi/180.d), nnear)+ind0
     endif else nnear = 0

     if nnear GT 0 then begin 

;        clean = readfits(flist[i], hclean, /silent, exten=0)

        cleanname = wise_name('clean', flist[i], cleanpath=cleanpath)
        if file_test(cleanname) eq 0 then message, 'cannot find file'
        clean = lfs_fits_access(cleanname, hclean, /silent, exten=0)
        dirty = lfs_fits_access(cleanname, /silent, exten=1)
        bmask = lfs_fits_access(cleanname, /silent, exten=2)
        if n_elements(clean) EQ 1 or n_elements(dirty) EQ 1 or n_elements(bmask) EQ 1 then begin
           print, 'missing extension, skipping...'
           delvarx, clean, dirty, bmask
           continue
        endif
        if keyword_set(warp) then begin
            wise_exposure_warp, clean, indstr[i].grad, w4=w4
            wise_exposure_warp, dirty, indstr[i].grad, w4=w4
        endif
; we cannot trust the zero point, so subtract it for now. 
        extast, hclean, aout
        nc = sxpar(hclean, 'NAXIS1')
        if ~wisezp then begin
           medclean = median(clean)
           clean -= medclean
           dirty -= medclean
           if total(wtlist) EQ 0 then begin 
             print, 'subtracted ', medclean, ' DN'
             sxaddpar, hdr, 'OFFSET', medclean
           endif
        endif

        wnan = where(finite(dirty) EQ 0, nwnan)
        if nwnan GT 0 then dirty[wnan] = median(dirty)   ;?????
        wnan = where(finite(clean) EQ 0, nwnan)
        if nwnan GT 0 then clean[wnan] = median(clean)   ;?????
     
        ad2xy, alist[near], dlist[near], aout, xx, yy  ; (xx, yy) are pixel centers

        w = where(xx gt npad and xx lt (nc-npad-1) and yy gt npad and yy lt (nc-npad-1), nw)
        if nwnan GT n_elements(clean)/2 then begin 
           message, /info, 'too many NaNs'
           nw = 0
        endif

        if nw gt 0 then begin 
           if total(~finite(clean)) GT 0 then stop

; -------- if we have made it this far, read the dirty image also
;           dirty = readfits(flist[i], hdr2, /silent, exten=1)

           nearw = near[w]
           intflux = interpolate(clean, xx[w], yy[w])
           intdirt = interpolate(dirty, xx[w], yy[w])
           intmask = bmask[round(xx[w]), round(yy[w])]

           wwt = where(wtlist[nearw] GE 0.5, nwwt)
           if (~wisezp) && (nwwt GT 100) then begin 
              dflux = median(intflux[wwt] - (fluxlist[nearw[wwt]]/wtlist[nearw[wwt]]), /even)
           endif else dflux = 0

; -------- compare to reference image; reject if too different
;           imref[sind[nearw]] is like dflux
           if keyword_set(imref) then begin 
              npsmall = nc/4
              xy2ad, xbox, ybox, aout, cra, cdec
              euler, cra, cdec, cl, cb, 1
              ad2xy, cl, cb, astr, imx, imy
              intref = interpolate(imref, imx, imy)
              intwt  = interpolate(wtref, imx, imy)
              wset, 1
              bclean   = bytscl(rebin(clean, npsmall, npsmall), min=-400, max=700)
              thiswt   = (intwt-1) > 0
              compmask = thiswt GE 0.5
              wwt = where(compmask, nwwt)
              if nwwt GT 1000 then begin 
                 dclean = median(clean[wwt]-intref[wwt], /even)
                 thisref = (thiswt NE 0)*(intref*(intwt > 0)-clean+dclean)/(thiswt + (thiswt eq 0))
                 bref    = bytscl(rebin(thisref, npsmall, npsmall), min=-400, max=700)
                 diffim  = clean-thisref
                 diffim -= median(diffim[wwt])
                 bdiff   = bytscl(rebin(diffim, npsmall, npsmall), min=-100, max=100)
                 barr    = [bclean, bref, bdiff]
                 tv, barr
                 wset, 0
                 littlediff = rebin(diffim, npsmall, npsmall)
                 littlemask = rebin(float(compmask), npsmall, npsmall)
                 lw = where(littlemask GE 0.5)
                 mysig = djsig(littlediff[lw], sigrej=10)
                 print, 'sigma:', mysig, '         dflux:', dflux
              endif else print, 'not enough overlap'
           endif else mysig = 0

           if mysig LT sigthresh then begin 
              intflux_dflux = intflux-dflux
              fluxlist[nearw] += intflux_dflux
              rawlist[nearw]  += (intdirt-dflux)
              wtlist[nearw]++
              maxlist[nearw] >= intflux_dflux   ; take the greater of the two
              minlist[nearw] <= intflux_dflux   ; take the lesser of the two
              andlist[nearw] AND= intmask
              orlist[nearw]  OR=  intmask
           endif else begin
              print, 'BAD!!!'
           endelse
           
; just test this for now
           if runlocal then begin 
              if (i mod 100) eq 0 then begin
                 flux = fltarr(sz, sz)
                 wt   = fltarr(sz, sz)
                 flux[sind] = fluxlist
                 wt[sind]   = wtlist
                 im  = flux/(wt + (wt eq 0))
                 offs = keyword_set(wisezp)*median(im[where(im NE 0)]) 
              
;             tv, bytscl(rebin(im, 500, 500), min=-200, max=900)
                 if (sz MOD 500) EQ 0 then $ 
                     tv, bytscl(rebin(im, 500, 500), min=-100+offs, max=300+offs) $ 
                 else $ 
                     tv, bytscl(congrid(im, 500, 500, /interp), min=-100+offs, max=300+offs)
              endif
           endif
        endif
     endif
  endfor

  print, 'mapping back to 2D'

  if healpix then begin
     flux = fltarr(npix)
     fraw = fltarr(npix)
     wt   = fltarr(npix)
     minim= fltarr(npix)
     maxim= fltarr(npix)
     amask= intarr(npix)
     omask= intarr(npix)
  endif else begin
     flux = fltarr(sz, sz)
     fraw = fltarr(sz, sz)
     wt   = fltarr(sz, sz)
     minim= fltarr(sz, sz)
     maxim= fltarr(sz, sz)
     amask= intarr(sz, sz)
     omask= intarr(sz, sz)
  endelse

  flux[sind]  = fluxlist
  fraw[sind]  = rawlist
  wt[sind]    = wtlist
  minim[sind] = minlist
  maxim[sind] = maxlist
  amask[sind] = andlist
  omask[sind] = orlist

  delvarx, fluxlist, rawlist, wtlist, minlist, maxlist, andlist, orlist, $
    sind, uv

  print, 'Time ', systime(1)-t0, ' seconds'
  print, 'Memory usage to this point:'
  help, /mem

; -------- check for NaNs (there should never be any)
  if total(~finite(flux)) GT 0 then stop
  if total(~finite(wt)) GT 0 then stop

  im  = flux/(wt + (wt eq 0))
  raw = fraw/(wt + (wt eq 0))
  delvarx, alist, dlist, flux, fraw

  return
end


pro doit

  lb=[344.0,23.0]


  lgal = lb[0]
  bgal = lb[1]
  

  goodlist = pair.ndiff GT 0

  wise_mosaic1b, im, raw, wt,nx=1000/2, sz=8000, lgal=lgal, bgal=bgal, dust=dust, hdr=hdr, goodlist=goodlist

  return
end



function wise_mosaic_beta, hdr

; -------- generate (x,y)
  nx = sxpar(hdr, 'NAXIS1')
  ny = sxpar(hdr, 'NAXIS2')

  xbox = (lindgen(nx, ny) mod nx)+0.d ; these are pixel *centers*
  ybox = (lindgen(nx, ny)  /  nx)+0.d
  
; -------- transform to (lambda, beta)
  extast, hdr, astr
  xy2ad, xbox, ybox, astr, l, b
  euler, l, b, lam, bet, 6

  return, bet
end

pro betafoo

  coeff = poly_fit(beta, im, 1)
  mdl = poly(beta, coeff)




  return
end
