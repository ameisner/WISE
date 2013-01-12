;+
; NAME:
;   wisetile_mosaic
;
; PURPOSE:
;   Make WISE mosaic on ISSA tile
;
; CALLING SEQUENCE:
;   wisetile_mosaic, im, imrat, imoffs, wt, nside=nside, nx=nx, sz=sz,
;   lgal=lgal, bgal=bgal, dust=dust, hdr=hdr, astr=astr, imref=imref,
;   wtref=wtref, index=index
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
;   2012-Feb-05 - Written by Douglas Finkbeiner, CfA
;
;----------------------------------------------------------------------
function binary_search, list, val
  
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


pro wisetile_mosaic, im, imrat, imoffs, wt, nside=nside, nx=nx, sz=sz, lgal=lgal, bgal=bgal, dust=dust, hdr=hdr, astr=astr, imref=imref, wtref=wtref, index=index, tilepath=tilepath, zsub=zsub, ecl=ecl, tileoffs=tileoffs

  t0 = systime(1)

  if ~ keyword_set(tilepath) then tilepath = '/n/home08/dfink/runwise/tile'

  maxrad = 9.0  ; max radius to corner of tile image

; -------- initialize mosaic projection WCS header (if not healpix)
  if keyword_set(nside) then begin 
     healpix = 1
     npix = 12LL * nside * nside
     print, 'HEALPIX nside:', nside, '   npix:', npix
  endif else begin 
     healpix = 0
     if keyword_set(astr) then begin 
        print, 'using the astrometry structure you passed'
     endif else begin 
        cd = 2.75 * (1000/nx) / 3600.d
        dummy = fltarr(sz, sz, /nozero)
;        hdr = wisetile_mosaic_gal_header(dummy, cd, lgal=lgal, bgal=bgal)
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
        euler, lbox[k*nrow:(k+1)*nrow-1], bbox[k*nrow:(k+1)*nrow-1], atemp, dtemp, $ 
          select = keyword_set(ecl) ? 4 : 2
        abox[k*nrow:(k+1)*nrow-1] = temporary(atemp)
        dbox[k*nrow:(k+1)*nrow-1] = temporary(dtemp)
     endfor
     print
     delvarx, lbox, bbox
  endif else begin 
     nslice = 16
     xbox = (lindgen(sz, sz/nslice) mod sz)+0.d  ; these are pixel *centers*
     ybox = (lindgen(sz, sz)  /  sz)+0.d
     abox = dblarr(sz, sz)
     dbox = dblarr(sz, sz)
     for k=0L, nslice-1 do begin     ; loop for memory efficiency
        print, '.', format='($,A)'
        nrow = sz/nslice  ; rows per slice
        xy2ad, xbox, ybox[*, k*nrow:(k+1)*nrow-1], astr, atemp, dtemp
        if galcoords then euler, atemp, dtemp, select=2 ; atemp,dtemp were Galactic l,b
        if eclcoords then euler, atemp, dtemp, select=4 ; atemp,dtemp were ecliptic lam,beta
        abox[*, k*nrow:(k+1)*nrow-1] = atemp
        dbox[*, k*nrow:(k+1)*nrow-1] = dtemp
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

  alist = abox[sind]  ; may sort this later
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

  fluxlist = fltarr(n_elements(dlist))
;  rawlist  = fltarr(n_elements(dlist))
  offslist  = fltarr(n_elements(dlist))
  ratlist  = fltarr(n_elements(dlist))
  wtlist   = fltarr(n_elements(dlist))

; -------- get list of files to process
; -------- read (ra, dec, filename) structure
  if size(index, /TNAME) EQ 'STRUCT' then indstr = index $ 
    else indstr = mrdfits(index, 1)
  if ~ keyword_set(indstr) then message, 'bad index file name'
  euler, lgal, bgal, racen, deccen, 2
  dangle = djs_diff_angle(indstr.ra, indstr.dec, racen, deccen)
  angsind = sort(dangle)

  indstr = indstr[angsind]
  flist = indstr.fname

; -------- loop over files and add them to mosaic. 
  npad = 0  ; for the 512x512 clean files!
;  print, concat_dir(tilepath, flist[0])
  hdr = headfits(concat_dir(tilepath, flist[0]))
  nc = sxpar(hdr, 'NAXIS1')  ; assume image is square

; -------- make (x,y) coords for WISETILE image
  xbox = (lindgen(nc, nc) mod nc)+0.d ; these are pixel *centers*
  ybox = (lindgen(nc, nc)  /  nc)+0.d

  for i=0L, n_elements(flist) -1 do begin 
;  for i=0L, 50 do begin 

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
        print, i, concat_dir(tilepath, flist[i]), format='(I8,"  ",A)'
        ind0 = binary_search(dlist, crval[1]-maxrad)
        ind1 = binary_search(dlist, crval[1]+maxrad)
        
        uvcen = float(ll2uv(transpose(crval)))
        near = where(transpose(uvcen)##uv[ind0:ind1, *] GT cos(maxrad*!dpi/180.d), nnear)+ind0
     endif else nnear = 0
     
     if nnear GT 0 then begin 

        cleanname = concat_dir(tilepath, flist[i])
        if file_test(cleanname) eq 0 then message, 'cannot find file'
        clean = readfits(cleanname, hclean, /silent, exten=0)

        if keyword_set(zsub) then begin
           zody = zody_image(hclean, 12.0, binfac=30);*0.4 ; fudge factor !!!!!
           zody -= median(zody)
           clean -= zody
        endif
        offs = keyword_set(zsub) ? zody : fltarr(nc, nc)

;        offs  = readfits(cleanname, /silent, exten=2)
; --------  what was the original intended meaning/use of rat ???
        rat   = fltarr(nc, nc)
;        dirty = readfits(cleanname, hclean, /silent, exten=1)

        extast, hclean, aout
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
           intoffs = interpolate(offs, xx[w], yy[w])
           intrat  = interpolate(rat, xx[w], yy[w])
;           intdirt = interpolate(dirty, xx[w], yy[w])
           wwt = where(wtlist[nearw] GE 0.5, nwwt)
           if nwwt GT 100 then begin 
              dflux = median(intflux[wwt] - (fluxlist[nearw[wwt]]/wtlist[nearw[wwt]]), /even)
           endif else dflux = 0
           
 dflux = 0
           if keyword_set(tileoffs) then dflux = -1*sxpar(hclean, 'OFFSET')

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
           endif
           
           fluxlist[nearw] += (intflux-dflux)
           offslist[nearw] += intoffs
           ratlist[nearw] += intrat
;           rawlist[nearw]  += (intdirt-dflux)
           wtlist[nearw]++

        endif
     endif
  endfor

  print, 'mapping back to 2D'

  if healpix then begin
     rat  = fltarr(npix)
     offs = fltarr(npix)
     flux = fltarr(npix)
     fraw = fltarr(npix)
     wt   = fltarr(npix)
  endif else begin
     flux = fltarr(sz, sz)
     fraw = fltarr(sz, sz)
     wt   = fltarr(sz, sz)
  endelse

  flux[sind] = fluxlist
  offs[sind] = offslist
  rat[sind] = ratlist
;  fraw[sind] = rawlist
  wt[sind]   = wtlist


  delvarx, fluxlist, wtlist, sind, uv

  print, 'Time ', systime(1)-t0, ' seconds'
  print, 'Memory usage to this point:'
  help, /mem

; -------- check for NaNs (there should never be any)
  if total(~finite(flux)) GT 0 then stop
  if total(~finite(wt)) GT 0 then stop

  im  = flux/(wt + (wt eq 0))
  imrat   = rat/(wt + (wt eq 0))
  imoffs  = offs/(wt + (wt eq 0))
;  raw = fraw/(wt + (wt eq 0))
  delvarx, alist, dlist, flux

  return
end


pro doit, nside=nside

  band = 'wise12'
  if ~ keyword_set(nside) then nside = 4096
  index = 'wisetile-index.fits'
  wisetile_mosaic, im, imrat, imoffs, wt, nside=nside, dust=dust, lgal=0, bgal=0, index=index

  outname = string(band, '_', nside, '.fits', format='(A,A,I4.4,A)')
  mwrfits, im, outname, /create
  mwrfits, imoffs, outname
  mwrfits, wt, outname

  return
end

pro akari_mosaics, ecl=ecl

  basedir = concat_dir(getenv('AKARI_DATA'), 'FITS/Release1.0')
  bands  = ['WideS', 'WideL', 'N60', 'N160']
  dtype = ['', 'Nsamp', 'Nscan', 'sigma']
  extname = ['flux density', 'Nsamp', 'Nscan', 'sigma']

  index = mrdfits('$WISE_DATA/wisetile-index-allsky.fits', 1)
  tnum = string(1 + lindgen(430), format='(I03)')
  nside = 4096

  t0 = systime(1)
  for i=0, n_elements(bands)-1 do begin
; ----- start with a dummy nside=4096 header
    hdr = headfits('$AKARI_DATA/sfd_4096.fits')
    outname = $ 
      concat_dir(getenv('AKARI_DATA'), 'akari-' + bands[i]  + '-4096.fits')
    for j=0, n_elements(dtype)-1 do begin
      tilepath = concat_dir(basedir, bands[i])
      fname = tnum + '_' + bands[i]
      if ~(dtype[j] EQ '') then begin
        tilepath = concat_dir(tilepath, dtype[j])
        fname += '_' + strlowcase(dtype[j])
      endif
      fname += '.fits'
      index.fname = fname
      wisetile_mosaic, im, imrat, imoffs, wt, nside=nside, dust=dust, $ 
        lgal=0, bgal=0, index=index, tilepath=tilepath, ecl=ecl
      sxaddpar, hdr, 'EXTNAME', extname[j]
      if j EQ 0 then begin
        writefits, outname, im, hdr
; ----- modify header for extensions
        sxaddpar, hdr, 'XTENSION', 'Image', ' IMAGE extension',before='SIMPLE'
        sxdelpar, hdr, 'SIMPLE'
        sxdelpar, hdr, 'EXTEND'
      endif else $ 
        writefits, outname, im, hdr, /append
    endfor
    sxaddpar, hdr, 'EXTNAME', 'weight'
    writefits, outname, wt, hdr, /append
  endfor
  print, 'dt = ', systime(1) - t0, 'seconds'

end
