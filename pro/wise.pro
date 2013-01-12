; good test images
pro testimages

; -------- no cirrus
  atlas=readfits('./2571p772_aa11/2571p772_aa11-w3-int-3.fits')

; -------- some cirrus
  atlas=readfits('./1440p848_aa11/1440p848_aa11-w3-int-3.fits')

; -------- bright star
  atlas=readfits('./0194p726_aa11/0194p726_aa11-w3-int-3.fits')

; -------- bright cirrus
  atlas=readfits('./0389p726_aa11/0389p726_aa11-w3-int-3.fits')


  return
end



function wise_bitmask, im, nedge=nedge

  satlevel = 1500

; -------- input image has missing pixels set to NaN
  nodata = finite(im) EQ 0

; -------- mask nedge pixels around edge of image
  if keyword_set(nedge) then begin 
     sz = size(im, /dimen)
     edge = bytarr(sz[0], sz[1])
     edge[0:nedge-1, *] = 1B
     edge[*, 0:nedge-1] = 1B
     edge[sz[0]-nedge:sz[0]-1, *] = 1B
     edge[*, sz[1]-nedge:sz[1]-1] = 1B
  endif else edge = 0B
  
; -------- mask saturated pixels
  satur = im GE satlevel
  
  kpix = 45     ; make this an odd number
  kern = shift(dist(kpix), kpix/2, kpix/2) LT (kpix/2+0.5)

  halo = dilate(satur, kern)

; -------- define mask bits.  Leave 4 for faint sources. 
  bitmask = nodata + $
    edge*2B + $
    halo*8B + $
    satur*16B


  return, bitmask
end




pro wise_find_stars, im, nstar, xstar, ystar, flux, badpixels=badpixels

;  psfvals = [0.8, 0.7]
; we have loosened this to include double stars, etc.  For a clean list
; of stars, use [0.8, 0.7] and ">"

  psfvals = [0.88, 0.78]

  ivar = im*0+1.0  ; ?!???
  if ~ keyword_set(badpixels) then badpixels = finite(im) EQ 0

  if NOT keyword_set(im)        then message, 'must set image'
  if NOT keyword_set(ivar)      then message, 'must set ivar'
  if NOT keyword_set(npad)      then npad = 2
  if NOT keyword_set(nsigma)    then nsigma = 20
  if NOT keyword_set(badpixels) then message, 'works better if you set badpixels'


; HARDWIRE sigma!
  sigma = 0.2
  cut = 8*sigma

  medscale = 5
  sz = size(im, /dimen)

; -------- median filter image, look at pixels high relative to med
  med = median(im, medscale)  ; slow step

  high = (im-med) GT cut
  ind = where(high, nhigh)

  if nhigh eq 0 then begin
     splog, 'no high pixels!'
     return
  endif

; -------- (x,y) pixel indices for pixel and its neighbors
  ix = ind mod sz[0]
  iy = ind / sz[0]

  ixm = (ix-1) > 0
  ixp = (ix+1) < (sz[0]-1)

  iym = (iy-1) > 0
  iyp = (iy+1) < (sz[1]-1)


  maxnb = im[ixm, iy] > im[ixp, iy] > im[ix, iym] > im[ix, iyp] > $
    im[ixm, iyp] > im[ixp, iyp] > im[ixp, iym] > im[ixm, iym]

  bd = (ivar EQ 0) OR badpixels
  badnb = bd[ixm, iy] OR bd[ixp, iy] OR bd[ix, iym] OR bd[ix, iyp] OR $
    bd[ixm, iyp] OR bd[ixp, iyp] OR bd[ixp, iym] OR bd[ixm, iym]

; -------- demand that the peak is larger than the neighbor
  peak = (im[ind] GT (maxnb*1.00001)) AND (badnb EQ 0)
  peak = peak AND (badpixels[ind] EQ 0)

; -------- and not near the edge of the image...
  peak = peak AND ((iy GT npad) AND (iy LT (sz[1]-npad-1)) AND $
    (ix GT npad) AND (ix LT (sz[0]-npad-1)))


; -------- real peaks (mostly) - still need to worry about
;          diff. spikes, badcols
  w = where(peak, npeak)

; -------- we compute these ratios to see if candidate consistent with PSF
;
; 02 12 22
; 01 11 21
; 00 10 20 
;
; diff = im-median(im, medscale*3)

; -------- bin down and median, for speed
  xbox = lindgen(sz[0], sz[0]) mod sz[0]
  ybox = lindgen(sz[0], sz[0])  /  sz[0]

  bord = 4
  ns = 256
  na = ns+2*bord
  thumb = fltarr(na, na)

  thumb[bord:ns+bord-1, bord:ns+bord-1] = rebin(im, ns, ns)
  thumb[0:bord-1, *] = rebin(thumb[bord, *], bord, na)
  thumb[bord+ns:na-1, *] = rebin(thumb[ns+bord-1, *], bord, na)
  thumb[*, 0:bord-1] = rebin(thumb[*, bord], na, bord)
  thumb[*, bord+ns:na-1] = rebin(thumb[*, ns+bord-1], na, bord)

  thumbmed = (median(thumb, 3))
  bigmed = interpolate(thumbmed, xbox/4.+bord-0.375, ybox/4.+bord-0.375) ; not general
  diff = im-temporary(bigmed)

  im00 = diff[ixm[w], iym[w]]
  im10 = diff[ix[w],  iym[w]]
  im20 = diff[ixp[w], iym[w]]

  im01 = diff[ixm[w], iy[w]]
  im11 = diff[ix[w],  iy[w]]
  im21 = diff[ixp[w], iy[w]]

  im02 = diff[ixm[w], iyp[w]]
  im12 = diff[ix[w],  iyp[w]]
  im22 = diff[ixp[w], iyp[w]]

  back1 = (im01+im21)/2.
  back2 = (im10+im12)/2.
  diag1 = (im00+im22)/2.
  diag2 = (im02+im20)/2.

  xstar = float(ix[w])
  ystar = float(iy[w])

; -------- call djs_phot to get fluxes and recenter positions
  flux = djs_phot(xstar, ystar, 3.5, [4.5, 6.], diff)
  badstar = ((back1 < back2)/im11 gt psfvals[0]) OR ((diag1 < diag2)/im11 GT psfvals[1]) OR (flux LE 0)
  wgood = where(badstar EQ 0, nstar)
  wbad  = where(badstar, nbad)

;  ww = where(sqrt((xstar-96.*2)^2+(ystar-222*2)^2) lt 6)
;  print, back1[ww]/im11[ww]
;  print, back2[ww]/im11[ww]

;  print, diag1[ww]/im11[ww]
 ; print, diag2[ww]/im11[ww]


  if nstar GT 0 then begin 
     xstar = xstar[wgood]
     ystar = ystar[wgood]
     flux  = flux[wgood]
  endif else begin 
     xstar = -1
     ystar = -1
     flux  = -1
  endelse

;  atv,im,/al,/s
;  atvplot,ix[w],iy[w],ps=6,color='green',syms=2
;  atvplot,ix[w[wh]],iy[w[wh]],ps=6,color='magenta',syms=2.5

;  mask = dilate(high, bytarr(3, 3)+1B)
;  mask = dilate(mask, bytarr(3, 3)+1B)
;  psub = im*(1-mask)

  return
end


pro goo

  flist=file_search('clean2/*fits')

  for i=0, 200 do begin 
     print, i
     atv,readfits(flist[i]),/al,/s
     wait, 1
  endfor

  return
end



function wise_atlas_clean, atlas, h, hout=hout, bitmask=bitmask, dirty=dirty

; -------- determine image size, subtract median
  sz = size(atlas, /dimen)
  npix = 4096
  if sz[0] NE (npix-1) or sz[1] NE (npix-1) then stop
  full = fltarr(npix, npix)
  imagemed = median(atlas)
  full[0:npix-2, 0:npix-2] = atlas   

; -------- pad out to 4096 with linear interpolation
  full[npix-1, *] = full[npix-2, *]*2 - full[npix-3, *]
  full[*, npix-1] = full[*, npix-2]*2 - full[*, npix-3]

; -------- bin down to ~ 2.4 pix per FWHM (barely well sampled at 12 mu)
  binned = rebin(temporary(full), 1024, 1024)
  binned -= imagemed

  bitmask  = wise_bitmask(binned, nedge=8)
  wmissing = where(bitmask AND 1B, nmissing)
  if nmissing GT 0 then binned[wmissing] = 0  ; set missing data to zero. 
  nx = 512
  dirty = rebin(fastconv(binned, 3, 1, /nodisp, /edge, /sil), nx, nx)

  wbright = where((bitmask AND 8B) NE 0, nbright)
  if nbright GT 0 then binned[wbright] = 0  ; set missing data to zero. 

; -------- define bad pixels for wise_find_stars
  badpixels = (bitmask AND 25B) NE 0

; -------- find stars
  wise_find_stars, binned, nstar, xstar, ystar, flux, badpixels=badpixels

; -------- remove stars
  psub = binned
  mask = badpixels
  wise_remove_star, psub, xstar, ystar, flux, mask=mask
  wise_remove_star, psub, xstar, ystar, flux, mask=mask ; call twice

  bitmask = bitmask OR (mask AND (~badpixels))*4B

; -------- filter a little more, smooth with a Guassian, rebin
  clean = rebin(fastconv(median(psub, 3), 3, 1, /nodisp, /edge, /sil), nx, nx)

; -------- rebin bitmask       THIS IS NOT RIGHT YET, should really dilate it!
  cmask = uintarr(nx, nx)  ; allow 16 bits
  for i=0, 7 do cmask = cmask OR (~rebin((bitmask AND (2B^i)) EQ 0B, nx, nx))*2U^i
  bitmask = cmask

; -------- update astrometry after rebinning
  extast, h, astr
  astr.naxis = [1, 1]*nx
  astr.cdelt *= (4096/nx)
  astr.crpix /= (4096/nx)

; -------- overwrite WCS header and return
  hout = h
  putast, hout, astr
  sxaddpar, hout, 'NAXIS1', nx
  sxaddpar, hout, 'NAXIS2', nx
  sxaddpar, hout, 'PXSCAL1', 1.375d * (4096/nx)
  sxaddpar, hout, 'PXSCAL2', 1.375d * (4096/nx)

  return, clean
end




pro playaround

  xbox = (lindgen(nx, nx) mod nx)+0.5
  ybox = (lindgen(nx, nx)  /  nx)+0.5
  xy2ad, xbox, ybox, astr, abox, dbox

  euler, abox, dbox, lbox, bbox, 1
  dust = dust_getval(lbox, bbox, /noloop, /interp)

; -------- postscript file
  loadct, 0
  dfpsplot, 'wise.ps', bits=8, xsize=6
  sys = sysvars(/print)
  !p.multi = [0, 1, 2]

  display, clean, abox[*, 0], dbox[0, *], xtit='RA [deg]', ytit='dec [deg]', $
    title='WISE 12 micron, cleaned', min=-1., max=3., top=255

  display, dust, abox[*, 0], dbox[0, *], xtit='RA [deg]', ytit='dec [deg]', $
    /noerase, title='SFD dust'

  restore_sysvars, sys
  dfpsclose

  return
end



pro wise_psub, im


  cut = 1.
  medscale = 5

  med = median(im, medscale)
  high = (im-med) GT cut
  mask = dilate(high, bytarr(3, 3)+1B)
  mask = dilate(mask, bytarr(3, 3)+1B)
  
  psub = im*(1-mask)




  return
end




pro makeclean

  flist = file_search('*/*w3-int-3.fits', count=nfile)
  base = strmid(flist, 0, 8)  ; dangerous

  splog, nfile, ' files found.'

  clean_dir = '$WISE_DATA/clean3'
  file_mkdir, clean_dir
  splog, 'Working in ', clean_dir
  
  for i=0, n_elements(flist)-1 do begin 

     outname = concat_dir(clean_dir, base[i]+'.fits')
     print, strcompress(string(i+1, ' of ', nfile, '  ', outname))

     if file_test(outname) eq 0 then begin 
        raw   = readfits(flist[i], hdr, /silent)
        clean = wise_atlas_clean(raw, hdr, hout=hout, bitmask=bitmask, dirty=dirty)
        if total(~finite(clean)) GT 0 then message, 'NaN found in clean!'
        writefits, outname, clean, hout
        mwrfits, dirty, outname, hout, /silent   ; append uncleaned image
        mwrfits, bitmask, outname, hout, /silent   ; append bitmask
     endif

  endfor

  return
end


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



function wise_mosaic_header, flux, cd

  sz = (size(flux, /dimen))[0]
; -------- define header for mosaic
  mkhdr, hdr, flux

  sxaddpar, hdr, 'EQUINOX', 2000.0
  sxaddpar, hdr, 'CTYPE1', 'RA---ZEA'
  sxaddpar, hdr, 'CTYPE2', 'DEC--ZEA'
  sxaddpar, hdr, 'CD1_1', -cd ; / Degrees / Pixel
  sxaddpar, hdr, 'CD2_1', 0.0
  sxaddpar, hdr, 'CD1_2', 0.0
  sxaddpar, hdr, 'CD2_2', cd  ;  / Degrees / Pixel
  sxaddpar, hdr, 'CRPIX1', sz/2 ; Reference Pixel in X
  sxaddpar, hdr, 'CRPIX2', sz/2 ; Reference Pixel in Y
  sxaddpar, hdr, 'CRVAL1', 0.00000000000 ;  / R.A. (degrees) of reference pixel
  sxaddpar, hdr, 'CRVAL2', 90.000; / Declination of reference pixel
  sxaddpar, hdr, 'LONPOLE', 180.000000000 ; / Native longitude of Celestial pole
  sxaddpar, hdr, 'LATPOLE', 90.0000000000 ; / Celestial latitude of native pole

  return, hdr
end



pro doakari
  
  fname = '~/AKARI/akaricorr2-090.fits'
  junk = readfits(fname, h)
  extast, h, astr

  sz = 1600
  cd = abs(astr.cd[0, 0])
  nx = 4096/(cd*3600.d/1.375)

  mosaic, im, raw, nx=nx, sz=sz, dust=dust, hdr=hdr, astr=astr


  return
end




pro domosaic

;    lgal    bgal
; p1  206    -26
; p2  231,   -30
; p3  344,    23
; p4  196    -10

  nx = 512
  grow = 10
  sz = nx*grow

  mosaic, im, raw, nx=nx, sz=sz, dust=dust, hdr=hdr, lgal=344., bgal=23.

  num = '3'
  writefits, 'imp'+num+'.fits', im, hdr
  writefits, 'rawp'+num+'.fits', raw, hdr
  writefits, 'dustp'+num+'.fits', dust, hdr


  return
end



pro doheal

  nside = 8192L

  mosaic, im, raw, nside=nside

  writefits, 'im8h.fits', im
  writefits, 'raw8h.fits', raw

  delvarx, im, raw
  healgen_lb, nside, lbox, bbox
  dust = dust_getval(lbox, bbox, /noloop, /interp)
  delvarx, lbox, bbox

  writefits, 'dust8h.fits', dust


  return
end




pro mosaic, im, raw, nside=nside, nx=nx, sz=sz, lgal=lgal, bgal=bgal, dust=dust, hdr=hdr, astr=astr

  t0 = systime(1)

  if keyword_set(nside) then begin 
     healpix = 1
     npix = 12LL * nside * nside
     print, 'HEALPIX nside:', nside, '   npix:', npix
  endif else begin 
     healpix = 0
     if keyword_set(astr) then begin 
        print, 'using the astrometry structure you passed'
     endif else begin 
        cd = 1.375 * (4096/nx) / 3600.d
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
     asmall = rebin(abox, sz/2, sz/2, /sample)
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
  rawlist  = fltarr(n_elements(dlist))
  wtlist   = fltarr(n_elements(dlist))

; -------- get list of files to process
  clean_dir = '$WISE_DATA/clean3'
  flist = file_search(clean_dir+'/*.fits')

; -------- loop over files and add them to mosaic. 
  npad = 4  ; for the 512x512 clean files!

  for i=0, n_elements(flist) -1 do begin 
;  for i=0, 240-1 do begin 

     print, i, flist[i], format='(I6,"  ",A)'
     hclean = headfits(flist[i])
     if n_elements(clean) EQ 1 then begin
        print, 'missing extension, skipping...'
        nnear = 0  ; force skip
     endif else begin 
        nc = (size(clean, /dimen))[0]
        extast, hclean, aout

; -------- check if crval is even in the target projection
        if healpix then goodcen = 1 else $ 
        begin 
           if galcoords then begin ; if GAL projection, convert to (l,b) first
              euler, aout.crval[0], aout.crval[1], thisl, thisb, 1
              ad2xy, thisl, thisb, astr, xcen, ycen
           endif else begin 
              if eclcoords then begin 
                 euler, aout.crval[0], aout.crval[1], thisl, thisb, 3 ; ecliptic
                 ad2xy, thisl, thisb, astr, xcen, ycen
              endif else begin 
                 ad2xy, aout.crval[0], aout.crval[1], astr, xcen, ycen
              endelse
           endelse
           goodcen = (xcen GE 0) AND (xcen LT sz) AND (ycen GE 0) AND (ycen LT sz)
        endelse
        
; -------- Get indices for relevant dec range     
        if goodcen then begin 
           ind0 = binary_search(dlist, aout.crval[1]-0.90)
           ind1 = binary_search(dlist, aout.crval[1]+0.90)
           
           uvcen = float(ll2uv(transpose(aout.crval)))
           near = where(transpose(uvcen)##uv[ind0:ind1, *] GT cos(1.15d/!radeg), nnear)+ind0
        endif else nnear = 0
     endelse

     if nnear GT 0 then begin 
        clean = readfits(flist[i], hclean, /silent, exten=0)
     
        ad2xy, alist[near], dlist[near], aout, xx, yy  ; (xx, yy) are pixel centers

        w = where(xx gt npad and xx lt (nc-npad-1) and yy gt npad and yy lt (nc-npad-1), nw)
        if nw gt 0 then begin 
           if total(~finite(clean)) GT 0 then stop

; -------- if we have made it this far, read the dirty image also
           dirty = readfits(flist[i], hdr2, /silent, exten=1)

           nearw = near[w]
           fluxlist[nearw] += interpolate(clean, xx[w], yy[w])
           rawlist[nearw]  += interpolate(dirty, xx[w], yy[w])
           wtlist[nearw]++
        endif
     endif
  endfor

  print, 'mapping back to 2D'

  if healpix then begin
     flux = fltarr(npix)
     fraw = fltarr(npix)
     wt   = fltarr(npix)
  endif else begin
     flux = fltarr(sz, sz)
     fraw = fltarr(sz, sz)
     wt   = fltarr(sz, sz)
  endelse

  flux[sind] = fluxlist
  fraw[sind] = rawlist
  wt[sind]   = wtlist


  delvarx, fluxlist, rawlist, wtlist, sind, uv

  print, 'Time ', systime(1)-t0, ' seconds'
  print, 'Memory usage to this point:'
  help, /mem

; -------- check for NaNs (there should never be any)
  if total(~finite(flux)) GT 0 then stop
  if total(~finite(wt)) GT 0 then stop

  im  = flux/(wt + (wt eq 0))
  raw = fraw/(wt + (wt eq 0))
  delvarx, alist, dlist, flux, fraw, wt

  return
end


pro pixmask, flist

  im = readfits(flist[5], h)
  extast, h, astr
  colors = ['red', 'green', 'blue', 'magenta']
  atv, im, /al, /s

; order D, H, O, P
  for i=0, 3 do begin 
     readcol, flist[i], ra, dec
     ad2xy, ra, dec, astr, x, y
     atvplot, x, y, color=colors[i], ps=3
  endfor

  


  return
end


function fake, im, dust
  
  sz = size(im, /dimen)

  sm1 = fastconv(im, 16.6, 2, /nodisp, /edge)
  sm = rebin(sm1, sz[0], sz[1]) ; half pixel shift here
  out = 20*rebin(dust, sz[0], sz[0])+(im-sm)

  return, out
end


function selectk, im, k0, k1

  sz = size(im, /dimen)
  if sz[0] NE sz[1] then stop
  if sz[0] GT 2050 then message, 'you should run this on a subimage'

  f = fft(im)  ; should really pad this

  d = dist(sz[0])
  goodk = (d GE k0 and d LE k1)
  
  out = float(fft(f*goodk, /inv))

  return, out
end

pro testtt


  a=im[4942:4942+2047,9385:9385+2047]
  f=fft(a)
  help,f

  
  atv,shift(abs(f),1024,1024)

  d=dist(2048)
  pk=fltarr(1449)
  for i=0L,n_elements(d)-1 do pk[round(d[i])]+=abs(f[i])^2
  splot,pk

  sz = size(im, /dimen)
  sm1 = fastconv(im, 16.6, 2, /nodisp, /edge)
  sm = rebin(sm1, sz[0], sz[1]) ; half pixel shift here
  ssm = sm[4942:4942+2047,9385:9385+2047]
  ff=fft(ssm)
  pks=fltarr(1449)
  for i=0L,n_elements(d)-1 do pks[round(d[i])]+=abs(ff[i])^2

  dd = dust[4942:4942+2047,9385:9385+2047]

  goodk = (d GT 20) and (d LT 30)


  dust1 = selectk(dust, 20, 30)



  return
end



pro fig4, psname, im, raw, fake, dust, xr, yr, min=min, max=max, dif=dif

;  min = -0.2
;  max = 7.8
  dfpsplot, psname, bits=8, ysize=8.
  sys = sysvars(/print)
  loadct, 0

  !p.multi = [0, 2, 2]
  !p.charsize = .8
  !x.margin = [5, 1.5]
;  dif = 2.0
  display, bytscl(-raw[xr[0]:xr[1], yr[0]: yr[1]], min=-max+dif, max=-min+dif), $
    title=textoidl('WISE 12\mum')
  display, bytscl(-im[xr[0]:xr[1], yr[0]: yr[1]], min=-max+dif, max=-min+dif), /noerase, $
    title='point sources removed'
  display, bytscl(-20*dust[xr[0]/2:xr[1]/2, yr[0]/2: yr[1]/2], min=-max, max=-min), /noerase, $
    title='SFD dust'
  display, bytscl(-fake[xr[0]:xr[1], yr[0]: yr[1]], min=-max, max=-min), /noerase, $
    title='SFD / WISE hybrid'
  restore_sysvars, sys
  dfpsclose

  return
end


pro zoom4, psname, fake, dust, xcen, ycen, min=min, max=max

;  min = -0.2
;  max = 7.8
  dfpsplot, psname, bits=8, ysize=8.
  sys = sysvars(/print)
  loadct, 0

  d = [11, 1024, 512, 256]/2
  d0 = 256
  dc0 = xcen/2-d0
  dc1 = xcen/2+d0-1


  !p.multi = [0, 2, 2]
  !p.charsize = .8
  !x.margin = [7, 1.5]

  xx = findgen(512)*11/60.  ; arcminutes
  xx -= mean(xx)
  xtit = 'dRA [arcmin]'
  ytit = 'ddec [arcmin]'
;  dif = 2.0
  display, bytscl(-40*rebin(dust[dc0:dc1, dc0:dc1], 512, 512), min=-max, max=-min), xx*4, xx*4, title='SFD', xtit=xtit, ytit=ytit
  display, bytscl(-rebin(fake[xcen-d[1]:xcen+d[1]-1, ycen-d[1]:ycen+d[1]-1], 512, 512), min=-max, max=-min), $
    xx*4, xx*4, title=textoidl('WISE 12\mum'), /noe, xtit=xtit, ytit=ytit
  display, bytscl(-rebin(fake[xcen-d[2]:xcen+d[2]-1, ycen-d[2]:ycen+d[2]-1], 512, 512), min=-max, max=-min), xx*2, xx*2, /noe, xtit=xtit, ytit=ytit
  display, bytscl(-rebin(fake[xcen-d[3]:xcen+d[3]-1, ycen-d[3]:ycen+d[3]-1], 512, 512), min=-max, max=-min), xx, xx, /noe, xtit=xtit, ytit=ytit

  restore_sysvars, sys
  dfpsclose

  return
end


pro dofigs

; [-0.2,7.8]
  fig4, 'fig1.ps', im, raw, f, dust, [0,511]+9810, [0,511]+11110, dif=2
  fig4, 'fig2.ps', im, raw, f, dust, [0,511]+9738, [0,511]+12549, min=-0.836, max=4.71, dif=2
  fig4, 'fig3.ps', im, raw, f, dust, [0,1023]+7762, [0,1023]+6600, min=12.8, max=52.0, dif=15

;  fig4, 'fig4.ps', im, raw, f, dust, [0,1023]+7762, [0,1023]+6600, min=12.8, max=52.0, dif=15
  zoom4, 'figzoom.ps', f, dust, 1118, 1087, min=-15., max=50.

  return
end


pro cleancenters, ra, dec

  clean_dir = '$WISE_DATA/L3a/*'
  flist = file_search(clean_dir+'/*w3-int-3.fits', count=nfile)
  
  ra  = fltarr(nfile)
  dec = fltarr(nfile)
  for i=0L, nfile-1 do begin 
     if (i mod 100) eq 0 then print, i
     h = headfits(flist[i])
     ra[i]  = sxpar(h, 'CRVAL1')
     dec[i] = sxpar(h, 'CRVAL2')
  endfor

  return
end


pro resfig

  im = readfits('im3.fits', hdr)
  
  dust=readfits('dust3.fits')
  d=rebin(dust,15360/2,15360/2)
  f=fake(im,d)

  row = 10154
  x = [7300, 7799]
  nx = x[1]-x[0]+1

  y1 = rebin(f[x[0]:x[1], row:row+3], nx, 1)
  y0 = dust[x[0]:x[1], row]

  xx = findgen(nx)*22./60.

  dfpsplot, 'resfig.ps', bits=8, xsize=6, ysize=6, /color
  sys = sysvars(/print)
  !p.charsize = 1.2

  plot, xx, y1/20, xtitle='arcmin', ytit='E(B-V) [mag]', /xst, thick=3
  djs_oplot, xx, y0, color='blue', thick=6

  restore_sysvars, sys
  dfpsclose

  return
end



pro makepol

; get hdr


  nx = 256
  sz = 20480
  cd = 1.375 * (4096/nx) / 3600.d
   
  dummy = fltarr(sz, sz, /nozero)
  hdr = wise_mosaic_header(dummy, cd)
  dummy = 0


  im = readfits('im4.fits', hdr_foo)
  dust = readfits('dust4.fits')
  f = fake(im, dust)

  extast, hdr, astr

  atv, f, he=hdr

  heilespol, a
  ad2xy, a.ra, a.dec, astr, x, y
  eps = a.pol/10.  ; deg
;  dra  = cos(a.PA/!radeg)/cos(a.dec/!radeg)*eps
;  ddec = sin(a.PA/!radeg)*eps
  dra  = sin(a.PA/!radeg)/cos(a.dec/!radeg)*eps
  ddec = cos(a.PA/!radeg)*eps

  ad2xy, a.ra+dra, a.dec+ddec, astr, x1, y1
  ad2xy, a.ra-dra, a.dec-ddec, astr, x2, y2




  atvplot, x, y, ps=1, syms=1

  for i=0L, n_elements(x1)-1 do if (x[i] gt 0) AND (y[i] GT 0) AND (y[i] lt 20000) and (x[i] lt 12000) then atvplot, [x1[i], x2[i]], [y1[i], y2[i]], color='yellow'



  return
end


pro ghost

  npix = 128
  x0 = 4086
  y0 = 4818
  g1 = im[x0:x0+npix-1, y0:y0+npix-1]
  x0 = 4309
  y0 = 50
  g2 = im[x0:x0+npix-1, y0:y0+npix-1]
  atv,sshift2d(g2,[0.5,-1.8])-g1*1.2,/al,/s


  return
end


pro fooo

  num = '1'
  im   = readfits('imp'+num+'.fits', hdr)
  dust = readfits('dustp'+num+'.fits', hdr)
  f = fake(im, dust)




  return
end


pro akari_mosaic, flist, im, wt, dust, nu=nu, astr=astr0

; -------- read the header for the first file
  if 0 then begin 
     h = headfits(flist[9])  
     nx = 1600                  ; size of big projection
  endif else begin 
     h = headfits(flist[0])  
     nx = 1024
  endelse

; -------- extract WCS astrometry and modify to make new large projection
  extast, h, astr0
  astr0.crpix = [nx/2.+0.5, nx/2.+0.5]

; -------- get (lam,bet) within this projection
  xbox = lindgen(nx, nx) mod nx
  ybox = lindgen(nx, nx)  /  nx
  xy2ad, xbox, ybox, astr0, lam, bet

  sz = [0, 0]
  im = fltarr(nx, nx)
  wt = fltarr(nx, nx)

; -------- loop over image tiles
  for i=0, n_elements(flist)-1  do begin 
     tile  = readfits(flist[i], h)
     sz[0] = sxpar(h, 'NAXIS1')
     sz[1] = sxpar(h, 'NAXIS2')

     xbox = lindgen(sz[0], sz[1]) mod sz[0]
     ybox = lindgen(sz[0], sz[1])  /  sz[1]
     pad = 10
     mask = (xbox gt pad and xbox lt sz[0]-pad-1 and $
             ybox gt pad and ybox lt sz[1]-pad-1)
     smask = fastconv(float(mask), pad, 1, /nodisp, /edge, /silent)
     
     extast, h, astr
     ad2xy, lam, bet, astr, x, y
     ind = where(x ge 0 and x lt sz[0]-1 and y ge 0 and y lt sz[1]-1)
     wt[ind] += smask[x[ind], y[ind]]
     im[ind] += (smask*tile)[x[ind], y[ind]] ; should interpolate!!!
  endfor

  im = im/(wt + (wt eq 0))

  euler, lam, bet, l, b, 5
  dust = predict_thermal(l, b, nu=nu, /noloop, /interp) ; nu in GHz

  return
end



pro akariwrite
  fac = 1 ; factor to convert from SFD E(B-V) to Akari

; -------- get the images
  if 0 then begin 
;     flist = file_search('~/AKARI_FITS2/*WideS.fits')
     flist = file_search('~/AKARI/*WideS.fits')
     lam = 90.
  endif else begin 
     flist = file_search('~/AKARI/*WideL.fits')
     lam = 140.
  endelse 
  akari_mosaic, flist, im, wt, dust, nu=3d5/lam, astr=astr
  mkhdr, h, im
  putast, h, astr

; -------- smooth akari image to SFD resolution
  sm = fastconv(im,6.1*4,1,/nodisp,/edge)

; -------- correction image  
  corr = sm-dust*fac
  w = where(wt ne 0)
  med = median(corr[w])
  corr -= med
  
  akaricorr = (im-corr)*(wt ne 0)

  lamstr = string(lam, format='(I3.3)')
  writefits, 'akaricorr-'+lamstr+'.fits', akaricorr, h
  writefits, 'dust-'+lamstr+'.fits', dust*(wt ne 0), h
  writefits, 'akari-'+lamstr+'.fits', im, h

  return
end


pro akari_filter

  fname = '~/wise/FIS_RSRF_070122.txt'
  readcol, fname, lam1, freq1, f60, f90, lam2, freq2, f140, f160


  djs_plot, lam1, f60, xr=[40, 220], /xst, color='light blue', $
    xtickinterval=40, yr=[0, 1.4]
  djs_oplot, lam1, f90, color='blue'
  djs_oplot, lam2, f140, color='red'
  djs_oplot, lam2, f160, color='dark red'



  return
end



pro fourier

;  sub = akaricorr[1124-255:1124,1163-255:1163]
  sub = akaricorr[180:180+1023, 200:200+1023]
  dustsub = dust2[180:180+1023, 200:200+1023]
  f = fft(sub)

  nx = 1024
  mask = fltarr(nx, nx)+1
  y = (1+findgen(3))*37
  mask[*, 512+y] = 0
  mask[*, 512-y] = 0

  mask = shift(mask, nx/2, nx/2)
  clean = float(fft(mask*f, /inv))


  return
end
