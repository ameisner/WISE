;+
; NAME:
;   wise_healpix_index
;
; PURPOSE:
;   Generate list of healpix indices and flux values for each image
;
; CALLING SEQUENCE:
;   wise_healpix_index, indlist, fluxlist, nside=, lb=, nimage=, cleanpath=
;
; INPUTS:
;   nside     - healpix nside
;   lb        - (l, b) of field center
;   nimage    - number of images to use
;   cleanpath - path to clean image files
;
; OUTPUTS:
;   indlist   - array (npix,nimage) of healpix indices
;   fluxlist  - flux values corresponding to indlist
;
; OPTIONAL OUTPUTS:
;   
; EXAMPLES:
;   
; COMMENTS:
;   healpix indices are ring-ordered. 
;
; REVISION HISTORY:
;   2012-Feb-25 - Written by Douglas Finkbeiner, CfA
;                  (some of it taken from wise_mosaic1b.pro)
;
;----------------------------------------------------------------------
pro wise_healpix_index, indlist, fluxlist, nside=nside, lb=lb, $
                        nimage=nimage, cleanpath=cleanpath, indstr=indstr, $ 
                        allsky=allsky, warp=warp, w4=w4

  t0 = systime(1)

  if ~ keyword_set(cleanpath) then cleanpath = '$WISE_DATA/clean'

  maxrad = 0.6  ; L1b

; -------- healpix 
  npix = 12LL * nside * nside
  print, 'HEALPIX nside:', nside, '   npix:', npix

; -------- compute (RA, dec) for each pixel in mosaic
  print, 'Preparing healpix coordinates...'

  healgen_lb, nside, lbox, bbox
  abox = fltarr(npix) ; try this with float -- should really use double!!
  dbox = fltarr(npix)
  nslice = 16
  for k=0L, nslice-1 do begin   ; loop for memory efficiency
     print, '.', format='($,A)'
     nrow = npix/nslice         ; rows per slice
     euler, lbox[k*nrow:(k+1)*nrow-1], bbox[k*nrow:(k+1)*nrow-1], atemp, dtemp, 2
     abox[k*nrow:(k+1)*nrow-1] = temporary(atemp)
     dbox[k*nrow:(k+1)*nrow-1] = temporary(dtemp)
  endfor
  print

  delvarx, lbox, bbox
  print, 'Sorting declinations...'
  sind = sort(dbox)

  alist = abox[sind]
  dlist = dbox[sind]
  delvarx, abox, dbox

; -------- compute unit vectors
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

  npsave = round(192*(nside/1024.)^2)

; -------- read metadata table, or use the one specified
  if ~keyword_set(indstr) then $ 
    indstr = wise_index_metadata(lb, nimage=nimage, allsky=allsky) $ 
    else nimage = keyword_set(nimage) ? nimage : n_elements(indstr)
  flist  = indstr.fname

; -------- allocate output arrays
  fluxlist = fltarr(npsave, nimage)
  indlist  = lonarr(npsave, nimage)-1

; -------- loop over files and add them to mosaic. 
  npad = 0  ; for the 512x512 clean files!

  for i=0L, nimage - 1 do begin 

     crval = [indstr[i].ra, indstr[i].dec]
     
; -------- Get indices for relevant dec range     
     print, i, fileandpath(flist[i]), format='(I8,"  ",A)'
;     ind0 = binary_search(dlist, crval[1]-0.90)
;     ind1 = binary_search(dlist, crval[1]+0.90)
     binary_search, dlist, crval[1]-0.90, ind0
     binary_search, dlist, crval[1]+0.90, ind1
     if (ind0 EQ -1) then ind0 = 0L
     if (ind1 EQ -1) then ind1 = n_elements(dlist)-1
  
     uvcen = float(ll2uv(transpose(crval)))
     near = where(transpose(uvcen)##uv[ind0:ind1, *] GT cos(maxrad*!dpi/180.d), nnear)+ind0
  
     if nnear GT 0 then begin 

        cleanname = wise_name('clean', flist[i], cleanpath=cleanpath)
        if file_test(cleanname) eq 0 then message, 'cannot find file'
        clean = lfs_fits_access(cleanname, hclean, /silent, exten=0)
        if n_elements(clean) LE 1 then continue
        if keyword_set(warp) then $ 
            wise_exposure_warp, clean, indstr[i].grad, w4=w4

        nc = (size(clean, /dimen))[0]

; we cannot trust the zero point, so subtract it for now. 
        extast, hclean, aout
        ;medclean = median(clean)
        ;clean -= medclean
        wnan = where(finite(clean) EQ 0, nwnan)
        if nwnan GT 0 then clean[wnan] = median(clean)   ;?????
     
        ad2xy, alist[near], dlist[near], aout, xx, yy  ; (xx, yy) are pixel centers

        w = where(xx gt npad and xx lt (nc-npad-1) and yy gt npad and yy lt (nc-npad-1), nw)
        if nwnan GT n_elements(clean)/2 then begin 
           message, /info, 'too many NaNs'
           nw = 0
        endif

        if nw gt 0 then begin 

           nearw = near[w]
           intflux = interpolate(clean, xx[w], yy[w])

; -------- store values 
           fluxlist[0:nw-1, i] = intflux
           indlist[0:nw-1, i] = sind[nearw]
           
        endif
     endif
  endfor

; -------- clean up
  delvarx, sind, uv

  print, 'Time ', systime(1)-t0, ' seconds'
  print, 'Memory usage to this point:'
  help, /mem

  return
end



pro quicklook, flux, ind, tot, wt, nside=nside

  if ~keyword_set(nside) then nside = 2048
  npix = 12LL * nside * nside
  sz = size(ind, /dimen)
  npsave = sz[0]
  nimage = sz[1]

  wt  = fltarr(npix)
  tot = fltarr(npix)
  for i=0L, nimage-1 do begin 
     w = where(ind[*, i] GE 0)
     wt[ind[w, i]]++
     tot[ind[w, i]] += flux[w, i]
  endfor

  
  return
end
