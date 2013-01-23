;+
; NAME:
;   wise_make_mosaic
;
; PURPOSE:
;   Make a mosaic centered on lb, with catalog-based PSF subtraction
;
; CALLING SEQUENCE:
;   wise_make_mosaic, lb, nimage, outname, cleanpath=cleanpath
;
; INPUTS:
;   lb        - Galactic (l,b) [deg] of mosaic center
;   nimage    - use nearest nimage exposures to (l,b)
;   outname   - name for output FITS file
;
; KEYWORDS:
;   cleanpath - path for clean file tree
;
; OUTPUTS:
;   <files>
;
; EXAMPLES:
;   wise_make_mosaic,[220,0.],200,'foo.fits',cleanpath='/n/panlfs/ameisner/clean.bright.v8'
;
; COMMENTS:
;   Identify point sources in each image based on the WISE W3 catalog
;   and then subtracts them
;   
; REVISION HISTORY:
;   2011-Dec-09 - Written by Douglas Finkbeiner, CfA
;   2012-Feb-13 - Generalized for arbitrary location of cleaned files - AM
;   2012-Feb-21 - Include bitmask in mosaic - DPF
;   2012-Feb-25 - Use healpix-based pairwise_compare
;
;----------------------------------------------------------------------
pro wise_make_mosaic, lb, nimage, outname, cleanpath=cleanpath, skip=skip, $ 
                      allsky=allsky, astr=astr, warp=warp

  t0 = systime(1)
  if ~ keyword_set(cleanpath) then $ 
      cleanpath = keyword_set(allsky) ? '/n/wise/ameisner/clean.allsky.v0' : $ 
                                        '/n/panlfs/ameisner/clean.prelim.v0'
  if keyword_set(astr) then begin
      racen = astr.crval[0]
      deccen = astr.crval[1]
      euler, racen, deccen, lgal, bgal, 1
      lb = [lgal, bgal]
  endif

  if ~ keyword_set(skip) then begin 

; -------- write the cleaned files
     wise_write_clean1b, lb, nimage=nimage, outpath=cleanpath, allsky=allsky

; -------- read the files
;     imcube = wise_read_imcube(lb, nimage=nimage, cleanpath=cleanpath)
     wise_healpix_index, indh, flux, nside=2048, lb=lb, nimage=nimage, cleanpath=cleanpath, allsky=allsky, warp=warp

; -------- compute variance of pairwise image differences
;     pair = wise_pairwise_compare(lb, imcube=imcube, cleanpath=cleanpath)
     pair = wise_pairwise_compare(lb, indh, flux, allsky=allsky)
     
; -------- remove the worst outliers
     wise_pairwise_reject, pair
     
; -------- generate a mosaic
     goodlist = pair.ndiff GT 0
     
  endif else begin 
     goodlist = bytarr(nimage)+1B
  endelse

  splog, 'Time so far:', systime(1)-t0, ' seconds'

  wise_mosaic1b, im, raw, wt, minim, maxim, amask, omask, $
    nx=1000/2, sz=8000, lgal=lb[0], bgal=lb[1], dust=dust, $
    hdr=hdr, goodlist=goodlist, cleanpath=cleanpath, allsky=allsky, $ 
    astr=astr, warp=warp


; -------- subtract min, max images and reweight
  imclean = ((wt GT 2)*(wt*im-(minim+maxim)))/((wt-2) > 1)
  art = (wt GT 2)*(im-imclean)

  writefits, outname, imclean, hdr
  mwrfits, raw, outname, hdr, /silent
  mwrfits, wt, outname, hdr, /silent
  mwrfits, dust, outname, hdr, /silent
  mwrfits, minim, outname, hdr, /silent
  mwrfits, maxim, outname, hdr, /silent
  mwrfits, amask, outname, hdr, /silent
  mwrfits, omask, outname, hdr, /silent
  mwrfits, art,   outname, hdr, /silent

  splog, 'Finished after:', systime(1)-t0, ' seconds'

  return
end
