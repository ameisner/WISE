;+
; NAME:
;   akari_starlist
; PURPOSE:
;   identify bright point or slightly resolved sources in Akari WideS, 
;   WideL FIS data
;
; CALLING SEQUENCE:
;   akari_starlist, h, xlist, ylist, long=, cat=
;
; INPUTS:
;   h      - header of fits file for which a source list is desired
;   
; KEYWORDS:
;   long   - use Akari WideL catalog (if not set default to using
;            WideS catalog)
;   cat    - variable into which to store a structure containing all
;            fields present in Akari catalog, but only for the relevant
;            source subset
;
; OUTPUTS:
;   xlist  - list of source x pixel values
;   ylist  - corresponding list of source y pixel values
;   
; EXAMPLES:
;   see starlist_single_tile below
;
; COMMENTS:
;   catalog downloaded from:
;
;   http://www.ir.isas.ac.jp/AKARI/Observation/PSC/Public/
;
;   AKARI-FIS_BSC_V1.fits.gz
;
;   note that i generated the last four flags {fwhm_wides, fwhm_widel,
;   centroid_wides, centroid_widel} with 1 = good, 0 = bad because
;   these cuts remove a lot of bogus sources
;
; REVISION HISTORY:
;   2012-Jun-26 - Written by Aaron Meisner
;
;----------------------------------------------------------------------
pro akari_starlist, h, xlist, ylist, long=long, cat=cat

  MSKPIX = 21

  COMMON AKARICATALOG, cat_short, cat_long
  if n_elements(cat_short) EQ 0 then begin
    cat = mrdfits('$AKARI_DATA/MP_1001_catp.fits', 1)
    wlong = where(finite(cat.flux140) AND $ 
                  (cat.fqual140 EQ 3) AND $
                  (cat.flux140 GT 5*cat.ferr140) AND (cat.flux140 GT 3.) AND $ 
                  (cat.fwhm_widel) AND (cat.centroid_widel)) 
    wshort = where(finite(cat.flux90) AND $ 
                  (cat.fqual90 EQ 3)  AND $
                  (cat.flux90 GT 5*cat.ferr90) AND (cat.flux90 GT 1.) AND $ 
                  (cat.fwhm_wides) AND (cat.centroid_wides))
    cat_short = cat[wshort]
    cat_long = cat[wlong]
  endif

  cat = keyword_set(long) ? cat_long : cat_short

  racen = sxpar(h, 'CRVAL1')
  deccen = sxpar(h, 'CRVAL2')
  XPIX = sxpar(h, 'NAXIS1')
  YPIX = sxpar(h, 'NAXIS2')

  maxsep = 12.5/sqrt(2) ; assume tiles 12.5 deg x 12.5 deg
  
  catra = cat.ra
  catdec = cat.dec

  dangle = djs_diff_angle(catra, catdec, racen, deccen)
  wclose = where(dangle LE maxsep, nclose)
  if (nclose EQ 0) then return

  extast, h, astr
  ad2xy, catra[wclose], catdec[wclose], astr, x, y
  ix = round(x)
  iy = round(y)
  wtile = where((ix GE -MSKPIX/2) AND (ix LT XPIX+MSKPIX/2) AND $ 
                (iy GE -MSKPIX/2) AND (iy LT YPIX+MSKPIX/2), nwtile)
  if nwtile EQ 0 then return
  xlist = x[wtile]
  ylist = y[wtile]

  if arg_present(cat) then cat = cat[wclose[wtile]]

end

; test star list code by overplotting stars on a single tile
pro starlist_single_tile, tnum, im, long=long, cat=cat

  band = keyword_set(long) ? 'WideL' : 'WideS'
  if ((tnum LT 1) OR (tnum GT 430)) then return
  akaridata = getenv('AKARI_DATA')
  tilebase = concat_dir(akaridata, 'FITS/Release1.0')
  tilepath = concat_dir(tilebase, band)
  fname = string(tnum, format='(I03)') + '_' + band + '.fits'
  tname = concat_dir(tilepath, fname)
  
  print, tname
  im = readfits(tname, h, /silent)
  euler, sxpar(h, 'CRVAL1'), sxpar(h, 'CRVAL2'), lgal, bgal, 1
  print, '(lgal, bgal) = (',lgal,',',bgal,')'
  akari_starlist, h, xlist, ylist, long=long, cat=cat
  if n_elements(xlist) EQ 0 then begin
    print, 'no stars on this tile'
    return
  endif
  atv, im, header=h ; want wcs
  print, n_elements(xlist) 
  atvplot, xlist, ylist, ps=6, color='red', syms=2
end

pro test_akari_flags, long=long

  cat = mrdfits('$AKARI_DATA/MP_1001_catp.fits', 1, /silent)
  fqual = keyword_set(long) ? cat.fqual140 : cat.fqual90
  flags = keyword_set(long) ? cat.flags140 : cat. flags90

  good = long(total(fqual EQ 3))
  nondet = long(total((fqual EQ 3) AND (flags EQ -1)))
  sidelobe = long(total((fqual EQ 3) AND ((flags AND 8) NE 0)))
  toofaint = long(total((fqual EQ 3) AND ((flags AND 2) NE 0)))
  cdsmode = long(total((fqual EQ 3) AND ((flags AND 1) NE 0)))

  bandname = keyword_set(long) ? 'WideL' : 'WideS'
  print, 'statistics for band: ', bandname
  print, 'total number of sources in catalog :', n_elements(cat)
  print, 'FQUAL = 3 (source detected in this band) : ', good, ' occurrences'
  print, 'FQUAL = 3, non-detection flag :', nondet,' occurrences'
  print, 'FQUAL = 3, side-lobe affected flag : ', sidelobe,' occurrences'
  print, 'FQUAL = 3, too faint flag : ', toofaint, ' occurrences'
  print, 'FQUAL = 3, CDS mode flag : ', cdsmode, ' occurrences'

end
