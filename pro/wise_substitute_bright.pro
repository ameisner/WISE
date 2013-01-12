;+
; NAME:
;   wise_substitute_bright
;
; PURPOSE:
;   combine bright source PSF annulus fits with standard WISE catalog
;
; CALLING SEQUENCE:
;   wise_substitute_bright, catmag, catra, catdec, fluxflag, fpath=fpath
;
; INPUTS:
;   catmag   - WISE catalog w3mpro values, after catalog has been
;            filtered with wise_filter_catalog and sorted by dec
;   catra    - WISE catalog ra values, with same stipulations as catmag
;   catdec   - WISE catalog dec values, with same stipulations as catmag
;
; KEYWORDS:
;   fpath    - location in which .fits binary tables of PSF annulus fits
;               are stored
; OUTPUTS:
;   fluxflag - bitflag based on annulus fits, 1=source is bogus, 0=source OK
;
; EXAMPLES:
;   see wise_starlist.pro
;
; REVISION HISTORY:
;   2012-Oct-8 - Written by Aaron Meisner
;
;----------------------------------------------------------------------
function assemble_fit_results, fpath

  spawn, 'ls '+fpath+'/*.fits', flist
  nfile = n_elements(flist)

  res = mrdfits(flist[0], 1, /silent)
  for i=1,nfile-1 do res = struct_append(res, mrdfits(flist[i], 1, /silent))
  return, res

end

function analyze_fit_metrics, fitra, fitdec, fitmag, fpath=fpath, good=good, $
                              res=res

  if ~keyword_set(fpath) then fpath = '$WISE_DATA/fluxes'
  res = assemble_fit_results(fpath)
  nsrc = n_elements(res)

  fitmag  = dblarr(nsrc)
  fitra   = dblarr(nsrc)
  fitdec  = dblarr(nsrc)
  xyshift = dblarr(2,nsrc)
  chi2    = fltarr(nsrc)
  pixgood = lonarr(nsrc)

  hasfit = res.nfit NE 0
  for i=0, nsrc-1 do begin
      if ~hasfit[i] then continue
      wfit          = where(res[i].flag EQ 0)
      fitmag[i]     = median((res[i].mag)[wfit])
      fitra [i]     = median((res[i].radec)[0,wfit])
      fitdec[i]     = median((res[i].radec)[1,wfit])
      xyshift[*,i]  = median((res[i].xyshift)[*,wfit]) 
      chi2[i]       = median((res[i].chimed)[wfit])
      pixgood[i]    = median((res[i].pixstats)[2,wfit])
  endfor

  euler, res.catra, res.catdec, lgal, bgal, 1 
  plane     = abs(bgal) LT 2
  f         = 10^(-fitmag/2.5) ; variable linear in flux
  chi2upper = 5.+25.*plane+0.26*f+0.2*f^2 ; empirical
  disp      = sqrt(xyshift[0,*]^2+xyshift[1,*]^2) 
  good      = (hasfit) AND (chi2 LT chi2upper) AND (disp LT 3) AND $
              (pixgood GE 250)

  return, res.index

end

pro wise_substitute_bright, catmag, catra, catdec, fluxflag, w4=w4, fpath=fpath

; ----- have not yet fit custom fluxes for bright W4 sources
  if keyword_set(w4) then begin
      fluxflag = bytarr(n_elements(catra))
      return
  endif

  if ~keyword_set(fpath) then fpath='$WISE_DATA/fluxes'
  wfit = analyze_fit_metrics(fitra, fitdec, fitmag, fpath=fpath, good=good)

  wgood = where(good, ngood, complement=wbad)

  if (ngood GT 0) then begin 
      catmag[wfit[wgood]] = fitmag[wgood]
      catra[wfit[wgood]]  = fitra[wgood]
      catdec[wfit[wgood]] = fitdec[wgood]
  endif

  fluxflag       = bytarr(n_elements(catmag))
  fluxflag[wfit[wbad]] = 1

end
