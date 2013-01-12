;+
; NAME:
;   wise_status
;
; PURPOSE:
;   Report on processing status of clean images and ISSA tiles
;
; CALLING SEQUENCE:
;   wise_status, cleanpath=cleanpath, tilepath=tilepath
;
; INPUTS:
;   cleanpath   - path to clean file tree
;   tilepath    - path to ISSA tiles
; 
; OUTPUTS:
;   plot
;
; COMMENTS:
;   Uses xargs for speed. 
;
; REVISION HISTORY:
;   2012-Feb-26 - Written by Douglas Finkbeiner, CfA
;   2012-Feb-29 - Added plot vs. time - DPF
;
;----------------------------------------------------------------------
pro wise_status, cleanpath=cleanpath, tilepath=tilepath, allsky=allsky, w4=w4
  
  if ~ keyword_set(cleanpath) then cleanpath = '/n/panlfs/ameisner/clean.prelim.v0'
  if ~ keyword_set(tilepath)  then tilepath = '~dfink/runwise/tile'

  window, 0, xsize=600, ysize=900

  splog, 'Getting file list...'
  pushd, cleanpath
  spawn,'\ls '+cleanpath+' | xargs -n 1 \ls', flist
  popd

  splog, 'Reading metadata table...'
  indstr = wise_index_metadata([0., 0.], /NOSORT, allsky=allsky, w4=w4)
  fname = fileandpath(indstr.fname)

  splog, 'Matching filename lists'
  matchlist, fname, flist, m1, m2

  euler, indstr.ra, indstr.dec, l, b, 1

  rand = (randomu(!pi, n_elements(l)) LT 0.1) 
  rand[m1] = 0
  w = where(rand, nw)

  pmulti = !p.multi
  !p.multi = [0, 1, 2]
  if nw GT 0 then begin 
     djs_plot, l[w], b[w], ps=3, xr=[360, 0], yr=[-90, 90], /xst, /yst, color='gray', xtit='l gal [deg]', ytit='b gal [deg]'
  endif else begin
     djs_plot, [l[0]], [b[0]], ps=3, xr=[360, 0], yr=[-90, 90], /xst, /yst, color='gray', xtit='l gal [deg]', ytit='b gal [deg]', /nodata
  endelse     
  djs_oplot, l[m1], b[m1], ps=3, xr=[360, 0], yr=[-90, 90], /xst, /yst, color=''

  print, n_elements(m1), ' completed'

; -------- tiles
  tlist = file_search(concat_dir(tilepath, 'wise_???.fits'), count=ntile)
  tpos = strpos(tlist, '.fits')
  tstr = strmid(tlist, tpos[0]-3, 3)
  for i=0L, ntile-1 do begin 
     h = headfits(tlist[i])
     euler, sxpar(h, 'CRVAL1'), sxpar(h, 'CRVAL2'), l, b, 1
     djs_oplot, [l], [b], ps=7, color='green'
     djs_xyouts, [l], [b], tstr[i], color='blue', align=0.5
  endfor

; -------- plot 2
  pra = 269.9   ; ra, dec of survey poles
  pdec = 66.55
  theta = djs_diff_angle(indstr.ra, indstr.dec, pra, pdec)

  mjd_l = floor(min(indstr.mjd))
  mjdrange = ceil(max(indstr.mjd)-mjd_l)
  xtit = 'MJD-'+string(mjd_l, format='(I05)')
  pad = 1

  djs_plot, indstr[m1].mjd-mjd_l, theta[m1], ps=3, color='', xtit=xtit, $ 
      ytit='survey pole angle', yr=[-5, 185], /yst, xr=[-pad, mjdrange+pad], $ 
      /xst
  if nw GT 1 then djs_oplot, indstr[w].mjd-mjd_l, theta[w], ps=3, color='gray'

  !p.multi = pmulti

  print, ntile, ' tiles finished at ', systime()

  return
end
