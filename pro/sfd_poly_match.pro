; NAME:
;   sfd_poly_match
;
; PURPOSE:
;   compute first order warp of WISE L1b image that best matches SFD i100
;
; CALLING SEQUENCE:
;   sfd_poly_match, im, msk, h, gain, offs, grad, nrej, r=r
;
; INPUTS:
;   im    - cleaned L1b image
;   msk   - artifact bitmask corresponding to cleaned L1b image
;   h     - header of cleaned L1b image
;
; OUTPUTS:
;   gain  - best fit "gain" in least squares fit, gain here
;           meaning conversion factor from i100 to W3 DN
;   offs  - best fit constant offset in least squares fit
;   grad  - best fit gradients, two-element array [xgrad, ygrad],
;           units DN/pix appropriate for 500x500 input image
;   nrej  - number of pixels rejected during iterative least squares
;   
; OPTIONAL OUTPUTS:
;   r     - Pearson R correlation coefficient, -1 <= r <= 1
;
; EXAMPLES:
;   see batch_poly_match below
;
; COMMENTS:
;   could definitely put more thought into best possible use of
;   bitmasks and the optimal amount by which to bin down before smoothing
;
;   warning: input im is modified (should probably change this...)
;   
; REVISION HISTORY:
;   2012-Oct-8 - Aaron Meisner
;
;----------------------------------------------------------------------
pro sfd_poly_match, im, msk, h, gain, offs, grad, nrej, r=r, w4=w4

; ----- dilate msk aggressively and then interpolate
  lpix = keyword_set(w4) ? 21 : 41
  latentkern = shift(dist(lpix), lpix/2, lpix/2) LT (lpix/2+0.5)
  spix = keyword_set(w4) ? 11 : 23
  satkern = shift(dist(spix), spix/2, spix/2) LT (spix/2+0.5)
  gpix = keyword_set(w4) ? 15 : 31
  ghostkern = shift(dist(gpix), gpix/2, gpix/2) LT (gpix/2+0.5)
  mask_dilate = bytarr(size(im, /DIM))
  mask_dilate OR= dilate((msk AND 8) NE 0, latentkern) ; latent = 2^3
  mask_dilate OR= dilate(msk AND 1, satkern) ; sat = 2^0
  mask_dilate OR= dilate((msk AND 2) NE 0, ghostkern) ; ghost = 2^1
  intx = djs_maskinterp(im, mask_dilate, iaxis=0)
  inty = djs_maskinterp(im, mask_dilate, iaxis=1)
  im = (intx+inty)/2

; ----- median filter
  im = median(im, 5)
; ----- bin down a lot, to ~30 asec pixels, assume 500x500 to start
;       for W3, 246x246 to start for W4
  par = psf_par_struc(w4=w4)
  pixclean = par.pclean
  binfac = keyword_set(w4) ? 3 : 5
  pixsmall = keyword_set(w4) ? 82 : 100
  npsmall = long(pixsmall)*pixsmall
  im = rebin(im, pixsmall, pixsmall)
; ----- smooth to SFD resolution, pixels now 2.754*2*5 = 27.54 asec (W3)
;       and 5.53*2*3 = 33.18 asec (W4)
;       13.2521 = sqrt((6.1*60)^2-(2.754*2*5)^2)/27.54
;       10.9853 = sqrt((6.1*60)^2-(5.53*2*3)^2)/33.18
  kpix = keyword_set(w4) ? 55 : 67 ; 5+ fwhm
  fwhm = keyword_set(w4) ?  10.9853 : 13.2521
  kern = psf_gaussian(NPIXEL=kpix, FWHM=[fwhm, fwhm], /NORM)
  smth = convol(im, kern, /edge_zero)
  wt = convol(fltarr(pixsmall, pixsmall)+1., kern, /edge_zero)
  smth = smth/wt ; wt should never be zero
; ----- get SFD on same footprint
  hrebin, fltarr(pixclean,pixclean, /nozero), h, _, hsmall, $ 
      OUTSIZE=[pixsmall,pixsmall]
  extast, hsmall, astr
  xbox = lindgen(pixsmall, pixsmall) MOD pixsmall
  ybox = lindgen(pixsmall, pixsmall) / pixsmall
  xy2ad, xbox, ybox, astr, ra, dec
  euler, ra, dec, lgal, bgal, 1
  i100 = dust_getval(lgal, bgal, map='i100', /interp, /noloop)
; ----- do hogg_iter_linfit, where "target" values are WISE pixel values
;       and model is gain*sfd+offset+xgrad*dx+ygrad*dy
  crpixsm = keyword_set(w4) ? 40.5  : 49.5
  nparam = 4
  wfit = lindgen(npsmall)
  A = fltarr(nparam, npsmall)
  A[0, *] = i100[wfit] ; for gain term
  A[1, *] = 1 ; for constant offset term
  A[2, *] = xbox[wfit]-crpixsm
  A[3, *] = ybox[wfit]-crpixsm
  YY = smth[wfit]
  iweight = fltarr(npsmall)+1.
  hogg_iter_linfit, A, YY, iweight, XX

; ---- make an image of what was kept/rejected for debugging
  iweight = reform(iweight, pixsmall, pixsmall)

  gain = XX[0]
  offs = XX[1]
  grad = [XX[2], XX[3]]/float(binfac) ; put back into cleaned L1b pixel scale
  nrej = long(total(iweight EQ 0))

; ----- construct the model
  synth = XX[0]*i100+XX[1]+(xbox-crpixsm)*XX[2]+(ybox-crpixsm)*XX[3]
; ----- compute correlation coefficient
  if arg_present(r) then begin
      if (nrej NE npsmall) then $
          r = correlate(synth[where(iweight NE 0)],smth[where(iweight NE 0)]) $
      else $ 
          r = -9999 ; if r is fit it is between +/- 1 by definition
  endif

end

function package_results, scan_id, frame_num, gain, offs, grad, nrej, r

  outstr = { scan_id   : scan_id,   $ 
             frame_num : frame_num, $
             gain      : gain,      $
             offs      : offs,      $
             grad      : grad,      $ 
             nrej      : nrej, $ 
             r         : r           }
  return, outstr

end

pro batch_poly_match, indstart, nproc, indstr, result=result, $ 
                      outpath=outpath, w4=w4, cleanpath=cleanpath

  write = ~arg_present(result)
  if ~keyword_set(outpath) then outpath='$WISE_DATA/warps'
  if ~keyword_set(cleanpath) then cleanpath = keyword_set(w4) ? $ 
    '/n/panlfs/ameisner/clean.w4.v0' : $ 
    '/n/wise/ameisner/clean.allsky.v1'
  print, 'cleanpath = ', cleanpath
  nfile = n_elements(indstr)
  indend = (long(indstart)+nproc-1) < (nfile-1)
  for i=long(indstart), indend do begin
     if (i MOD 50) EQ 0 then print, i
     cname = wise_name('clean', indstr[i].fname, $ 
         cleanpath=cleanpath)
     im  = lfs_fits_access(cname, h, /silent)
     if n_elements(im) EQ 1 then continue
     msk = lfs_fits_access(cname, exten=2, /silent)
     t0 = systime(1)
     sfd_poly_match, im, msk, h, gain, offs, grad, nrej, r=r, w4=w4
     dt = systime(1)-t0
     if (i MOD 50) EQ 0 then print, 'dt = ', dt, ' s'
     outstr = package_results(indstr[i].scan_id, $ 
         indstr[i].frame_num, gain, offs, grad, nrej, r)
     delvarx, r
     if n_elements(result) EQ 0 then result=outstr else $
         result = struct_append(result, outstr)
  endfor

  if write then begin
      outname = 'warp_'+string(indstart, format='(I07)')+'.fits'
      outname = concat_dir(outpath, outname)
      mwrfits, result, outname
  endif
 
end

pro warp_one_tile, tnum, indstart, nproc, outpath=outpath, w4=w4

  tstr = mrdfits('$WISE_DATA/wisetile-index-allsky.fits', 1)
  tind = where(tnum EQ fix(strmid(tstr.fname,5,3))) ; assume valid tnum
  euler, tstr[tind].ra, tstr[tind].dec, lgal, bgal, 1
  lb = [lgal,bgal]
  indstr = wise_index_metadata(lb, angle=9.0, /allsky, w4=w4)
  
  batch_poly_match, indstart, nproc, indstr, outpath=outpath, w4=w4

end

pro warp_by_mjd, indstart, nproc, outpath=outpath, w4=w4, cleanpath=cleanpath

   if ~keyword_set(outpath) then outpath = '$WISE_DATA/warps'
   indstr = wise_index_metadata([0.,0.], /nosort, /allsky, w4=w4)
   batch_poly_match, indstart, nproc, indstr, outpath=outpath, w4=w4, $ 
       cleanpath=cleanpath

end

function assemble_poly_fits, fpath=fpath

; routine to make one structure out of the multiple .fits outputs from
; above batch processing routines

  if ~keyword_set(fpath) then fpath = '$WISE_DATA/warps'

  spawn, 'ls '+fpath+'/warp_???????.fits', flist
  if ~file_test(flist[0]) then return, -1

  nfile = n_elements(flist)
  print, 'reading polynomial warps from ', nfile, ' files'
  res = mrdfits(flist[0], 1, /silent)
  for i=1, nfile-1 do res = struct_append(res, mrdfits(flist[i], 1, /silent))

  return, res

end

pro make_healpix_maps, gain_map, r_map, grad_map, rej_map, nside=nside, ecl=ecl

; reading in the .fits summary files is by far the slowest part of this...

  if ~keyword_set(nside) then nside = 32
  nptot = 12L*nside*nside

  res    = assemble_poly_fits()
  indstr = wise_index_metadata([0.,0.], /nosort, /allsky)

  id_all = indstr.scan_id + string(indstr.frame_num, format='(I03)')
  id     = res.scan_id + string(res.frame_num, format='(I03)')

  matchlist, id, id_all, ind, ind_all

  nwarp = n_elements(res)
  ra = dblarr(nwarp)
  dec = dblarr(nwarp)
  ra[ind] = indstr[ind_all].ra
  dec[ind] = indstr[ind_all].dec

  select = keyword_set(ecl) ? 3 : 1
  euler, ra, dec, l, b, select

  ang2pix_ring, nside, (90.-b)/!radeg, l/!radeg, pix

  pix_unique = pix[uniq(pix, sort(pix))]
  nunique = n_elements(pix_unique)

  print, 'nunique = ', nunique

  matchlist, pix_unique, pix, ind_pix, ind_res
  sind = sort(ind_pix)
  ind_pix = ind_pix[sind]
  ind_res = ind_res[sind]
  bdy = uniq(ind_pix)

  gain_map = fltarr(nptot)
  r_map = fltarr(nptot)
  grad_map = fltarr(nptot)
  rej_map = fltarr(nptot)
  for i=0L, nunique-1 do begin
      if (i MOD 200) EQ 0 then print, i
      indlower = (i EQ 0) ? 0 : bdy[i-1] + 1
      indupper = bdy[i]
      wthispix = lindgen(indupper-indlower+1) + indlower
; ---- take median by pixel
      gain_map[pix_unique[i]] = median(res[ind_res[wthispix]].gain)
      r_map[pix_unique[i]] = median(res[ind_res[wthispix]].r)
      rej_map[pix_unique[i]] = $ 
          float(median(res[ind_res[wthispix]].nrej))/10000. ; wrong for W4!!
      thisgrad = res[ind_res[wthispix]].grad
      grad_map[pix_unique[i]] = median(sqrt(thisgrad[0,*]^2+thisgrad[1,*]^2))
  endfor

end

pro write_healpix_maps, nside=nside, ecl=ecl

    make_healpix_maps, gain_map, r_map, grad_map, rej_map, nside=nside, ecl=ecl
   
    pstr = keyword_set(ecl) ? 'ecl' : 'gal'
    outpath = '$WISE_DATA/warps/healpix'
    outname = 'warp_'+strtrim(string(nside),1)+'.'+pstr+'.fits'
    outname = concat_dir(outpath, outname)

    JyperDN = 2.9045E-06 ; [Jy per DN] wise prelim data release docs Sec 2.3f 
    calfac = JyperDN/((2.75/3600./!radeg)^2)*1E-6 ; [MJy/sr/DN]
    pix_deg = 3600./(2.754*2) ; pixels per deg
    grad_DN_deg = grad_map*pix_deg ; gadient in DN
    grad_MJysr_deg = grad_DN_deg*calfac ; gradient in MJy/sr per deg
    rat = gain_map*calfac ; 100um to 12um ratio when both in same units

    writefits, outname, rat
    writefits, outname, grad_MJysr_deg, /append
    writefits, outname, r_map, /append
    writefits, outname, rej_map, /append

end

pro index_append_grad, indstr, w4=w4

  if keyword_set(w4) then fpath = '/nfs_pan1/ameisner/warps'
  res    = assemble_poly_fits(fpath=fpath)
  par = psf_par_struc(w4=w4, /allsky)
  indstr = mrdfits(par.indexfile, 1)
  nfile  = n_elements(indstr)

  id_all = indstr.scan_id + string(indstr.frame_num, format='(I03)')
  id     = res.scan_id + string(res.frame_num, format='(I03)')

  matchlist, id, id_all, ind, ind_all

  grad_out = dblarr(2, nfile)
  grad_out[*, ind_all] = res[ind].grad

  gradstr = replicate({grad:[0.d,0.d]}, nfile)
  gradstr.grad = grad_out

  indstr = struct_addtags(indstr,gradstr)

end
