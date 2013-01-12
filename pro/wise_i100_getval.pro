;+
; NAME:
;   wise_i100_getval
;
; PURPOSE:
;   use WISE-enhanced SFD i100 tiles to quote i100 values for a
;   specified list of coordinates
;
; CALLING SEQUENCE:
;   i100 = wise_i100_getval(ra, dec, sfd=sfd)
;
; INPUTS:
;   ra    - right ascension
;   dec   - declination
;
; OPTIONAL INPUTS:
;   
; KEYWORDS:
;   sfd   - keyword to pick up SFD i100 values corresponding to input coords
;   flags - keyword to retrieve star artifact bitmask
;   tnum  - keyword to retrieve tile number from which each
;           coordinate pair was sampled
;
;
; OUTPUTS:
;   i100  - WISE-enhanced i100 values (MJy/sr) corresponding to input coords
;
; OPTIONAL OUTPUTS:
;   
; EXAMPLES:
;   
; COMMENTS:
;   flag bits inherit same definitions as given in the header of 
;   wise_l1b_clean.pro
;   
; REVISION HISTORY:
;   2012-Aug-1 - Written by Aaron Meisner
;
;----------------------------------------------------------------------
function wise_i100_getval, ra, dec, sfd=sfd, flags=flags, tnum=tnum

  tfile = '$WISE_DATA/wisetile-index-allsky.fits'
  tstr = mrdfits(tfile, 1)  
  ncoord = n_elements(ra)

  matchlength = 12.5/sqrt(2) ; deg
  _ = djs_angle_match(ra, dec, tstr.ra, tstr.dec, dtheta=matchlength, $ 
      mindx=tindx, mmax=1)
 
  sind = sort(tindx)
  tindx = tindx[sind]
  bdy = uniq(tindx)

  ntile = n_elements(bdy)
  i100_out = fltarr(ncoord)
  tnum = intarr(ncoord)
  flags = intarr(ncoord)
  for i=0, ntile-1 do begin
      indlower = (i EQ 0) ? 0 : bdy[i-1] + 1
      indupper = bdy[i]
      thisra = ra[sind[indlower:indupper]]
      thisdec = dec[sind[indlower:indupper]]
      _ = i100_to_ebv(tindx[indlower]+1, /noneg, astr=astr_tile, $ 
                      i100=tile_i100, bitmask=bitmask)
      ad2xy, thisra, thisdec, astr_tile, x, y
      i100 = interpolate(tile_i100, x, y)
      flags[sind[indlower:indupper]] = bitmask[round(x), round(y)]
      i100_out[sind[indlower:indupper]] = i100
      tnum[sind[indlower:indupper]] = tindx[indlower]+1
  endfor

  if arg_present(sfd) then begin
      euler, ra, dec, lgal, bgal, 1
      sfd = dust_getval(lgal, bgal, map='i100', /interp, /noloop)
  endif

  return, i100_out

end

pro read_fibers, ra, dec, number

  datafile = '/n/home09/ameisner/wise/pro/fiber_coords2.dat.txt'
  readcol, datafile, ra, dec, number, F='F, F, L'

end

pro dgl_lookup

  read_fibers, ra, dec, number
  i100 = wise_i100_getval(ra, dec, sfd=sfd, flags=flags, tnum=tnum)

  outstr = replicate({ra:0., dec:0., number:0L, wise:0., sfd:0., flag:0}, $ 
                     n_elements(ra))
  outstr.ra = ra
  outstr.dec = dec
  outstr.number = number
  outstr.wise = i100
  outstr.sfd = sfd
  outstr.flag = flags

  mwrfits, outstr, '/n/panlfs/ameisner/dgl/dgl.fits'

end

pro sspp_lookup, sfd=sfd, i100=i100

  fspec = '/n/panlfs/ameisner/peek.fits'
  specstr = mrdfits(fspec, 1)
  nspec = n_elements(specstr)

  wgood = where(specstr.goodflag)

  ra = specstr[wgood].ra
  dec = specstr[wgood].dec
  i100 = wise_i100_getval(ra, dec, sfd=sfd, $ 
                          flags=flags, tnum=tnum)

  addstr = replicate({wise:0.}, nspec)
  addstr[wgood].wise = i100
  outstr = struct_addtags(specstr, addstr)
  mwrfits, outstr, '/n/panlfs/ameisner/dgl/sspp/sspp.fits'

end
