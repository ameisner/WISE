;+
; NAME:
;   wise_check_atlas_files
;
; PURPOSE:
;   Verify that all WISE atlas files and directories are present
; 
; CALLING SEQUENCE:
;   wise_check_atlas_files
;   
; OUTPUTS:
;   messages to stdout
;
; EXAMPLES:
;   just run it. 
;
; COMMENTS:
;   This procedure identifies missing files and generates a query file
;    in IPAC table format, ready to upload to the IRSA web form. 
;
; REVISION HISTORY:
;   2011-Jun-6 - Written by Douglas Finkbeiner, CfA
;
;----------------------------------------------------------------------
pro wise_check_atlas_files

  tilefile = '$WISE_DIR/etc/wise_tiles.dat'
  readcol, tilefile, ii, aa, ra, dec, format='(L,A,D,D)'
  ntile = n_elements(ra) 

  path = '/n/panlfs/dfink/wisedata/L3a'
  badra = fltarr(1000)
  baddec = fltarr(1000)
  k = 0

  for i=0L, ntile-1 do begin 
     badtile = 0
     if (i mod 250) eq 0 then print, i
     tiledir = aa[i]+'_aa11'
     thisdir = concat_dir(path, tiledir)
     if file_test(thisdir) then begin 
        for j=0L, 3 do begin 
           bandstr = (['1', '2', '3', '4'])[j]
           intname = concat_dir(thisdir, aa[i]+'_aa11-w'+bandstr+'-int-3.fits')
           uncname = concat_dir(thisdir, aa[i]+'_aa11-w'+bandstr+'-unc-3.fits')
           covname = concat_dir(thisdir, aa[i]+'_aa11-w'+bandstr+'-cov-3.fits')
           if ~ file_test(intname) then badtile++
           if ~ file_test(uncname) then badtile++
           if ~ file_test(covname) then badtile++
        endfor
     endif else begin 
        print, 'Missing directory: ', tiledir
        badtile++
     endelse
     if badtile NE 0 then begin
        print, i, badtile, '   ', aa[i], ra[i], dec[i]
        badra[k] = ra[i]
        baddec[k] = dec[i]
        k++
     endif 
  endfor
  
  badra = badra[0:k-1]
  baddec = baddec[0:k-1]

; -------- write out a new query file to fill in missing data

  head1 = '|    ra      |   dec       |'
  head2 = '|     r      |     r       |'

  outfile = 'wisequery_fix.dat'
  outname = concat_dir('.', outfile)
  print, 'Writing ', outname
  openw, wlun, outname, /get_lun
  printf, wlun, head1
  printf, wlun, head2
  for i = 0L, k-1 do $
    printf, wlun, badra[i], baddec[i], format='(2F13.8)'
  free_lun, wlun
  
  return
end
