;+
; NAME:
;   wise_build_hdr
;
; PURPOSE:
;   gather and store astrometry information from all headers
;
; CALLING SEQUENCE:
;   wise_build_hdr, START, NFILE
;
; INPUTS:
;   START  - starting index in L1b index structure
;   NFILE  - number of files for which to gather header information
; 
; EXAMPLES:
;   wise_build_hdr, 400000, 100000
; COMMENTS:
;   writes out structure containing a string array per L1b file
;   
; REVISION HISTORY:
;   2011-Mar-17 - Written by Aaron Meisner
;
;----------------------------------------------------------------------
pro wise_build_hdr, START, NFILE, allsky=allsky

  indstr = keyword_set(allsky) ? $ 
    mrdfits('$WISE_DATA/index-allsky-L1b.fits', 1) : $ 
    mrdfits('$WISE_DATA/index-metadata-L1b.fits', 1)

  fname  = indstr.fname
  NTOT   = n_elements(fname)

  LAST    = (START + NFILE - 1) < (NTOT - 1)

  for i  = long(START), long(LAST) do begin
    
    hastr = match_astr_params(fname[i])
    if (n_elements(outstr) EQ 0) then begin
      
      outstr = replicate({header:hastr}, LAST - START + 1)
      
    endif else begin
      outstr[i - START].header = hastr
    endelse
    if ((i MOD 100) EQ 0) then print, i

  endfor
  outname = 'hdr_' + strtrim(START, 1) + '.fits'
  mwrfits, outstr, outname

end

pro indstr_append, indstr_new

  hdr_file_dir = '/n/home09/ameisner/wise/pro' ; not anymore
  hdr_files = ['hdr_0.fits', 'hdr_100000.fits', 'hdr_200000.fits', $ 
               'hdr_300000.fits', 'hdr_400000.fits', 'hdr_500000.fits', $ 
               'hdr_600000.fits', 'hdr_700000.fits']
  
  nfile = n_elements(hdr_files)
  
  for i = 0, nfile - 1 do begin

    thisfile = concat_dir(hdr_file_dir, hdr_files[i])
    thisdata = mrdfits(thisfile, 1)
    if (n_elements(outstr) EQ 0) then begin
      outstr = thisdata
    endif else begin
      outstr = struct_append(outstr, thisdata)
    endelse

  endfor

  indstr = mrdfits('$WISE_DATA/index-metadata-L1b.fits', 1)

;------ check if index structure already has a header tag, if not
;       append it now
  if (tag_indx(indstr, 'HEADER') EQ -1) then begin

    indstr_new = struct_addtags(indstr, outstr)
;------ write file out somewhere random w/ space, can be moved around later
    mwrfits, indstr_new, '/n/panlfs/ameisner/index-meta-astr-L1b.fits'
  endif

end
