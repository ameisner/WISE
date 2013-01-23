;+
; NAME:
;   previous exposure
;
; PURPOSE:
;  find the exposures N*11.1s prior to specified exposure (1 <= N <= 5)
;   
; CALLING SEQUENCE:
;   fprev = previous_exposure(mjd, allsky=, w4=, forward=, N=)
;
; INPUTS:
;   mjd     - MJD of exposure of interest
;
; KEYWORDS:
;   allsky  - set for all-sky release
;   w4      - set for W4 (default is W3)
;   forward - set to retrieve exposures FOLLOWING input mjd (forward in time)
;   N       - maximum range of time for which to return file names, in units
;            of WISE exposure cycles (1 cycle = 11.1s)
;
; OUTPUTS:
;   outstr - structure with field fname containing names of 
;            previous exposure(s) and corresponding field N giving
;            number of 11.1s time intervals by which each file name 
;            preceded/followed specified MJD
;
; EXAMPLES:
;   see wise_correct_latent.pro
;
; COMMENTS:
;   file names in output structure are always time-ordered
;
;   the intended usage of this routine is that input mjd be
;   an actual WISE exposure time, not some arbitrary value
;
; REVISION HISTORY:
;   2011-Oct-18 - Aaron Meisner
;----------------------------------------------------------------------
function previous_exposure, mjd, allsky=allsky, w4=w4, forward=forward, N=N, $ 
                            band=band

; ----- N=0 is kind of meaningless here so don't do anything in that case
  if (size(N, /TYPE) NE 0) && ((abs(N) GT 5) OR (N EQ 0)) then begin
      print, 'previous_exposure.pro only verified for N=1, 2, 3, 4, 5'
      return, -1
  endif

  if ~keyword_set(N) then N = 1
  forward = keyword_set(forward)
; ----- preferably forward vs. backward in time should be specified with
;       forward keyword, rather than sign of N but handle this anyway...
  if (N LT 0) then begin
      N = -N
      forward = ~forward
  endif
  dt = 11.1 ; length of one WISE exposure cycle, seconds

  par = psf_par_struc(allsky=allsky, w4=w4, band=band)
; ----- assume index files sorted by MJD !!!
  indexfile = par.indexfile

  COMMON MJD, fname_sorted, mjd_sorted, nexp, indexfile_cache
  if (n_elements(mjd_sorted) EQ 0) || (indexfile NE indexfile_cache) then begin
    indstr = mrdfits(indexfile, 1)
    fname_sorted = indstr.fname
    mjd_sorted = indstr.mjd
    nexp = n_elements(mjd_sorted)
    indexfile_cache = indexfile
  endif

  binary_search, mjd_sorted, mjd + double(1.0e-5), thisind
  if (thisind EQ -1) then thisind = nexp - 1
  ind_u = (thisind + (forward ? N : 0)) < (nexp-1)
  ind_l = (thisind - (forward ? 0 : N)) > 0

  if ((ind_u-ind_l) EQ 0) then return, -1

; ----- this next step works because, without exception, the time lags
;       between exposures fall within these narrow, discrete ranges:
;    N = 1 : 10.713600 s < dt < 11.145600 s
;    N = 2 : 21.513600 s < dt < 22.291200 s
;    N = 3 : 32.227200 s < dt < 33.436800 s 
;    N = 4:  43.027200 s < dt < 44.496000 s
;    N = 5:  53.827200 s < dt < 55.641600 s
  dN = $ 
round(86400.d*(forward ? 1 : -1)*(mjd_sorted[(ind_l+forward):(ind_u-(~forward))]-mjd)/dt)

  w = where((dN GT 0) AND (dN LE N), nw)
  if (nw EQ 0) then return, -1
  outstr = replicate({fname:'',N:0B}, nw)
  outstr.N = dN[w]
  outstr.fname = (fname_sorted[(ind_l+forward):(ind_u-(~forward))])[w]

  return, outstr

end
