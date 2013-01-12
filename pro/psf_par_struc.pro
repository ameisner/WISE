;+
; NAME:
;   psf_par_struc
;
; PURPOSE:
;   repository for PSF model related parameters
;
; CALLING SEQUENCE:
;   par = psf_par_struc(w4=, feat=, everything=, allsky=)
;   
; KEYWORDS:
;   w4   - set if W4 parameters desired rather than W3 parameters
;   feat - particular PSF feature of interest, if specified should be
;          one of: 'wings', 'core', 'ghost', 'latent'
;   everything - set to retrieve all possible parameters for band/release
;                of interest
;   allsky - set to retrieve all-sky release parameters
;
; OUTPUTS:
;   par  - structure containing PSF model parameters
;   
; EXAMPLES:
;   
; COMMENTS:
;
; REVISION HISTORY:
;   2011-Oct-14 - Aaron Meisner
;----------------------------------------------------------------------
function psf_par_struc, w4=w4, feat=feat, everything=everything, allsky=allsky

  forward_function psf_par_struc

  if keyword_set(feat) && (feat NE 'wings') AND (feat NE 'core') AND $ 
    (feat NE 'ghost') AND (feat NE 'latent') AND (feat NE 'latent2') then begin
      print, $ 
'PSF feature must be one of the following: wings, core, ghost, latent, latent2'
      return, psf_par_struc(w4=w4, everything=everything, allsky=allsky)
  endif

  w4 = keyword_set(w4)
; ----- sidelength of L1b image in pixels
  impix = w4 ? 508 : 1016

; ----- coordinate of image center
  crpix = w4 ? 253.5 : 507.5 

; ----- pixel scale, asec/pix
  pscl  = w4 ? 5.53  : 2.754

; ----- threshold for tolerable number of NaNs in an L1b image
  nsatmax = w4 ? 125000 : 500000

; ----- zero point
  magzp = w4 ? 12.90 : (keyword_set(allsky) ? 17.645 : 17.6)

; ----- file containing L1b exposure metadata table
  indexfile = w4 ? '$WISE_DATA/index-allsky-w4.fits' : $ 
                (keyword_set(allsky) ? '$WISE_DATA/index-allsky-L1b.fits' : $ 
                                       '$WISE_DATA/index-metadata-L1b.fits')
; ----- file containing compact source catalog
  catfile = w4 ?  $ 
               '$WISE_DATA/w4_catalog-allsky.fits' : $ 
            (keyword_set(allsky) ? $ 
               '$WISE_DATA/w3_catalog-allsky.fits' : $  
               '$WISE_DATA/w3_catalog.fits' )

; ----- file containing latent image
  latim = w4 ? '$WISE_DATA/latent-taper-w4.fits' : $ 
        (keyword_set(allsky) ? '$WISE_DATA/latent-nonlinear-allsky.v1.fits' : $
                               '$WISE_DATA/latent-nonlinear.fits')

; ----- file containing image of second latent
  latim2 = w4 ? '$WISE_DATA/latent2-w4.fits' : '$WISE_DATA/latent2-w3.fits'

; ----- file containing PSF polynomial coefficients
  fpsf = w4 ? '$WISE_DATA/psf_coeff-taper-w4.fits' : $ 
           (keyword_set(allsky) ? '$WISE_DATA/psf_coeff-big.fits' : $
                                  '$WISE_DATA/psf_coeff.fits')

; ----- file containing ghost image cutout, not yet applicable to W4
  fghost = w4 ? '' : '$WISE_DATA/ghost-w3-big.fits'

; ----- sidelength of square cutout over which PSF wings are modeled
  szwings = w4 ? 285 : 325

; ----- size of core region expected to be ruined by saturation in
;       in ultra-bright stars used to model PSF wings
  szcore = w4 ? 27 : 35 ; W3 value ?

; ----- latent model sidelength
  szlat = w4 ? 191 : 325

; ----- 2nd latent model sidelength
  szlat2 = w4 ? 189 : 271

; ----- faintest w?mpro mag that can induce latent
  latmax = w4 ? 1.2 : 4.0

; ----- faintest w?mpro mag for which to set mask bit for Nth latent
  lmskmax = w4 ? [1.2,-0.8,-2.8,-3.4] : [4.,1.5,0.0,-1.5]

; ----- brightest w?mpro mag for which latent nonlinearity modeled
  latmin = w4 ? -4.8 : -2.0

; ----- faintest w?mpro mag for which to apply second latent correction
  l2max = w4 ? -1.0 : -0.5

; ----- latent bitmask values for 1st, 2nd, 3rd, 4th latents
  latflag = [8, 256, 512, 1024]

; ----- ghost model x sidelength
  xgpix = w4 ? 45 : 135

; ----- ghost model y sidelength
  ygpix = w4 ? 33 : 87

; ----- y offset between center of PSF core model and center of ghost model
  ygoffs = w4 ? 104 : 206

; ----- size of full PSF model including core+wings+ghost, could be calculated
;       from szwings, ygpix, ygoffs but this will be convenient
  psfpix = w4 ? 285 : 499

; ----- size of PSF cutout to be subtracted from faint sources
  pfaint = 115 ; applies to both W3, W4

; ----- magnitudes by which to pad definition of bright versus faint
  bpad = w4 ? 0.00 : 0.15

; ----- w?mpro magnitude of faintest source for which sinc interpolation
;       performed before PSF subtraction
  interplim = w4 ? 8.25 : 11.5

; ----- order of per-pixel polynomial for each PSF model component
  corder = w4 ? 2 : 3
  worder = 1
  gorder = 3
  lorder = 0

; ----- radius inside of which PSF wings are replaced with values from
;       PSF core fit
  radcore = w4 ? 6 : 13 ; W3 value ?

; ----- sky annulus parameters for PSF wings
  skyrad = w4 ? szwings/2+23+[0,5] : szwings/2+[0, 5]

; ----- sky annulus parameters for latent
  skyrad_lat = szlat/2+[0, 5]

; ----- sky annulus parameters for second latent
  skyrad_lat2 = szlat2/2+[0, 5]

; ----- sidelength of cleaned L1b images in pixels
  pclean = w4 ? 246 : 500 

; ----- general L1b properties included in all cases
  par = { impix      : impix,      $
          crpix      : crpix,      $
          pscl       : pscl,       $ 
          nsatmax    : nsatmax,    $ 
          indexfile  : indexfile,  $
          pclean     : pclean,     $ 
          catfile    : catfile      }

  if ~keyword_set(feat) AND ~keyword_set(everything) then begin
      return, par
  endif

; ----- if you set the everything keyword, you're getting everything
;       and any value of the feat keyword will be ignored
  if keyword_set(everything) then begin
      add = { szwings     : szwings,     $
              szlat       : szlat,       $
              szlat2      : szlat2,      $
              szcore      : szcore,      $
              psfpix      : psfpix,      $    
              pfaint      : pfaint,      $
              xgpix       : xgpix,       $
              ygpix       : ygpix,       $
              ygoffs      : ygoffs,      $
              skyrad      : skyrad,      $ 
              skyrad_lat  : skyrad_lat,  $
              skyrad_lat2 : skyrad_lat2, $
              corder      : corder,      $ 
              worder      : worder,      $
              gorder      : gorder,      $ 
              lorder      : lorder,      $
              radcore     : radcore,     $ 
              magzp       : magzp,       $ 
              latmax      : latmax,      $
              latmin      : latmin,      $
              l2max       : l2max,       $
              lmskmax     : lmskmax,     $
              latflag     : latflag,     $
              latim       : latim,       $ 
              latim2      : latim2,      $
              fpsf        : fpsf,        $ 
              bpad        : bpad,        $ 
              interplim   : interplim,   $
              fghost      : fghost        }

      return, struct_addtags(par, add)
  endif

  if keyword_set(feat) then begin
      case feat of
          'wings' :   begin
                        add = { szx       : szwings,     $
                                szy       : szwings,     $
                                skyrad    : skyrad,      $
                                yoffs     : 0,           $
                                order     : worder        }
                      end
          'core'  :   begin
                        add = { szx       : szcore,      $
                                szy       : szcore,      $
                                skyrad    : skyrad,      $
                                yoffs     : 0,           $ 
                                order     : corder        }
                      end
          'ghost' :   begin
                        add = { szx       : xgpix,       $
                                szy       : ygpix,       $
                                skyrad    : skyrad,      $
                                yoffs     : ygoffs,      $ 
                                order     : gorder        }
                      end
          'latent':   begin
                        add = { szx       : szlat,       $
                                szy       : szlat,       $
                                skyrad    : skyrad_lat,  $
                                yoffs     : 0,           $ 
                                order     : lorder        }
                      end
          'latent2':  begin
                        add = { szx       : szlat2,       $
                                szy       : szlat2,       $
                                skyrad    : skyrad_lat2,  $
                                yoffs     : 0,           $ 
                                order     : lorder        }
                      end
      endcase
      return, struct_addtags(par, add)
  endif

end