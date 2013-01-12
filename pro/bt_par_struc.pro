;+
; NAME:
;   bt_par_struc
;
; PURPOSE:
;   act as common repository for important blue tip fitting parameters
;   
; CALLING SEQUENCE:
;   par = bt_par_struc(PS= )
;   
; KEYWORDS:
;   PS     - set if you want Pan-STARRs parameters, otherwise you'll
;            get SDSS parameters
;
; OUTPUTS:
;   par    - structure with blue tip fitting parameters for desired survey
;   
; REVISION HISTORY:
;   2012-Jul-17 - Written by Aaron Meisner
;----------------------------------------------------------------------
function bt_par_struc, PS=PS

  bands   = keyword_set(PS) ? ['g','r','i','z','y'] : ['u','g','r','i','z']
; ----- reddening coefficients, normalized to A_r
;                       u        g      r       i         z        y
  acoeff  = keyword_set(PS) ? [1.39674, 1., 0.740643, 0.582122, 0.478644] : $ 
                     [1.85514, 1.44551, 1., 0.743107, 0.552735]
; ----- blue floor values, still need to fit value for PS zy
;                        ug      gr     ri     iz    zy
  F       = keyword_set(PS) ? [0.031, 0.038, 0.027, 0.039] : $
                       [0.220, 0.050, 0.050, 0.063]
; ----- fitting region lower, upper boundaries relative to estimate of x0
  dxlower = keyword_set(PS) ? [0.180, 0.170, 0.170, 0.170] : $ 
                       [0.150, 0.130, 0.170, 0.170] 
  dxupper = keyword_set(PS) ? [0.370, 0.230, 0.230, 0.230]: $ 
                       [0.650, 0.420, 0.230, 0.230]
; ----- typical low-reddening values of bluetip width parameter sigma_P
  sigma   = keyword_set(PS) ? [0.050, 0.025, 0.030, 0.030] : $ 
                       [0.050, 0.050, 0.025, 0.030]
; ----- approximate zero-reddening x0 value, the location of blue tip
  x0      = keyword_set(PS) ? [0.230, 0.040, -0.04, -0.055] : $ 
                       [0.800, 0.230, 0.070, -0.02]
; ----- catalog marker for non-detection, 0 for PS, 30 for SDSS
  nondet  = keyword_set(PS) ? 0. : 30.
; ----- Schlafly & Finkbeiner 2011 conversion factor from SFD E(B-V) to A_r
  ebv2ar  = keyword_set(PS) ? 2.271 : 2.285

  par = { bands      : bands,   $
          acoeff     : acoeff,  $
          F          : F,       $
          dxlower    : dxlower, $
          dxupper    : dxupper, $
          sigma      : sigma,   $
          x0         : x0,      $
          nondet     : nondet,  $
          ebv2ar     : ebv2ar     }

  return, par

end
