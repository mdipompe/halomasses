;+
;  NAME:
;    bias2mhalo
;  PURPOSE:
;    Convert bias to halo mass
;
;  USE:
;    mhalo=bias2mhalo(bias,z,bias_err=bias_err,weights=weights,h0=h0,omega_m=...,model='tinker10')
;
;  INPUT:
;    bias - single bias value
;    redshift - redshift(s).  Weighted equally if weights not supplied
;
;  OPTIONAL INPUT:
;    bias_err - error on the bias.  Assumed symmetric
;    weights - redshift weighing.  Should be array of same size as redshifts
;    h0 - little h (H_0/100). Defaults to 0.702
;    omega_m - omega_matter.  Defaults to 0.275
;    omega_l - omega_lambda. Defaults to 0.725
;    omega_b - omega_baryon. Defaults to 0.046
;    sigma_8 - power spectrum normalization sigma8.  Defaults to 0.81
;    spec_ind - power spectrum spectral index.  Defaults to 0.968
;    model - halo collapse model.  Options are:
;            'tinker10' (Tinker et al. 2010)  Default
;            'tinker05' (Tinker et al. 2005)
;            'smt' (Sheth, Mo, & Turman 2001)
;
;  OUTPUT:
;    Array containing mhalo, upper mhalo error, lower mhalo error.
;
;  NOTES:
;    Hacked from codes by Ryan Hickox
;
;  HISTORY:
;    8-26-15 - Written - MAD (UWyo)
;-
FUNCTION bias2mhalo, bias, redshift, bias_err=bias_err, $
                     h0=h0, omega_m=omega_m, omega_b=omega_b, $
                     omega_l=omega_l, sigma_8=sigma_8, $
                     model=model,weights=weights

;MAD Set default cosmology
IF ~keyword_set(h0) THEN h0=0.702
IF ~keyword_set(omega_m) THEN omega_m=0.275
IF ~keyword_set(omega_b) THEN omega_b=0.046
IF ~keyword_set(omega_l) THEN omega_l=0.725
IF ~keyword_set(sigma_8) THEN sigma_8=0.81
IF ~keyword_set(spec_ind) THEN spec_ind=0.968

;MAD Default model is Tinker et al. (2010)
IF ~keyword_set(model) THEN model='tinker10'

;MAD Set errors to 0 if needed
IF ~keyword_set(bias_err) THEN bias_err=0.

;MAD Generate array of log(mhalo) to test/interpolate
;MAD (from 6 to 16, where sigma(M) best estimated)
masses=6.+findgen(10000)*(0.001)

;MAD Loop over z
FOR i=0L,n_elements(redshift)-1 DO BEGIN
   ;MAD Get bias for each test Mhalo
   biases=mhalo2bias(masses, redshift[i],$
                     h0=h0, omega_m=omega_m, $
                     omega_b=omega_b, omega_l=omega_l, $
                     sigma_8=sigma_8, spec_ind=spec_ind, $
                     model=model)
   
   ;MAD Interpolate to desired bias
   quadterp, biases, masses, [bias-bias_err,bias,bias+bias_err],mhalo
   IF (n_elements(mhalo_errlo) EQ 0) THEN mhalo_errlo=mhalo[1]-mhalo[0] ELSE $
      mhalo_errlo=[mhalo_errlo,mhalo[1]-mhalo[0]]
   IF (n_elements(mhalo_errhi) EQ 0) THEN mhalo_errhi=mhalo[2]-mhalo[1] ELSE $
      mhalo_errlo=[mhalo_errhi,mhalo[2]-mhalo[1]]
   IF (n_elements(mhalo_fin) EQ 0) THEN mhalo_fin=mhalo[1] ELSE $
      mhalo_fin=[mhalo_fin,mhalo[1]]
ENDFOR

IF keyword_set(weights) THEN BEGIN
   mhalo_fin=total(mhalo_fin*weights)/total(weights)
   mhalo_errlo=total(mhalo_errlo*weights)/total(weights)
   mhalo_errhi=total(mhalo_errhi*weights)/total(weights)
ENDIF

return,[mhalo_fin,mhalo_errlo,mhalo_errhi]
END
