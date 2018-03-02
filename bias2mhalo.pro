;+
;  NAME:
;    bias2mhalo
;  PURPOSE:
;    Convert bias to halo mass
;
;  USE:
;    mhalo=bias2mhalo(bias,z,bias_err=bias_err,weights=weights,model='tinker10')
;
;  INPUT:
;    bias - single bias value
;    zs - z(s).  Returns an array of mhalo if weights not supplied,
;         otherwise returns a weighted average
;
;  OPTIONAL INPUT:
;    bias_err - error on the bias.  Assumed symmetric
;    weights - redshift weighing.  Should be array of same size as
;              redshifts. If not supplied, returns bias at each z.
;    model - halo collapse model.  Options are:
;            'tinker10' (Tinker et al. 2010)  Default
;            'tinker05' (Tinker et al. 2005)
;            'smt' (Sheth, Mo, & Turman 2001)
;    power_spec - string name of file containing matter power
;                 spectrum.
;    Delta - Overdensity (relative to mean background) cut used to
;            define halo mass.  Setting only works for Tinker et
;            al. (2010) model. Defaults to 200.
;
;  OPTIONAL KEYWORDS
;    use_gf - use growth factor to scale sigma_m from z=0.  This is
;             done automatically if you supply more than one z to save
;             time.
;    silent - don't print anything
;
;  OUTPUT:
;    Array containing mhalo, lower mhalo error, upper mhalo error.
;
;  NOTES:
;    Partially hacked from codes by Ryan Hickox
;
;  HISTORY:
;    8-26-15 - Written - MAD (UWyo)
;    9-18-15 - Added power_spec input option - MAD (Dartmouth)
;   10-29-15 - Added Delta keyword - MAD (Dartmouth)
;    3-01-18 - Added use of common block for cosmology parameters
;              (from load_cosmology.pro)
;-
FUNCTION bias2mhalo, bias, zs, bias_err=bias_err, $
                     model=model,weights=weights,power_spec=power_spec, $
                     use_gf=use_gf, Delta=Delta,silent=silent

  COMMON cosmological_parameters

  ;MAD Check weights
  IF keyword_set(weights) THEN BEGIN
     IF n_elements(weights) NE n_elements(zs) THEN Begin
        print,'BIAS2MHALO: Not an equal number of weights as redshifts...quitting'
        return,'-1'
     ENDIF
  ENDIF
   
  ;MAD Default model is Tinker et al. (2010) with Delta=200.
  IF ~keyword_set(model) THEN model='tinker10'
  IF ~keyword_set(Delta) THEN Delta=200.

  ;MAD Set errors to 0 if needed
  IF ~keyword_set(bias_err) THEN bias_err=0.

  ;MAD Generate array of log(mhalo) to test/interpolate
  ;MAD (from 6 to 16, where sigma(M) best estimated)
  masses=6.+findgen(1000)*(0.01)

  ;MAD Loop over z
  FOR i=0L,n_elements(zs)-1 DO BEGIN
     ;MAD Get bias for each test Mhalo
     IF (n_elements(zs) EQ 1) THEN BEGIN
        biases=mhalo2bias(masses, zs[i],$
                          model=model,power_spec=power_spec,$
                          use_gf=use_gf,Delta=Delta,silent=silent)
     ENDIF ELSE BEGIN
        biases=mhalo2bias(masses, zs[i],$
                          model=model,power_spec=power_spec,$
                          /use_gf,silent=silent)
     ENDELSE

     ;MAD Interpolate to desired bias
     quadterp, biases, masses, [bias-bias_err,bias,bias+bias_err],mhalo
     IF (n_elements(mhalo_errlo) EQ 0) THEN mhalo_errlo=mhalo[1]-mhalo[0] ELSE $
        mhalo_errlo=[mhalo_errlo,mhalo[1]-mhalo[0]]
     IF (n_elements(mhalo_errhi) EQ 0) THEN mhalo_errhi=mhalo[2]-mhalo[1] ELSE $
        mhalo_errhi=[mhalo_errhi,mhalo[2]-mhalo[1]]
     IF (n_elements(mhalo_fin) EQ 0) THEN mhalo_fin=mhalo[1] ELSE $
        mhalo_fin=[mhalo_fin,mhalo[1]]
  ENDFOR

  ;IF ~keyword_set(weights) THEN weights=fltarr(n_elements(redshift))+1.
  
  IF keyword_set(weights) THEN BEGIN
     mhalo_fin=total(mhalo_fin*weights)/total(weights)
     mhalo_errlo=total(mhalo_errlo*weights)/total(weights)
     mhalo_errhi=total(mhalo_errhi*weights)/total(weights)
  ENDIF
   
  return,transpose([[mhalo_fin],[mhalo_errlo],[mhalo_errhi]])
END
