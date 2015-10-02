;+
;  NAME:
;    lifetimes
;  PURPOSE:
;    Get quasar lifetime by abundance matching
;
;  USE:
;    lifetimes,log_mhalo,density,zrange,halodensity,lifetime,model='tinker10',$
;              h0=h0,omega_m=omega_m,omega_l=omega_l,omega_b=omega_b,$
;              spec_ind=spec_ind
;
;  INPUT:
;    log_mhalo - log halo mass of population
;    density - space density of population
;    zrange - two element array of redshift range under consideration
;
;  OPTIONAL INPUT:
;    err_mhalo - error on halo mass (can be two element array, lower
;                and upper)
;    h0 - little h (H_0/100). Defaults to 0.702
;    omega_m - omega_matter.  Defaults to 0.275
;    omega_l - omega_lambda. Defaults to 0.725
;    omega_b - omega_baryon. Defaults to 0.046
;    spec_ind - power spectrum spectral index.  Defaults to 0.968
;    model - mass function model.  Options are:
;            'tinker10' (Tinker et al. 2010)  Default
;            'tinker08' (Tinker et al. 2008)
;            'smt' (Sheth, Mo, & Turman 2001)
;
;  OUTPUT:
;    halodensity - the space density of halos with mhalo (Mpc^-3) (3
;                  values - density, lower error, upper error)
;    lifetime - Lifetime of objects (Myr) (3 values - lifetime, lower
;               error, upper error)
;
;  NOTES:
;    Hacked from codes by Ryan Hickox
;
;  HISTORY:
;    10-2-15 - Written - MAD (UWyo)
;-
PRO lifetimes, log_mhalo, density, zrange, halodensity, lifetime, $
               err_mhalo=err_mhalo, $
               h0=h0, omega_m=omega_m, omega_b=omega_b, $
               omega_l=omega_l, $
               spec_ind=spec_ind, $
               model=model,power_spec=power_spec

;MAD Set default cosmology
IF ~keyword_set(h0) THEN h0=0.702
IF ~keyword_set(omega_m) THEN omega_m=0.275
IF ~keyword_set(omega_b) THEN omega_b=0.046
IF ~keyword_set(omega_l) THEN omega_l=0.725
IF ~keyword_set(spec_ind) THEN spec_ind=0.968
IF ~keyword_set(model) THEN model='tinker10'

;MAD Set error to 0 if needed
IF ~keyword_set(err_mhalo) THEN err_mhalo=[0.,0.]
IF (n_elements(err_mhalo) EQ 1) THEN err_mhalo=[err_mhalo,err_mhalo]

;MAD Get mean redshift
meanz=mean(zrange)

;MAD Set up range of halo masses, get dn/dlog(M)
testmhalo=findgen(200)/20+8
testdense=(halo_mass_function(testmhalo,meanz,$
                              h0=h0,omega_m=omega_m,omega_b=omega_b,$
                              omega_l=omega_l,spec_ind=spec_ind,$
                              model=model,power_spec=power_spec))*(10.^testmhalo)*ALOG(10)

;MAD Interpolate to desired value/errors
linterp,testmhalo,testdense,[log_mhalo-err_mhalo[0],log_mhalo,$
                             log_mhalo+err_mhalo[1]],halodense

;MAD Multiply out factors of h
halodense=halodense*(h0^3.)

;MAD Define duty cycle, cosmic time range
duty_cycle=density/halodense
tzlo=cosmocalc(zrange[0])
tzlo=tzlo.t_l
tzhi=cosmocalc(zrange[1])
tzhi=tzhi.t_l

;MAD Get lifetime in Myr
lifetime = duty_cycle*(tzhi-tzlo)*1000.

;MAD Set output
halodensity=[halodense[1],halodense[0]-halodense[1],halodense[1]-halodense[2]]
lifetime=[lifetime[1],lifetime[2]-lifetime[1],lifetime[1]-lifetime[0]]

return
END
