;+
;  NAME:
;    halo_mass_function
;  PURPOSE:
;    Get N(M) for a halo mass (or masses)
;
;  USE:
;    res=halo_mass_function
;
;  INPUT:
;    log_mhalo - Log of the halo mass(es) (in M_sun/h). 
;    redshift - redshift.
;
;  OPTIONAL INPUT:
;    h0 - little h (H_0/100). Defaults to 0.702
;    omega_m - omega_matter.  Defaults to 0.275
;    omega_l - omega_lambda. Defaults to 0.725
;    omega_b - omega_baryon. Defaults to 0.046
;    spec_ind - power spectrum spectral index.  Defaults to 0.968
;    model - mass function model.  Options are:
;            'tinker10' (Tinker et al. 2010)  Default
;            'tinker08' (Tinker et al. 2008)
;            'smt' (Sheth, Mo, & Turman 2001)
;    power_spec - filename of matter power spectrum (gets with CAMB
;                 if not supplied)
;
;  OUTPUT:
;    Normalized N(M)
;
;  NOTES:
;    Hacked from codes by Ryan Hickox
;
;  HISTORY:
;    8-26-15 - Written - MAD (UWyo)
;-
FUNCTION halo_mass_function, log_mhalo, redshift,$
                             h0=h0, omega_m=omega_m, omega_b=omega_b, $
                             omega_l=omega_l, $
                             spec_ind=spec_ind, $
                             model=model,power_spec=power_spec

IF (n_elements(redshift) GT 1) THEN BEGIN
   print,'HALO_MASS FUNCTION: One redshift at a time please!'
   return,'-1'
ENDIF

;MAD Set default cosmology
IF ~keyword_set(h0) THEN h0=0.702
IF ~keyword_set(omega_m) THEN omega_m=0.275
IF ~keyword_set(omega_b) THEN omega_b=0.046
IF ~keyword_set(omega_l) THEN omega_l=0.725
IF ~keyword_set(spec_ind) THEN spec_ind=0.968
IF ~keyword_set(model) THEN model='tinker10'

;MAD Set Hubble parameter scaling E(z)
E_z=sqrt((omega_m*(1.+redshift)^3.) + omega_l)

;MAD Scale omega_m, omega_l with z
omega_m_z=omega_m*((1.+redshift)^3.)/(E_z^2.)

;MAD Get sigma_m
sig_m=sigma_m(log_mhalo,redshift,h0=h0,omega_m=omega_m,omega_b=omega_b,$
                 omega_l=omega_l,spec_ind=spec_ind,power_spec=power_spec)

;MAD Use approximation of NFW 97 for delta (valid in universe with
;MAD Lambda, while Tinker delta_c=1.69 only for Omega_m=1)
delta_c=0.15*((12.*!dpi)^(2./3.))*((omega_m_z)^(0.0055))	

nu=delta_c/sig_m

IF (model EQ 'tinker10') THEN BEGIN
   ;MAD For delta=200
   alpha=0.368
   beta=0.589
   gam=0.864
   phi=-0.729
   eta=-0.243

   term1=1+((beta*nu)^((-2)*phi))
   term2=nu^(2.*eta)
   term3=EXP(((-1.)*gam*(nu^2.))/2.)

   f=alpha*term1*term2*term3
ENDIF
IF (model EQ 'tinker08') THEN BEGIN
   delta=200
   logalpha=-(0.75/alog10(delta/75))^1.2
   alpha=10.^logalpha
   
   abig0=0.186
   asmall0=1.47
   b0=2.57
   c=1.19
   
   abig=abig0*(1.+redshift)^(-0.14)
   asmall=asmall0*(1.+redshift)^(-0.06)
   b=b0*(1.+redshift)^(-alpha)
   
   term1=abig
   term2=(sig_m/b[0])^(-asmall[0])+1.
   term3=exp(-c[0]/sig_m^2.)
   
   f=term1*term2*term3/nu
ENDIF
IF (model EQ 'smt') THEN BEGIN
   A=0.3222
   q=0.3
   smalla=0.707
   
   nuprime=sqrt(smalla)*nu
   
   term1=(2.*A)/nu
   term2=(1.+(1./(nuprime^(2.*q))))
   term3=((nuprime^2.)/(2.*!dpi))^(1./2.)
   term4=EXP((-1.)*(nuprime^2.)/2.)

   f=term1*term2*term3*term4
ENDIF


;MAD Get N(M) in physical units
rho_crit=2.7745d11
rho_matter=rho_crit*omega_m

mhalo=10.^(log_mhalo)
   
d_sig_M=sigma_m(log_mhalo+0.001,redshift,h0=h0,omega_m=omega_m,omega_b=omega_b,$
                 omega_l=omega_l,spec_ind=spec_ind,power_spec=power_spec)
d_lognu=ALOG10(delta_c/d_sig_M)-ALOG10(nu)
dlognu_dlogM=d_lognu/0.001

n_mhalo=(nu*f*rho_matter*dlognu_dlogM)/(mhalo^2.)

return,n_mhalo

END
