FUNCTION halo_mass_function, log_mhalo, redshift,$
                             h0=h0, omega_m=omega_m, omega_b=omega_b, $
                             omega_l=omega_l, sigma_8=sigma_8, $
                             spec_ind=spec_ind, $
                             model=model

IF (n_elements(redshift) GT 1) THEN BEGIN
   print,'HALO_MASS FUNCTION: One redshift at a time please!'
   return,'-1'
ENDIF

;MAD Set default cosmology
IF ~keyword_set(h0) THEN h0=0.702
IF ~keyword_set(omega_m) THEN omega_m=0.275
IF ~keyword_set(omega_b) THEN omega_b=0.046
IF ~keyword_set(omega_l) THEN omega_l=0.725
IF ~keyword_set(sigma_8) THEN sigma_8=0.81
IF ~keyword_set(spec_ind) THEN spec_ind=0.968

IF ~keyword_set(model) THEN model='tinker10'

;MAD Set Hubble parameter scaling E(z)
E_z=sqrt((omega_m*(1.+redshift)^3.) + omega_l)

;MAD Scale omega_m, omega_l with z
omega_m_z=omega_m*((1.+redshift)^3.)/(E_z^2.)

;MAD Get sigma_m
sigma_m=sigma_m(log_mhalo,redshift,h0=h0,omega_m=omega_m,omega_b=omega_b,$
                omega_l=omega_l,sigma_8=sigma_8,spec_ind=spec_ind)

;MAD Use approximation of NFW 97 for delta (valid in universe with
;MAD Lambda, while Tinker delta_c=1.69 only for Omega_m=1)
delta_c=0.15*((12.*!dpi)^(2./3.))*((omega_m_z)^(0.0055))	

nu=delta_c/sigma_m

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
   term2=(sigma_M/b[0])^(-asmall[0])+1.
   term3=exp(-c[0]/sigma_M^2.)
   
   f=term1*term2*term3/nu
ENDIF

;MAD Get N(M) in physical units
rho_crit=2.7745d11
rho_matter=rho_crit

mhalo=10.^(log_mhalo)
   
d_sigma_M=sigma_M(log_mhalo+0.001,redshift)
d_lognu=ALOG10(delta_c/d_sigma_M)-ALOG10(nu)
dlognu_dlogM=d_lognu/0.001

n_mhalo=(nu*f*rho_matter*dlognu_dlogM)/(mhalo^2.)

return,n_mhalo

END
