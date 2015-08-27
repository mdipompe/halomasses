;+
;  NAME:
;    mhalo2bias
;  PURPOSE:
;    Calculate bias for a given DM halo mass
;
;  USE:
;    bias=mhalo2bias(log_mhalo,z,h0=h0,omega_m=...,model='tinker10')
;
;  INPUT:
;    log_mhalo - Log of the halo mass(es) (in M_sun/h).  Cannot be an
;                array if redshift is
;    redshift - redshift(s).  Cannot be an array if log_mhalo is.
;
;  OPTIONAL INPUT:
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
;    bias (array corresponding to z array or Mhalo array if supplied)    
;
;  NOTES:
;    Hacked from codes by Ryan Hickox
;
;  HISTORY:
;    8-26-15 - Written - MAD (UWyo)
;-
FUNCTION mhalo2bias, log_mhalo, redshift,$
                     h0=h0, omega_m=omega_m, omega_b=omega_b, $
                     omega_l=omega_l, sigma_8=sigma_8, spec_ind=spec_ind, $
                     model=model

IF (n_elements(log_mhalo) GT 1 AND n_elements(redshift) GT 1) THEN BEGIN
   print,'MHALO2BIAS: You can''t supply an array of redshifts AND halo masses. I quit...'
   return,'-1'
ENDIF

;MAD Set default cosmology
IF ~keyword_set(h0) THEN h0=0.702
IF ~keyword_set(omega_m) THEN omega_m=0.275
IF ~keyword_set(omega_b) THEN omega_b=0.046
IF ~keyword_set(omega_l) THEN omega_l=0.725
IF ~keyword_set(sigma_8) THEN sigma_8=0.81
IF ~keyword_set(spec_ind) THEN spec_ind=0.968

;MAD Default model is Tinker et al. (2010)
IF ~keyword_set(model) THEN model='tinker10'

IF (model NE 'tinker10' AND model NE 'tinker05' AND model NE 'smt') THEN BEGIN
   print,'MHALO2BIAS: Unsupported model supplied.'
   print,'            Options are ''tinker10'', ''tinker05'', or ''smt'''
   print,'            Returning ''-1'' and quitting.
   return,'-1'
ENDIF

;MAD Set Hubble parameter scaling E(z)
E_z=sqrt((omega_m*(1.+redshift)^3.) + omega_l)

;MAD Scale omega_m, omega_l with z
omega_m_z=omega_m*((1.+redshift)^3.)/(E_z^2.)
omega_l_z=omega_l/(E_z^2.)

;MAD Get growth factor (using Carroll et al. 1992 approx, 
;MAD accurate to ~1%.  Better than we know halo mass and bias usually)
;MAD Normalized to 1 at z=0
gz=(5./2.)*(omega_m_z/((omega_m_z^(4./7.))-omega_l_z+$
   (1.+(omega_m_z/2.))*(1.+(omega_l_z/70.))))
g0=(5./2.)*(omega_m/((omega_m^(4./7.))-omega_l+$
   (1.+(omega_m/2.))*(1.+(omega_l/70.))))
D_z=(gz/g0)/(1.+redshift)

;MAD Get sigma(M) (van den Bosch 2002, good within 0.005 
;MAD for 10^6 < M < 10^16)
gamma=omega_m*h0*EXP((-1.)*omega_b*(1.+(SQRT(2.*h0)/omega_m)))
c=3.804d-4
x=(c*gamma*((10^(log_mhalo/3.))))/((omega_m)^(1./3.))
g1=64.087*((1.+1.074*(x^(0.3))-1.581*(x^(0.4))+0.954*(x^(0.5))-0.185*(x^(0.6)))^(-10.))
x=(32.*gamma)
g2=64.087*((1.+1.074*(x^(0.3))-1.581*(x^(0.4))+0.954*(x^(0.5))-0.185*(x^(0.6)))^(-10.))
f=(g1^2.)/(g2^2.)
s0=SQRT(f*(sigma_8^2.))
;MAD Tweak for non-unity spectral index
s0=s0*(10^((log_mhalo-14.09)*(1.00-spec_ind)/9.2))/(1.+(1.00-spec_ind)/9.2)
sigma_m=s0*D_z


;MAD Use approximation of NFW 97 for delta (valid in universe with
;MAD Lambda, while Tinker delta_c=1.69 only for Omega_m=1)
delta_c=0.15*((12.*!dpi)^(2./3.))*((omega_m_z)^(0.0055))	

nu=delta_c/sigma_m

IF (model EQ 'tinker10') THEN BEGIN
   print,'MHALO2BIAS: Using Tinker et al. (2010) model...'
   ;MAD Assume delta=200
   y=ALOG10(200.)
   Abig=1.+0.24*y*exp(-(4/y)^4.)
   asmall=0.44*y-0.88
   Bbig=0.183
   bsmall=1.5
   Cbig=0.019+0.107*y+0.19*exp(-(4/y)^4.)
   csmall=2.4
   
   term1 = Abig*((nu^asmall)/((nu^asmall)+(delta_c^asmall)))
   term2 = Bbig*(nu^bsmall)
   term3 = Cbig*(nu^csmall)
   
   bias=1.-term1+term2+term3
ENDIF ELSE BEGIN
   IF (model EQ 'tinker05') THEN BEGIN
      print,'MHALO2BIAS: Using Tinker et al. (2005) model...'
      a=0.707
      b=0.35
      c=0.80
   ENDIF
   IF (model EQ 'smt') THEN BEGIN
      print,'MHALO2BIAS: Using Sheth, Mo, & Turman (2001) model...'
      a=0.707
      b=0.5
      c=0.6
   ENDIF

   denom_1=SQRT(a)*delta_c
   numer_1=a*SQRT(a)*(nu^2.)
   numer_2=SQRT(a)*b*((a*(nu^2.))^(1.-c))
   numer_3=-(a*(nu^2.))^c
   denom_3=((a*(nu^2.))^c)+b*(1.-c)*(1.-c/2.)
   numer_4=numer_3/denom_3	

   bias=1.+(numer_1+numer_2+numer_4)/denom_1
ENDELSE

return,bias
END
