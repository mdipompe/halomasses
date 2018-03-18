;+
;  NAME:
;    halo_mass_function
;  PURPOSE:
;    Get N(M) for a halo mass (or masses)
;
;  USE:
;    res=halo_mass_function(log_mhalo, z, [model='model',...])
;
;  INPUT:
;    log_mhalo - Log of the halo mass(es) (in M_sun/h). 
;    redshift - redshift.
;
;  OPTIONAL INPUT:
;    model - mass function model.  Options are:
;            'tinker10' (Tinker et al. 2010)  Default
;            'tinker08' (Tinker et al. 2008)
;            'smt' (Sheth, Mo, & Turman 2001)
;    power_spec - filename of matter power spectrum (gets with CAMB
;                 if not supplied)
;    Delta - Overdensity (relative to mean background) cut used to
;            define halo mass. Setting only works for Tinker et
;            al. (2008, 2010) models.  Defaults to 200.
;
;  OPTIONAL KEYWORDS:
;    silent - don't print anything, except errors
;    logM - set to return dN/dlogM (instead of dn/dM)
;
;  OUTPUT:
;    Normalized N(M)
;
;  NOTES:
;    Hacked from codes by Ryan Hickox
;
;  HISTORY:
;    8-26-15 - Written - MAD (UWyo)
;   10-29-15 - Added Delta keyword - MAD (Dartmouth)
;    12-1-16 - Edited to use common block from load_cosmology
;     3-4-18 - Added evolution of params with z in Tinker10 model
;-
FUNCTION halo_mass_function, log_mhalo, redshift,$
                             model=model,power_spec=power_spec,$
                             Delta=Delta,silent=silent,logM=logM

  COMMON cosmological_parameters
  
  IF (n_elements(redshift) GT 1) THEN message,'One redshift at a time please!'

  ;MAD Set default params
  IF ~keyword_set(model) THEN model='tinker10'
  IF ~keyword_set(Delta) THEN Delta=200.

  ;MAD Scale omega_m, omega_l with z
  omega_m_z=evolve_density(redshift,type='matter')

  ;MAD Get sigma_m
  sig_m=sigma_m(log_mhalo,redshift,power_spec=power_spec,silent=silent)

  ;MAD Use approximation of NFW 97 for delta (valid in universe with
  ;MAD Lambda, while Tinker delta_c=1.69 only for Omega_m=1)
  delta_c=0.15*((12.*!dpi)^(2./3.))*((omega_m_z)^(0.0055))	

  nu=delta_c/sig_m
  
  IF (model EQ 'tinker10') THEN BEGIN
     ;MAD Find fit params for Delta (interpolate Tinker et al. (2010) Table 4)
     Deltas=[200.,300.,400.,600.,800.,1200.,1600.,2400.,3200.]
     alphas=[0.368,0.363,0.385,0.389,0.393,0.365,0.379,0.355,0.327]
     betas=[0.589,0.585,0.544,0.543,0.564,0.623,0.637,0.673,0.702]
     gams=[0.864,0.922,0.987,1.09,1.20,1.34,1.50,1.68,1.81]
     phis=[-0.729,-0.789,-0.910,-1.05,-1.20,-1.26,-1.45,-1.50,-1.49]
     etas=[-0.243,-0.261,-0.261,-0.273,-0.278,-0.301,-0.301,-0.319,-0.336]

     beta=spline(deltas,betas,delta)
     gam=spline(deltas,gams,delta)
     phi=spline(deltas,phis,delta)
     eta=spline(deltas,etas,delta)
          
     ;MAD Evolve params with redshift (Tinker 2010 Eq 9-12)
     ;MAD but if z > 3, use z=3 params
     IF (redshift GT 3) THEN z_use=3. ELSE z_use=redshift
     beta=beta*(1.+z_use)^(0.20)
     phi=phi*(1.+z_use)^(-0.08)
     eta=eta*(1.+z_use)^(0.27)
     gam=gam*(1.+z_use)^(-0.01)

     ;MAD some conditions must be met for normalization to work
     IF (gam LE 0.) THEN BEGIN
        print,'Gamma must be > 0, setting to 0.001'
        gam=1d-3
     ENDIF
     IF (eta LE -0.5) THEN BEGIN
        print,'Eta must be > -0.5, setting to -0.499'
        eta=-0.499
     ENDIF
     IF (eta-phi LE -0.5) THEN BEGIN
        print,'eta-phi must be > -0.5, setting phi to eta+0.499'
        phi=eta+0.499
     ENDIF
     IF (beta LE 0) THEN BEGIN
        print,'beta must be > 0, setting to 0.001'
        beta=1d-3
     ENDIF

     ;MAD Calculate normalization (this is hacked from the HMFcalc codes,
     ;MAD I honestly don't know where it comes from...
     ;https://github.com/steven-murray/hmf/blob/master/hmf/fitting_functions.py
     xx=where(deltas EQ delta,cnt)
     IF (cnt EQ 0. or redshift NE 0.) THEN BEGIN
        term1=2.^(eta-phi-0.5)
        term2=beta^(-2.*phi)
        term3=gam^(-0.5-eta)
        term4_1=(2.^phi)*(beta^(2.*phi))
        term4_2=gamma(eta+0.5)
        term4_3=(gam^phi)*gamma(0.5+eta-phi)
        alpha=1./(term1*term2*term3*(term4_1*term4_2+term4_3))
     ENDIF ELSE BEGIN
        alpha=alphas[xx]
     ENDELSE
     
     term1=1+((beta*nu)^((-2)*phi))
     term2=nu^(2.*eta)
     term3=EXP(((-1.)*gam*(nu^2.))/2.)
     
     f=alpha*term1*term2*term3
  ENDIF
  IF (model EQ 'tinker08') THEN BEGIN
          abigs=[0.1858659,$
            0.1995973,$
            0.2115659,$
            0.2184113,$
            0.2480968,$
            0.2546053,$
            0.2600000,$
            0.2600000,$
            0.2600000]
     
     asmalls=[1.466904,$
              1.521782,$
              1.559186,$
              1.614585,$
              1.869936,$
              2.128056,$
              2.301275,$
              2.529241,$
              2.661983]
     
     bs=[2.571104,$
         2.254217,$
         2.048674,$
         1.869559,$
         1.588649,$
         1.507134,$
         1.464374,$
         1.436827,$
         1.405210]
     
     cs=[1.193958,$
         1.270316,$
         1.335191,$
         1.446266,$
         1.581345,$
         1.795050,$
         1.965613,$
         2.237466,$
         2.439729]
     
     deltas=[200, 300, 400, 600, 800, 1200, 1600, 2400, 3200.]
     
     abig0=spline(deltas,abigs,delta)
     asmall0=spline(deltas,asmalls,delta)
     b0=spline(deltas,bs,delta)
     c=spline(deltas,cs,delta)

     logalpha=-(0.75/alog10(Delta/75))^1.2
     alpha=10.^logalpha
     
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
  rho_crit=rho_crit(0.,/physical)
  rho_matter=rho_crit*omega_m

  mhalo=10.^(log_mhalo)
   
  d_sig_M=sigma_m(log_mhalo+0.001,redshift,power_spec=power_spec,silent=silent)
  d_lognu=ALOG10(delta_c/d_sig_M)-ALOG10(nu)
  dlognu_dlogM=d_lognu/0.001
  
  n_mhalo=(nu*f*rho_matter*dlognu_dlogM)/(mhalo^2.)
  
  IF (n_elements(logM) NE 0) THEN n_mhalo=n_mhalo*(mhalo*alog(10))
  
  return,n_mhalo

END
