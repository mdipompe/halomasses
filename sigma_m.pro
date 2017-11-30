;+
;  NAME:
;    sigma_m
;  PURPOSE:
;    Calculate RMS of mass density field for given halo
;    mass/redshift.  Will call CAMB4IDL to get matter power spectrum
;    if necessary.
;
;  USE:
;    result=sigma_m(log_mhalo,z,h0=h0...,CAMB parameters...)
;
;  INPUT:
;    log_mhalo - Log of the halo mass(es) (in M_sun/h).
;    redshift - redshift.  Cannot be an array currently.
;               Note that if you have many redshifts and want to do it
;               quickly, you can just do z=0 and multiply by growth
;               factor, which is slightly less accurate.
;
;  OPTIONAL INPUT:
;    h0 - little h (H_0/100). Defaults to 0.702
;    omega_m - omega_matter.  Defaults to 0.275
;    omega_l - omega_lambda. Defaults to 0.725
;    omega_b - omega_baryon. Defaults to 0.046
;    spec_ind - power spectrum spectral index.  Defaults to 0.96
;    power_spec - string name of file with power spectrum.  Two
;                 columns, k (h/Mpc) and P(k) (h^-3 Mpc^3).
;                 If not supplied, runs camb with CAMB4IDL to get it.
;    maxk - maximum k value to run CAMB to.  Will interpolate this out
;           to 10^5.  Defaults to 50.
;    paramfile - CAMB paramfile of defaults.  Cosmology and others
;                will be set in call below, but you might want to play
;                with others.  Defaults to default_params.ini.
;
;  OPTIONAL KEYWORDS:
;    silent - don't print anything (except errors)
;
;  OUTPUT:
;    sigma_M
;    If camb is called, several camb files will be written in current directory
;
;  NOTES:
;    
;  HISTORY:
;    8-27-15 - Written - MAD (UWyo)
;-
FUNCTION sigma_m, log_mhalo, redshift, $
                  h0=h0, omega_m=omega_m, omega_b=omega_b, $
                  omega_l=omega_l, spec_ind=spec_ind,$
                  power_spec=power_spec,maxk=maxk,$
                  paramfile=paramfile,silent=silent
                  
IF (n_elements(redshift) GT 1) THEN BEGIN
   print,'SIGMA_M: Currently only supports single redshifts...I quit'
   print,'         (Note that you can just do z=0 and multiply by'
   print,'         growth factor (growth_factor.pro) as a good'
   print,'         approximation.)'
   return,'-1'
ENDIF

;MAD Set default cosmology
IF ~keyword_set(h0) THEN h0=0.702
IF ~keyword_set(omega_m) THEN omega_m=0.275
IF ~keyword_set(omega_l) THEN omega_l=0.725
IF ~keyword_set(omega_b) THEN omega_b=0.046
IF ~keyword_set(spec_ind) THEN spec_ind=0.96

IF ~keyword_set(power_spec) THEN BEGIN
   ;MAD setup redshifts and files for CAMB
   tfile='transfer_' + strtrim(redshift,2) + '.dat'
   mpfile='matterpower_' + strtrim(redshift,2) + '.dat'

   ;MAD Set some camb defaults
   IF ~keyword_set(maxk) THEN maxk=50.
   IF ~keyword_set(paramfile) THEN paramfile='default_params.ini'

   res=camb4idl(/runcamb, paramfile=paramfile,output_root='camb', $
                get_scalar='T',get_transfer='T', $
                get_tensor='F',do_lensing='T', $
                do_nonlinear=0, $
                scalar_spectral_index=spec_ind, $
                transfer_kmax=maxk, $
                transfer_num_redshifts=n_elements(redshift), $
                transfer_redshift=redshift, $
                transfer_filename=tfile, transfer_matterpower=mpfile, $
                use_physical='F', $
                hubble=h0*100., $
                omega_baryon=omega_b, $
                omega_cdm=omega_m-omega_b, $
                omega_lambda=omega_l)
   readcol,'camb_'+mpfile,wavenum,pk,format='D,D',silent=silent
ENDIF ELSE BEGIN
   readcol,power_spec,wavenum,pk,format='D,D',silent=silent
ENDELSE

mhalo=10.^log_mhalo
IF (n_elements(z) GT 1) THEN rho_crit=rho_crit(0.,/phys) ELSE $
   rho_crit=rho_crit(z,/phys)
rho_mean=rho_crit*omega_m

;MAD Interpolate power spectrum to improve integration
exponents=cgScaleVector(Findgen(10000),-5,5)
newk=10.^exponents
newpk=interpol(pk,wavenum,newk)
newk=newk[where(newpk GT 0)]
newpk=newpk[where(newpk GT 0)]

;MAD Integrate to get sigma_M. (See eq 3 of van den Bosch 2002)
Rf=((3./4.)*(mhalo/(!dpi*rho_mean)))^(1./3)
sigm=fltarr(n_elements(mhalo))
FOR i=0L,n_elements(mhalo)-1 DO BEGIN
   WHat=(3./((newk*Rf[i])^3.))*(SIN(newk*Rf[i])-((newk*Rf[i])*COS(newk*Rf[i])))
   integ=newpk*(WHat^2.)*(newk^2.)
   res=int_tabulated(newk,integ,/double)
   sigm[i]=((1./(2.*(!dpi^2.)))*res)^(1./2.)
ENDFOR


;MAD Check integral accuracy by normalizing to known sigma_8, then
;MAD integrating and seeing if you get that value back.
;IF keyword_set(check_acc) THEN BEGIN
;   sigma_8=0.8
;   Rf=8.
;   WHat=(3./((newk*Rf)^3.))*(SIN(newk*Rf)-((newk*Rf)*COS(newk*Rf)))
;   integ=((newpk^2.)*newk)*(WHat^2.)*(newk^2.)
;   res1=((1./(2.*(!dpi^2.)))*int_tabulated(newk,integ,/double))^(1./2.)
;   correc=((sigma_8)/res1)^2.
;   integ=(correc*(newpk^2.)*newk)*(WHat^2.)*(newk^2.)
;   res2=((1./(2.*(!dpi^2.)))*int_tabulated(newk,integ,/double))^(1./2.)
;   acc=100.*((res2-(sigma_8))/(sigma_8))
;   print,'SIGMA_8: Power spectrum integration accurate to'
;   print,'        ' + strtrim(abs(acc),2) + '%'
;ENDIF

return,sigm
END
