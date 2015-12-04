;+
;  NAME:
;    sigma8
;
;  PURPOSE:
;    Given a matter power spectrum, measure sigma8
;
;  USE:
;    res=sigma8(power_spec,rho_crit=rho_crit,omega_m=omega_m)
;
;  INPUT:
;    power_spec - can be 2xn array of [k,P(k)], or string name of text
;                 file with these two columns.
;
;  OPTIONAL INPUT:
;    rho_crit - critical density (with h=1) Defaults to 2.7745d11
;    omega_m - Omega_matter, defaults to 0.273
;
;  OUTPUT:
;    sigma8 - RMS mass fluctuations at radius 8 Mpc/h
;
;  HISTORY:
;    11-7-15 - Written - MAD (Dartmouth)
;-
FUNCTION sigma8,power_spec,rho_crit=rho_crit,omega_m=omega_m

IF ~keyword_set(rho_crit) THEN rho_crit=2.7745d11
IF ~keyword_set(omega_m) THEN omega_m=0.273
rho_mean=rho_crit*omega_m

check=size(power_spec,/type)
IF (check EQ 4 OR check EQ 5) THEN BEGIN
   wavenum=power_spec[0,*]
   pk=power_spec[1,*]
ENDIF ELSE IF (check EQ 7) THEN BEGIN
   readcol,power_spec,wavenum,pk,format='D'
ENDIF ELSE BEGIN
   print,'SIGMA8: I don''t understand your power spectrum input'
   return,'-1'
ENDELSE

;MAD Interpolate power spectrum to help integration
exponents=cgScaleVector(Findgen(10000),-5,5)
newk=10.^exponents
newpk=interpol(pk,wavenum,newk)
newk=newk[where(newpk GT 0)]
newpk=newpk[where(newpk GT 0)]

;MAD Integrate to get sigma_M. (See eq 3 of van den Bosch 2002)
Rf=8.
WHat=(3./((newk*Rf)^3.))*(SIN(newk*Rf)-((newk*Rf)*COS(newk*Rf)))
integ=newpk*(WHat^2.)*(newk^2.)
res=int_tabulated(newk,integ,/double)
sig8=((1./(2.*(!dpi^2.)))*res)^(1./2.)

return,sig8
END
