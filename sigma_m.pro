FUNCTION sigma_m, log_mhalo, redshift, $
                  h0=h0, omega_m=omega_m, omega_b=omega_b, $
                  omega_l=omega_l,sigma_8=sigma_8, spec_ind=spec_ind

;MAD Set default cosmology
IF ~keyword_set(h0) THEN h0=0.702
IF ~keyword_set(omega_m) THEN omega_m=0.275
IF ~keyword_set(omega_l) THEN omega_l=0.725
IF ~keyword_set(omega_b) THEN omega_b=0.046
IF ~keyword_set(sigma_8) THEN sigma_8=0.81
IF ~keyword_set(spec_ind) THEN spec_ind=0.968

;MAD Get growth factor
D_z=growth_factor(redshift,h0=h0,omega_m=omega_m,omega_l=omega_l)

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
;MAD multiply by growth factor to get from z=0 to redshift
sigm=s0*D_z

return,sigm
END
