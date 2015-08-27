FUNCTION growth_factor,redshift,$
                       h0=h0, omega_m=omega_m, omega_l=omega_l

IF ~keyword_set(h0) THEN h0=0.702
IF ~keyword_set(omega_m) THEN omega_m=0.275
IF ~keyword_set(omega_l) THEN omega_l=0.725

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
Dz=(gz/g0)/(1.+redshift)


return,dz
END
