;+
;  NAME:
;    growth_factor
;
;  PURPOSE:
;    Calculate the growth factor of linear perterbations at given z,
;    using Carroll et al. (1992) approximation.
;
;  USE:
;    res=growth_factor(z,h0=h0,omega_m=omega_m,omega_l=omega_l)
;
;  INPUT:
;    redshift - redshift(s).
;
;  OPTIONAL INPUT:
;
;  OUTPUT:
;    Growth factor at z.
;
;  NOTES:
;    Hacked from codes by Ryan Hickox
;
;  HISTORY:
;    8-26-15 - Written - MAD (UWyo)
;    12-5-17 - Modified to use common block cosmological_parameters
;              set in load_cosmology.pro, evolve density functions
;-
FUNCTION growth_factor,redshift

  COMMON cosmological_parameters

  ;MAD Evolve densities to z
  omega_m_z=evolve_density(omega_m,redshift,type='matter')
  omega_l_z=evolve_density(omega_l,redshift,type='lambda')

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
