FUNCTION rho_crit,z,physical=physical

  COMMON cosmological_parameters

  ;MAD Put grav constant into (Mpc/h)/((M_sun/h) * s)
  G=(6.67408d-11)*((1.99d30/h))*(1./(3.0857d22/h))^3.
  
  ;MAD Convert little h to big H0, into 1/s
  h0=h*100*(1./(3.0857d19))
  Hz=sqrt((h0^2.)*((omega_m*(1.+z)^3.) + omega_l))

  ;MAD Get rho_crit at z
  rho_crit=(3.*Hz^2.)/(8.*!dpi*G)

  IF (n_elements(physical) NE 0) THEN rho_crit=rho_crit/((1.+z)^3.)
  
  return,rho_crit
END
