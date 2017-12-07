FUNCTION evolve_density,z,type=type

  ;MAD Take current density parameters and adjust them to given z
  
  COMMON cosmological_parameters

  IF (n_elements(type) EQ 0) THEN $
     message,'Need to specify the density type (matter or lambda)'

  E_z=sqrt((omega_m*(1.+z)^3.) + omega_l)

  IF (type EQ 'matter') THEN $
     omega_z=omega_m*((1.+z)^3.)/(E_z^2.)
  IF (type EQ 'lambda') THEN $
     omega_z=omega_l/(E_z^2.)

  return,omega_z
END
