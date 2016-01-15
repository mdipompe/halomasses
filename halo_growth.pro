;+
;  NAME:
;    halo_growth
;
;  PURPOSE:
;    Takes a halo mass at a redshift and extrapolates it across cosmic
;    time, using formulae (eqn 2) of  Fakhouri, Ma & Boylan-Kolchin (2010)
;
;  USE:
;    out=halo_growth(m_halo,redshift,z_extrap=z_extrap,h=h,omega_m=omega_m,
;                    lambda=lambda,outfile=outfile,model=model)
;
;  INPUTS:
;    m_halo - the starting halo mass
;    redshift - the starting redshift
;    
;  KEYWORDS:
;
;  OPTIONAL INPUTS:
;    z_out - list of redshifts to output masses for (defaults to 0 to
;               5 in increments of 0.01).  Maximum value of 10.
;    h - little h (H0/100)
;    omega_m - matter density
;    lambda - dark energy density
;    outfile - filename to write hardcopy
;    model - evolution parameterization ('mean' or 'median', defaults to median)
;
;  OUTPUTS:
;    out - structure containing tags z and mhalo
;
;  NOTES:
;    
;  HISTORY:
;    1-12-13 - Written - MAD (Dartmouth) 
;-
FUNCTION halo_growth,m_halo,redshift,z_out=z_out,$
                     h=h,omega_m=omega_m,lambda=lambda,outfile=outfile,model=model

  ;MAD Set defaults
  IF ~keyword_set(h) THEN h=0.702
  IF ~keyword_set(omega_m) THEN omega_m=0.273
  IF ~keyword_set(lambda) THEN  lambda=0.727
  IF ~keyword_set(model) THEN model='median'

  ;MAD Make z array
  IF keyword_set(z_out) THEN BEGIN
     IF (max(z_out) GT 10.) THEN message,'Maximum requested z is too high (>10)'
     IF (max(z_out) GT 5.) THEN zvals=findgen(1001)/100. ELSE $
        zvals=findgen(501)/100.
  ENDIF ELSE BEGIN
     zvals=findgen(501)/100.
  ENDELSE
  
  ;MAD Set constants for mean/median growth rates
  IF (model EQ 'mean') THEN BEGIN
     A=46.1
     B=1.11
  ENDIF ELSE IF (model EQ 'median') THEN BEGIN
     A=25.3
     B=1.65
  ENDIF ELSE BEGIN
     message,'Unknown evolution model supplied (use ''mean'' or ''median'')'
     return,-1
  ENDELSE
  
  ;MAD Sort output redshifts, if needed
  IF keyword_set(z_out) THEN z_out=z_out[bsort(z_out)]

  ;MAD Find nearest z to input redshift, set to initial redshift
  indx=closest(zvals,redshift)
  redshift=zvals[indx]

  ;MAD Determine which, if any, extrapolated redshifts are
  ;MAD below/above input
  lowz=where(zvals LE redshift, n_low)
  highz=where(zvals GE redshift, n_high)

  ;MAD Extrapolate below input redshift (to higher mass halos)
  IF (n_low GT 1) THEN BEGIN
     usez=reverse(zvals[lowz])
     tarray=fltarr(n_elements(usez))
     tmp=cosmocalc(usez[0],h=h,om=omega_m,lambda=lambda)
     tarray[0]=tmp.t_l*1d9
     marray=fltarr(n_elements(usez))
     marray[0]=m_halo
     FOR i=1L,n_elements(usez)-1 DO BEGIN
        IF (usez[i] EQ 0) THEN BEGIN
           tarray[0]=0.
        ENDIF ELSE BEGIN
           tmp=cosmocalc(usez[i],h=h,om=omega_m,lambda=lambda)
           tarray[i]=tmp.t_l*1d9
        ENDELSE
        deltat=tarray[i-1]-tarray[i]
        deltam=(A*(((10.^marray[i-1])/h/1d12)^1.1)*$
                (1+B*usez[i])*sqrt((omega_m*((1+usez[i])^3.))+lambda))*deltat
        marray[i]=alog10((10.^marray[i-1])+(deltam*h))   
     ENDFOR     
  ENDIF
  ;MAD Extrapolate above input redshift (to lower mass halos)
  IF (n_high GT 1) THEN BEGIN
     usez=zvals[highz]
     tarray=fltarr(n_elements(usez))
     tmp=cosmocalc(usez[0],h=h,om=omega_m,lambda=lambda)
     tarray[0]=tmp.t_l*1d9
     marray2=fltarr(n_elements(usez))
     marray2[0]=m_halo
     FOR i=1L,n_elements(usez)-1 DO BEGIN
        tmp=cosmocalc(usez[i],h=h,om=omega_m,lambda=lambda)
        tarray[i]=tmp.t_l*1d9
        deltat=tarray[i]-tarray[i-1]
        deltam=(A*(((10.^marray2[i-1])/h/1d12)^1.1)*$
                (1+B*usez[i])*sqrt((omega_m*((1+usez[i])^3.))+lambda))*deltat
        marray2[i]=alog10((10.^marray2[i-1])-(deltam*h))   
     ENDFOR     
  ENDIF

  ;MAD Remove duplicates of initial mass if needed
  IF (n_low LE 1) THEN mvals=marray2
  IF (n_high LE 1) THEN mvals=reverse(marray)
  IF (n_low GT 1 AND n_high GT 1) THEN $
     mvals=[reverse(marray[0:n_elements(marray)-1]),marray2[1:n_elements(marray2)-1]]

  ;MAD Extrapolate to desired output values, if needed
  IF keyword_set(z_out) THEN BEGIN
     linterp,zvals,mvals,z_out,m_out
  ENDIF ELSE BEGIN
     m_out=mvals
     z_out=zvals
  ENDELSE
  
  ;MAD Make output structure
  out={z:0.,mhalo:0.}
  out=replicate(out,n_elements(m_out))
  out.z=z_out
  out.mhalo=m_out

  ;MAD Write out file, if needed
  IF keyword_set(outfile) THEN BEGIN
     openw,1,outfile
     FOR i=0L,n_elements(out)-1 DO BEGIN
        printf,1,strtrim(out[i].z,2)+'   '+strtrim(out[i].mhalo,2),format='(A)'
     ENDFOR
     close,1
  ENDIF
  
  return,out
END
