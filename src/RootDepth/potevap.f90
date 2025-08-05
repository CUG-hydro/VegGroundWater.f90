subroutine POTEVAP_Priestly_Taylor(i, j, tempk, rad, presshp, pet)
  integer :: i, j
  real, parameter :: cp = 1013.*1e-6
  real :: tempk, tempc, rad, presskp, presshp, pet
  real :: alpha, delta, gamma, lambda

  tempc = tempk - 273.15     !C
  presskp = presshp*0.1       !kPa
  rad = rad*24.*3600.*1.e-6 !MJ/day/m2

  alpha = 1.26
  delta = 0.2*(0.00738*tempc + 0.8072)**7.-0.000116
  lambda = 2.501 - 0.002361*tempc
  gamma = (cp*presskp)/(0.622*lambda)

  pet = alpha*rad*delta/(delta + gamma)
  pet = pet/lambda
!if(i.eq.2761.and.j.eq.571)write(6,*)'mirar pt',tempc,presskp,rad,delta,lambda,gamma
end subroutine potevap_priestly_taylor

subroutine POTEVAP_Penman_Monteith(i, j, tempk, rad, rshort, press, qair, wind, lai, veg, hveg, pet)
  real, parameter :: cp = 1013., vk = 0.41, Rd = 287.
  integer :: i, j
  real :: tempk, rad, rshort, press, qair, wind, lai, veg, pet
  real :: zm, hveg, hdisp, z0m, zh, z0h, hhveg
  real :: slai, frad, fswp, fvpd, g_d, ra, rs
  real :: delta, pressesat, pressvap, lambda, gamma, tempc, vpd, dens
  real, dimension(30) :: rl
  real, dimension(0:30) :: zdis
  data rl/150., 150., 500., 500., 175., 240., 110., 100., 250., 150. &
    , 80., 225., 225., 250., 180., 180., 240., 500., 240., 500. &
    , 175., 250., 250., 175., 225., 150., 110., 180., 250., 250./
  data zdis/0.1, 0.1, 0.1, 15., 20., 15., 20., .2, 1., .1, .5, .1, 1., 1., 20., .7, .7, 1. &
    , 10.2, 20.7, 9.2, 7.2, 6.5, 7.4, 3.6, 1.4, .2, .2, .2, .2, 1.1/

  tempc = tempk - 273.15
  pressesat = 610.8*exp(17.27*tempc/(tempc + 237.3)) !Pa
!pressvap = pressesat*rh !Pa
  pressvap = qair*press/(0.622 + qair) !Pa
!delta = 4098. * 610.8 * exp( 17.27*tempc / (tempc+237.3) )  / (tempc+237.3)**2. !Pa/K
!pressvap = pressesat*rh !Pa
  delta = 4098.*pressesat/(tempc + 237.3)**2. !Pa/K
  vpd = pressesat - pressvap
  lambda = (2.501 - 0.002361*tempc)*1.e6 !J/kg
  gamma = (cp*press)/(0.622*lambda)
  dens = press/(Rd*tempk*(1.+0.608*qair))

  if (i == 300 .and. j .eq. 300) write (6, *) 'mirar 1 forcing', tempk, rad, rshort, press, qair, wind, lai, veg, hveg
  if (i == 3 .and. j .eq. 3) write (6, *) 'mirar 2 forcing', tempk, rad, rshort, press, qair, wind, lai, veg, hveg

!hveg=min(10.,hhveg)
!hveg=zdis(nint(veg))
!aerodynamic resistance
  zm = max(10., hveg)
!      hdisp = 0.7*hveg
!      if (hveg.lt.10.)then
!             zm=10.
!      else
!             zm=10+hdisp
!      endif
  z0m = 0.1*hveg
  zh = max(2., hveg)
!      if (hveg.lt.2.)then
!             zh=2.
!      else
!             zh=2.+hveg
!      endif
  z0h = 0.1*z0m
  ra = log((zm - hdisp)/z0m)*log((zh - hdisp)/z0h)/(vk**2.*wind)

!bulk surface resistance
  frad = min(1., (0.004*rshort + 0.05)/(0.81*(1.+0.004*rshort)))
  fswp = 1.

  if (hdisp .gt. 2.) then
    g_d = 0.0003
  else
    g_d = 0.
  end if

  fvpd = exp(-g_d*vpd)
  slai = 0.5*lai

  if (slai*frad*fswp*fvpd .eq. 0.) then
    rs = 5000.
  else
    rs = min(rl(nint(veg))/(slai*frad*fswp*fvpd), 5000.)
  end if

!pet
  pet = (delta*rad + dens*cp*vpd/ra)/(delta + gamma*(1.+rs/ra))
  pet = 3.*3600.*pet/lambda

  if (i == 300 .and. j .eq. 300) write (6, *) 'mirar 1 results', delta, pressesat, pressvap, lambda, gamma, dens, ra, rs, pet
  if (i == 3 .and. j .eq. 3) write (6, *) 'mirar 2 results', delta, pressesat, pressvap, lambda, gamma, dens, ra, rs, pet
end subroutine potevap_penman_monteith

subroutine POTEVAP_Shutteworth_Wallace(i, j, deltat, tempk, rad, rshort, press, qair, wind, lai, veg, hhveg &
                                       , delta, gamma, lambda, ra_a, ra_c, rs_c, R_a, R_s &
                                       , pet_s, pet_c, pet_w, pet_i, floodflag)
  real, parameter :: cp = 1013., vk = 0.41, Rd = 287.
  integer :: i, j, floodflag
  real :: deltat, tempk, rad, rshort, press, qair, wind, lai, veg, pet
  real :: zm, hveg, z0m, zh, z0h, hhveg, z0g
  real :: slai, frad, fswp, fvpd, g_d, ra, rs
  real :: smoi, smoiwp, smoifc
  real :: delta, pressesat, pressvap, lambda, gamma, tempc, vpd, dens
  real :: Rn_s, za, z0c, c_d, d0, ustar, K_h, n, dp, Z0, ra_a, ra_s, uc, wleaf, rb, ra_c, rs_c, rs_s
  real :: R_a, R_c, R_s, C_c, C_s, pet_c, pet_s, pet_w, pet_i
  real, dimension(30) :: rl
  real, dimension(0:30) :: zdis, z0gr, wmax
  real, dimension(2, 0:30) :: bioparms
  data rl/150., 150., 500., 500., 175., 240., 110., 100., 250., 150. &
    , 80., 225., 225., 250., 180., 180., 240., 500., 240., 500. &
    , 175., 250., 250., 175., 225., 150., 110., 180., 250., 250./
  data zdis/0.1, 0.1, 0.1, 15., 20., 15., 20., .2, 1., .1, .5, .1, 1., 1., 20., .7, .7, 1. &
    , 10.2, 20.7, 9.2, 7.2, 6.5, 7.4, 3.6, 1.4, .2, .2, .2, .2, 1.1/
  data bioparms/ &
    .001, 0. & !  0  Ocean
    , .001, 0. & !  1  Lakes, rivers, streams (inland water)
    , .001, 0. & !  2  Ice cap/glacier
    , .02, .001 & !  3  Evergreen needleleaf tree
    , .02, 0.001 & !  4  Deciduous needleleaf tree
    , .02, 0.08 & !  5  Deciduous broadleaf tree
    , .02, 0.05 & !  6  Evergreen broadleaf tree
    , .01, 0.01 & !  7  Short grass
    , .01, 0.01 & !  8  Tall grass
    , .001, 0.01 & !  9  Desert
    , .01, 0.01 & ! 10  Semi-desert
    , .01, 0.01 & ! 11  Tundra
    , .02, 0.01 & ! 12  Evergreen shrub
    , .02, 0.01 & ! 13  Deciduous shrub
    , .02, 0.04 & ! 14  Mixed woodland
    , .005, 0.01 & ! 15  Crop/mixed farming
    , .005, 0.01 & ! 16  Irrigated crop
    , .01, 0.01 & ! 17  Bog or marsh
    !LDAS LSPs, but emissivity based on above
    , .01, .001 & ! 18  Evergreen needleleaf forest
    , .02, .05 & ! 19  Evergreen broadleaf forest
    , .02, .001 & ! 20  Deciduous needleleaf forest
    , .02, .08 & ! 21  Deciduous broadleaf forest
    , .01, .01 & ! 22  Mixed cover
    , .02, .04 & ! 23  Woodland
    , .02, .01 & ! 24  Wooded grassland
    , .02, .01 & ! 25  Closed shrubland
    , .02, .01 & ! 26  Open shrubland
    , .01, .01 & ! 27  Grassland
    , .005, .01 & ! 28  Cropland
    , .001, .01 & ! 29  Bare ground
    , .02, 0./     ! 30  Urban and built up

  z0gr(:) = bioparms(1, :)
  wmax(:) = bioparms(2, :)

  tempc = tempk - 273.15
  pressesat = 610.8*exp(17.27*tempc/(tempc + 237.3)) !Pa
!pressvap = pressesat*rh !Pa
  pressvap = qair*press/(0.622 + qair) !Pa
!delta = 4098. * 610.8 * exp( 17.27*tempc / (tempc+237.3) )  / (tempc+237.3)**2.
!!Pa/K
!pressvap = pressesat*rh !Pa
  delta = 4098.*pressesat/(tempc + 237.3)**2. !Pa/K
  vpd = pressesat - pressvap
  lambda = (2.501 - 0.002361*tempc)*1.e6 !J/kg
  gamma = (cp*press)/(0.622*lambda)
  dens = press/(Rd*tempk*(1.+0.608*qair))

!make sure that hveg is not zero
  hveg = max(hhveg, 0.1)

  IF (nint(veg) .le. 1) then

    pet_w = (delta*rad + gamma*6.43*(1.+0.536*wind)*vpd/(24.*3600.))/(delta + gamma)
    pet_w = max(deltat*pet_w/lambda, 0.)

    pet_s = 0.
    pet_c = 0.
    pet_i = 0.

  ELSE
    pet_w = 0.

!!Radiation

!net radiation on the ground
    Rn_s = rad*exp(-0.5*lai)

!!!Resistances
!ra_a aerodynamic resistance from canopy to reference height
    za = hveg + 2. !ref. height

!roughness for a closed canopy z0c
    if (hveg .le. 1.) then
      z0c = 0.13*hveg
    elseif (hveg .gt. 1. .and. hveg .lt. 10.) then
      z0c = 0.139*hveg - 0.009*hveg**2.
    else
      z0c = 0.05*hveg
    end if

!mean drag coefficient for individual leafs
    if (hveg .eq. 0.) then
      c_d = 1.4e-3
    else
      c_d = (-1.+exp(0.909 - 3.03*z0c/hveg))**4./4.
    end if

!zero plane displacement height d0
    if (lai .ge. 4.) then
      d0 = max(hveg - z0c/0.3, 0.)
    else
      d0 = 1.1*hveg*log(1.+(c_d*lai)**0.25)
    end if

    if (d0 .gt. hveg) write (6, *) 'big problem!', d0, hveg, i, j
!if(i.eq.29.and.j.eq.19)write(6,*)'mirar lai,hveg,d0',lai,hveg,d0

!reference height
!   za = 10. + d0
!ground roughness length
    if (floodflag .eq. 0) then
      z0g = z0gr(nint(veg))
    else
      z0g = z0gr(1)
    end if

!roughness lengtt of canopy z0
    z0 = min(0.3*(hveg - d0), z0g + 0.3*hveg*(c_d*lai)**0.5)
    z0 = max(z0, z0g)

!friction  velocity ustar
!if(j.gt.8498)write(6,*)'mirar',i,j,wind,d0,z0
!    ustar = vk * wind / log( (za-d0)/z0 )
    ustar = vk*wind/log(10./z0)

!Eddy diffusion coefficient at the top of the canopy
    K_h = vk*ustar*(hveg - d0)

!eddy diff. decay constant for vegetation, n

    if (hveg .le. 1.) then
      n = 2.5
    elseif (hveg .gt. 1. .and. hveg .lt. 10.) then
      n = 2.306 + 0.194*hveg
    else
      n = 4.25
    end if

!preferred roughness length Z0
    Z0 = 0.13*hveg
!preferred zero plane displacement dp
    dp = 0.63*hveg
!and finally
    ra_a = log((za - d0)/(hveg - d0))/(vk*ustar) &
           + hveg*(exp(n*(1.-(Z0 + dp)/hveg)) - 1.)/(n*K_h)
!if(i.eq.29.and.j.eq.19)write(6,*)'mirar hveg,za,z0,ustar,K_h,n,Z0,dp,ra_a',hveg,za,z0,ustar,K_h,n,Z0,dp,log( (za-d0)/(hveg-d0) ),exp( n * ( 1. - (Z0+dp)/hveg ) )

!!ra_s aerodynamic resistance from soil to canopy
    ra_s = hveg*exp(n)*(exp(-n*z0g/hveg) - exp(-n*(Z0 + dp)/hveg))/(n*K_h)

!if(i.eq.29.and.j.eq.19)write(6,*)'mirar ra_a,ra_s',ra_a,ra_s
!!Bulk boundary layer resistance of canopy, ra_c
!uc, wind at canopy top
    uc = ustar*log((hveg - d0)/z0)/vk

!wleaf
    select case (nint(veg))
    case (4, 5, 13, 20, 21)
      wleaf = wmax(nint(veg))*(1.-exp(-0.6*lai))
    case default
      wleaf = wmax(nint(veg))
    end select

!rb
    rb = 100.*(wleaf/uc)**0.5/((1.-exp(-n/2.))*n)

!
    if (lai .gt. 0.1) then
      ra_c = rb*0.5/lai
    else
      ra_c = 0.
    end if

!!Bulk stomatal resistance of canopy rs_c
    frad = min(1., (0.004*rshort + 0.05)/(0.81*(1.+0.004*rshort)))
!this is how it was in the runs for PNAS
    fswp = 1.

    if (d0 .gt. 2.) then
      g_d = 0.0003
    else
      g_d = 0.
    end if

    fvpd = exp(-g_d*vpd)

    slai = 0.5*lai

    if (slai*frad*fswp*fvpd .eq. 0.) then
      rs_c = 5000.
    else
      rs_c = min(rl(nint(veg))/(slai*frad*fswp*fvpd), 5000.)
    end if

!if(i.eq.29.and.j.eq.19)write(6,*)'mirar uc,wleaf,rb,ra_c,rs_c',uc,wleaf,rb,ra_c,rs_c

!if(j.eq.801.and.(i.eq.295.or.i.eq.415))write(6,*)'mirar',i,j,tempk,rad,rshort,press,qair,wind,lai,veg,hhveg

!!Surface resistance of substrate soil rs_s
!now in extraction
!     if(floodflag.eq.0)then
!         rs_s = 500.
!     else
!         rs_s = 0.
!     endif

!!!!!!!!!!!!!!!!
    R_a = (delta + gamma)*ra_a
!      R_c = (delta + gamma) * ra_c + gamma*rs_c in the call
!      R_s = (delta + gamma) * ra_s + gamma*rs_s in the call

!      C_c = 1. / ( 1. + R_a*R_c / (R_s * (R_c+R_a) ) )
!      C_s = 1. / ( 1. + R_a*R_s / (R_c * (R_s+R_a) ) )

!      if(lai.lt.0.001)then
!             C_c=0.
!             C_s=1.
!      endif

!!PET
!pet_c = C_c * ( delta*rad + ( dens*cp*vpd - delta*ra_c*Rn_s ) / (ra_a+ra_c) ) / (delta + gamma*(1.+rs_c/(ra_a+ra_c)) )
!pet_c = max( deltat * pet_c / lambda , 0.)
    pet_c = (delta*rad + (dens*cp*vpd - delta*ra_c*Rn_s)/(ra_a + ra_c))
    pet_s = (delta*rad + (dens*cp*vpd - delta*ra_s*(rad - Rn_s))/(ra_a + ra_s))

!pet_s = C_s * ( delta*rad + ( dens*cp*vpd - delta*ra_s*(rad-Rn_s) ) / (ra_a+ra_s) ) / (delta + gamma*(1.+rs_s/(ra_a+ra_c)) )
!pet_s = max( deltat * pet_s / lambda , 0.)
!pet_s = ( delta*rad + ( dens*cp*vpd - delta*ra_s*(rad-Rn_s) ) / (ra_a+ra_s) ) / (delta + gamma*(1.+rs_s/(ra_a+ra_c)) )

!if(i.eq.29.and.j.eq.19)write(6,*)'mirar pet_s',pet_s,rad,Rn_s,ra_s,ra_a

!for PET from interception loss rs_c=rs_s = 0.
    R_c = (delta + gamma)*ra_c
    R_s = (delta + gamma)*ra_s
    C_c = 1./(1.+R_a*R_c/(R_s*(R_c + R_a)))
    if (lai .lt. 0.001) C_c = 0.

    pet_i = C_c*(delta*rad + (dens*cp*vpd - delta*ra_c*Rn_s)/(ra_a + ra_c))/(delta + gamma)
    pet_i = max(deltat*pet_i/lambda, 0.)

    pet_w = 0.

!if(j.eq.801.and.(i.eq.295.or.i.eq.415))write(6,*)'mirar 2',i,j,ra_a,ra_c,rs_c,ra_s,rs_s,pet_c
!if(i.eq.51.and.j.eq.51)write(6,*)'mirar pet_c,pet_s,pet_i',pet_c,pet_s,pet_i

    if ((pet_c .ne. pet_c) .or. (pet_s .ne. pet_s) .or. (pet_i .ne. pet_i)) then
      write (6, *) 'something wrong with pet', pet_c, pet_s, pet_i, ra_a, ra_c, rs_c, ra_s, rs_s
      write (6, *) 'forcings', i, j, tempk, rad, rshort, press, qair, wind, lai, veg, hhveg
    end if
  END IF

!pet = 3.*3600.* ( C_c*pet_c + C_s*pet_s + pet_w ) / lambda
end subroutine POTEVAP_Shutteworth_Wallace
