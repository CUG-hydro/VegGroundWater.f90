subroutine EXTRACTION(i, j, nzg, slz, dz, deltat, soiltxt, wtd, smoi, smoiwtd &
                      , delta, gamma, lambda, lai, ra_a, ra_c, rs_c_factor, R_a, R_s, petfactor_s, petfactor_c, pet_s &
                      , pet, watdef, dsmoi, dsmoideep &
                      , inactivedays, maxinactivedays, fieldcp, hhveg, fdepth, icefac)

  real, parameter :: potleaf = -153. !now equal to wilting point
!real, parameter :: potleaf = -102.!-204.  !for now constant
!real, parameter :: potleaf = -204.  !for now constant
  real, parameter :: potwilt = -153. !matric potential at wilting point
  real, parameter :: potfc = -3.366 !matric potential at field capacity
  integer :: i, j, nzg, nsoil, nsoil1, k, alarm, iwtd, kwtd, maxinactivedays, kroot
  real, dimension(nzg, nstyp) :: fieldcp
  integer :: soiltxt(2)
  integer, dimension(0:nzg + 1) :: inactivedays
  real, dimension(nzg + 1) :: slz
  integer, dimension(nzg) :: rootmask
  integer*1, dimension(nzg) :: icefac
  real, dimension(nzg) :: smoi, dz, dz2, vctr4, rootactivity, easy, dsmoi, maxwat
  real :: deltat, pet, transpwater, totwater, watdef, extract, kf, pot, toteasy, easydeep, dz3
  real :: wtd,smoiwtd,dzwtd,rootactivitydeep,dsmoideep,smoimin,maxeasy,soilfactor,fieldc,hveg,hhveg,zz,fdepth,psisat,smoisat,smoifc
   real :: delta,gamma,lambda,lai,ra_a,ra_c,rs_s,rs_c_factor,rs_c,R_a,R_s,petfactor_s,petfactor_c,pet_s,fswp,rootsmoi,rootfc,R_c,C_c,C_s

  hveg = 2.*hhveg/3.
!initialize
  easy = 0.
  easydeep = 0.
  dzwtd = 0.
  rootmask = 0
  dz2 = dz
  dz3 = 0.

!take water from layers
!    transpwater = pet * 1.e-3
!calculate where the water table is
  do k = 1, nzg
    if (wtd .lt. slz(k)) exit
  end do
  iwtd = k
  kwtd = k - 1

  if (kwtd .ge. 1 .and. kwtd .lt. nzg) dz2(kwtd) = slz(iwtd) - wtd

!calculate lowest layer of the root zone

  do k = 1, nzg
    if (inactivedays(k) .le. maxinactivedays) exit
  end do
  kroot = k - 1

  do k = max(kwtd, kroot, 1), nzg
!if(i.eq.50.and.j.eq.66)write(6,*)'mirar',k,inactivedays(k),inactivedays(k+1),wtd

!check if this layer has roots or can have roots growing from the layer above
!         if(inactivedays(k).gt.maxinactivedays.and.inactivedays(k+1).gt.maxinactivedays)cycle
    if (inactivedays(k) .le. maxinactivedays) rootmask(k) = 1

    vctr4(k) = 0.5*(slz(k) + slz(k + 1))

    if (slz(k) .lt. -0.30) then
      nsoil = soiltxt(1)
    else
      nsoil = soiltxt(2)
    end if

!calculte the easiness function for extraction for each layer
!      if(abs(theta_sat(nsoil)-smoi(k)).lt.1.e-6.and.k.ne.nzg)then
!           easy(k)=0.
!      if(smoi(k).le.slwilt(nsoil))then
!          easy(k)=0.
!      else

! calculate moisture potential
    smoisat = theta_sat(nsoil)*max(min(exp((vctr4(k) + 1.5)/fdepth), 1.), 0.1)
    psisat = slpots(nsoil)*min(max(exp(-(vctr4(k) + 1.5)/fdepth), 1.), 10.)
    pot = psisat*(smoisat/smoi(k))**slbs(nsoil)

    if (icefac(k) .eq. 0) then
      soilfactor = 1.
    else
      soilfactor = 0.
    end if

    easy(k) = max(-(potleaf - pot)*soilfactor/(hveg - vctr4(k)), 0.)
!      endif
  end do

  dsmoi = 0.
  dsmoideep = 0.
  watdef = 0.
!eliminate small root activity
!       where(easy.lt.0.01)easy=0.
!       if(easydeep.lt.0.01)easydeep=0.

!to grow roots anew, the layer has to be easiest to get water from than the
!current active layers with roots
  maxeasy = maxval(easy, rootmask == 1)

!eliminate small root activity
  where (easy .lt. 0.001*maxeasy) easy = 0.

!if(i.eq.50.and.j.eq.66)write(6,*)'mirar 2',maxeasy,easy
  do k = max(kroot, 1), nzg
    if (inactivedays(k) .gt. maxinactivedays .and. easy(k) .lt. maxeasy) easy(k) = 0.
  end do

  toteasy = sum(easy*dz2)
  if (toteasy .eq. 0.) then
!              watdef = transpwater
!              return
!         endif
    rootactivity = 0.
  else
    rootactivity = min(max((easy*dz2)/toteasy, 0.), 1.)
  end if
!if(i.eq.50.and.j.eq.66)write(6,*)'mirar 3',rootactivity

!eliminate small root activity
!         where(rootactivity.lt.0.01.and.rootmask==1)easy=0.
!         if(rootactivitydeep.lt.0.01)easydeep=0.
!recalculate
!         toteasy=sum(easy*dz) + easydeep*dzwtd
!         rootactivity = min ( max(  ( easy*dz ) /  toteasy  , 0. ), 1. )
!         rootactivitydeep = min ( max(  ( easydeep*dzwtd ) /  toteasy  , 0. ), 1. )

!if(i.eq.50.and.j.eq.66)write(6,*)'mirar 4',rootactivity

  do k = 1, nzg
    if (easy(k) .eq. 0.) then
      inactivedays(k) = inactivedays(k) + 1
    else
      inactivedays(k) = 0
    end if
  end do

  inactivedays = min(inactivedays, maxinactivedays + 1)

!if(i.eq.50.and.j.eq.66)write(6,*)'mirar 5',inactivedays
!         if(i.eq.89.and.j.eq.193)write(6,*)(rootactivity(k),k=1,nzg),sum(rootactivity)
!         if(i.eq.89.and.j.eq.193)write(6,*)(easy(k),k=1,nzg),sum(easy*dz),toteasy

  rootsmoi = 0.
  rootfc = 0.

  do k = max(kwtd, 1), nzg

    if (slz(k) .lt. -0.30) then
      nsoil = soiltxt(1)
    else
      nsoil = soiltxt(2)
    end if

    smoisat = theta_sat(nsoil)*max(min(exp((vctr4(k) + 1.5)/fdepth), 1.), 0.1)
    psisat = slpots(nsoil)*min(max(exp(-(vctr4(k) + 1.5)/fdepth), 1.), 10.)
    smoimin = smoisat*(psisat/potwilt)**(1./slbs(nsoil))
    smoifc = smoisat*(psisat/potfc)**(1./slbs(nsoil))

!         smoimin = max(smoimin,slwilt(nsoil))
    maxwat(k) = max((smoi(k) - smoimin)*dz(k), 0.)  !max water than can be taken from a layer

    rootsmoi = rootsmoi + max(rootactivity(k)*(smoi(k) - smoimin), 0.)
    rootfc = rootfc + max(rootactivity(k)*(smoifc - smoimin), 0.)
  end do

  if (rootsmoi .le. 0) then
    fswp = 0. ! 水分限制因子
  elseif (rootsmoi/rootfc .le. 1) then
    fswp = rootsmoi/rootfc
  else
    fswp = 1.
  end if

  if (fswp .eq. 0.) then
    rs_c = 5000. ! stomatal resistance
  else
    rs_c = min(rs_c_factor/fswp, 5000.)
  end if

  nsoil = soiltxt(2)
  rs_s = 33.5 + 3.5*(theta_sat(nsoil)/smoi(nzg))**2.38 !soil resistance

  R_c = (delta + gamma)*ra_c + gamma*rs_c
  R_s = R_s + gamma*rs_s

  C_c = 1./(1.+R_a*R_c/(R_s*(R_c + R_a)))
  C_s = 1./(1.+R_a*R_s/(R_c*(R_s + R_a)))

!if(i.eq.29.and.j.eq.19)write(6,*)'mirar C_c,C_s,lai',C_c,C_c,lai
  if (lai .lt. 0.001) then
    C_c = 0.
!             C_s=1.
  end if

!calculate transpiration and soil evaporation. Both depend on stomatal
!resistence, thats why the final step has to be computed here
  pet = C_c*petfactor_c/(delta + gamma*(1.+rs_c/(ra_a + ra_c)))
  pet = max(deltat*pet/lambda, 0.)

  pet_s = C_s*petfactor_s/(delta + gamma*(1.+rs_s/(ra_a + ra_c)))
  pet_s = max(deltat*pet_s/lambda, 0.)

!if(i.eq.29.and.j.eq.19)write(6,*)'mirar final',pet_s,rs_s,ra_a,ra_c

  transpwater = pet*1.e-3

  if (toteasy .eq. 0.) then
    watdef = transpwater
    return
  end if

  do k = max(kwtd, 1), nzg
!calculate hyd. conductivity
!         kf =   Ksat(nsoil) * (smoi(k)  / theta_sat(nsoil)) ** (2. * slbs(nsoil) + 3.)
!         maxwat = min( maxwat , kf*deltat )

!extract water
    extract = max(rootactivity(k)*transpwater, 0.)   !water to be extracted from this layer this timestep

    if (extract .le. maxwat(k)) then
      dsmoi(k) = extract
    else
      dsmoi(k) = maxwat(k)
      watdef = watdef + (extract - maxwat(k))
    end if
  end do

  dsmoi = max(dsmoi, 0.)
  if (abs(watdef - transpwater) .gt. 1.e9) write (6, *) 'algo no esta bien', i, j, transpwater*1.e3, watdef*1.e3

!now total rootactiviy is dsmoi/totwater normalized by soil layer depth, return dsmoi (total water taken from each layer) to do calculation later and update soil moisture

!       smoi = smoi - dsmoi/dz
end subroutine extraction
