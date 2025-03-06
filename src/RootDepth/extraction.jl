function extraction(i, j, nzg, slz, dz, deltat, soiltxt, wtd, smoi, smoiwtd,
  delta, gamma, lambda, lai, ra_a, ra_c, rs_c_factor, R_a, R_s, petfactor_s, petfactor_c, pet_s,
  pet, watdef, dsmoi, dsmoideep,
  inactivedays, maxinactivedays, fieldcp, hhveg, fdepth, icefac)

  # Constants
  potleaf = -153.0  # now equal to wilting point
  potwilt = -153.0  # matric potential at wilting point
  potfc = -3.366    # matric potential at field capacity

  # Initialize arrays and variables
  easy = zeros(Float64, nzg)
  easydeep = 0.0
  dzwtd = 0.0
  rootmask = zeros(Int, nzg)
  dz2 = copy(dz)
  dz3 = 0.0
  rootactivity = zeros(Float64, nzg)
  vctr4 = zeros(Float64, nzg)
  maxwat = zeros(Float64, nzg)

  hveg = 2.0 * hhveg / 3.0

  # Calculate where the water table is
  k = 1
  for k = 1:nzg
    if wtd < slz[k]
      break
    end
  end
  iwtd = k
  kwtd = k - 1

  if kwtd >= 1 && kwtd < nzg
    dz2[kwtd] = slz[iwtd] - wtd
  end

  # Calculate lowest layer of the root zone
  k = 1
  for k = 1:nzg
    if inactivedays[k] <= maxinactivedays
      break
    end
  end
  kroot = k - 1

  for k = max(kwtd, kroot, 1):nzg
    if inactivedays[k] <= maxinactivedays
      rootmask[k] = 1
    end

    vctr4[k] = 0.5 * (slz[k] + slz[k+1])

    if slz[k] < -0.30
      nsoil = soiltxt[1]
    else
      nsoil = soiltxt[2]
    end

    # Calculate moisture potential
    smoisat = theta_sat(nsoil) * max(min(exp((vctr4[k] + 1.5) / fdepth), 1.0), 0.1)
    psisat = slpots(nsoil) * min(max(exp(-(vctr4[k] + 1.5) / fdepth), 1.0), 10.0)
    pot = psisat * (smoisat / smoi[k])^slbs(nsoil)

    if icefac[k] == 0
      soilfactor = 1.0
    else
      soilfactor = 0.0
    end

    easy[k] = max(-(potleaf - pot) * soilfactor / (hveg - vctr4[k]), 0.0)
  end

  dsmoi .= 0.0
  dsmoideep = 0.0
  watdef = 0.0

  # Find maximum easiness among active root layers
  maxeasy = maximum(easy[findall(rootmask .== 1)])

  # Eliminate small root activity
  for k = 1:nzg
    if easy[k] < 0.001 * maxeasy
      easy[k] = 0.0
    end
  end

  # Check inactive days
  for k = max(kroot, 1):nzg
    if inactivedays[k] > maxinactivedays && easy[k] < maxeasy
      easy[k] = 0.0
    end
  end

  toteasy = sum(easy .* dz2)
  if toteasy == 0.0
    rootactivity .= 0.0
  else
    rootactivity .= min.(max.((easy .* dz2) ./ toteasy, 0.0), 1.0)
  end

  # Update inactive days
  for k = 1:nzg
    if easy[k] == 0.0
      inactivedays[k] += 1
    else
      inactivedays[k] = 0
    end
  end

  inactivedays .= min.(inactivedays, maxinactivedays + 1)

  rootsmoi = 0.0
  rootfc = 0.0

  for k = max(kwtd, 1):nzg
    if slz[k] < -0.30
      nsoil = soiltxt[1]
    else
      nsoil = soiltxt[2]
    end

    smoisat = theta_sat(nsoil) * max(min(exp((vctr4[k] + 1.5) / fdepth), 1.0), 0.1)
    psisat = slpots(nsoil) * min(max(exp(-(vctr4[k] + 1.5) / fdepth), 1.0), 10.0)
    smoimin = smoisat * (psisat / potwilt)^(1.0 / slbs(nsoil))
    smoifc = smoisat * (psisat / potfc)^(1.0 / slbs(nsoil))

    maxwat[k] = max((smoi[k] - smoimin) * dz[k], 0.0)  # max water that can be taken from a layer

    rootsmoi += max(rootactivity[k] * (smoi[k] - smoimin), 0.0)
    rootfc += max(rootactivity[k] * (smoifc - smoimin), 0.0)
  end

  if rootsmoi <= 0
    fswp = 0.0
  elseif rootsmoi / rootfc <= 1
    fswp = rootsmoi / rootfc
  else
    fswp = 1.0
  end

  if fswp == 0.0
    rs_c = 5000.0
  else
    rs_c = min(rs_c_factor / fswp, 5000.0)
  end

  nsoil = soiltxt[2]
  rs_s = 33.5 + 3.5 * (theta_sat(nsoil) / smoi[nzg])^2.38

  R_c = (delta + gamma) * ra_c + gamma * rs_c
  R_s = R_s + gamma * rs_s

  C_c = 1.0 / (1.0 + R_a * R_c / (R_s * (R_c + R_a)))
  C_s = 1.0 / (1.0 + R_a * R_s / (R_c * (R_s + R_a)))

  if lai < 0.001
    C_c = 0.0
  end

  # Calculate transpiration and soil evaporation
  pet = C_c * petfactor_c / (delta + gamma * (1.0 + rs_c / (ra_a + ra_c)))
  pet = max(deltat * pet / lambda, 0.0)

  pet_s = C_s * petfactor_s / (delta + gamma * (1.0 + rs_s / (ra_a + ra_c)))
  pet_s = max(deltat * pet_s / lambda, 0.0)

  transpwater = pet * 1.0e-3

  if toteasy == 0.0
    watdef = transpwater
    return
  end

  for k = max(kwtd, 1):nzg
    # Extract water
    extract = max(rootactivity[k] * transpwater, 0.0)  # water to be extracted from this layer this timestep

    if extract <= maxwat[k]
      dsmoi[k] = extract
    else
      dsmoi[k] = maxwat[k]
      watdef += (extract - maxwat[k])
    end
  end

  dsmoi .= max.(dsmoi, 0.0)
  if abs(watdef - transpwater) > 1.0e9
    println("algo no esta bien ", i, " ", j, " ", transpwater * 1.0e3, " ", watdef * 1.0e3)
  end

  # Now total rootactiviy is dsmoi/totwater normalized by soil layer depth
  # Return dsmoi (total water taken from each layer) to do calculation later and update soil moisture
end
