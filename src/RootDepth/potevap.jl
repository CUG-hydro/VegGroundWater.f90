# Translated from potevap.f90

function potevap_priestly_taylor(i, j, tempk, rad, presshp, pet)
    cp = 1013.0e-6
    tempc = tempk - 273.15
    presskp = presshp * 0.1
    rad = rad * 24.0 * 3600.0 * 1.0e-6

    alpha = 1.26
    delta = 0.2 * (0.00738 * tempc + 0.8072)^7.0 - 0.000116
    λ = 2.501 - 0.002361 * tempc
    γ = (cp * presskp) / (0.622 * λ)

    pet = alpha * rad * delta / (delta + γ)
    pet = pet / λ
    return pet
end

function potevap_penman_monteith(i, j, tempk, rad, rshort, press, qair, wind, lai, veg, hveg, pet)
    cp = 1013.0
    vk = 0.41
    Rd = 287.0
    rl = [150.0, 150.0, 500.0, 500.0, 175.0, 240.0, 110.0, 100.0, 250.0, 150.0,
          80.0, 225.0, 225.0, 250.0, 180.0, 180.0, 240.0, 500.0, 240.0, 500.0,
          175.0, 250.0, 250.0, 175.0, 225.0, 150.0, 110.0, 180.0, 250.0, 250.0]

    tempc = tempk - 273.15
    pressesat = 610.8 * exp(17.27 * tempc / (tempc + 237.3))
    pressvap = qair * press / (0.622 + qair)
    delta = 4098.0 * pressesat / (tempc + 237.3)^2
    vpd = pressesat - pressvap
    λ = (2.501 - 0.002361 * tempc) * 1.0e6
    γ = (cp * press) / (0.622 * λ)
    dens = press / (Rd * tempk * (1.0 + 0.608 * qair))

    zm = max(10.0, hveg)
    hdisp = 0.0
    z0m = 0.1 * hveg
    zh = max(2.0, hveg)
    z0h = 0.1 * z0m
    ra = log((zm - hdisp) / z0m) * log((zh - hdisp) / z0h) / (vk^2 * wind)

    frad = min(1.0, (0.004 * rshort + 0.05) / (0.81 * (1.0 + 0.004 * rshort)))
    fswp = 1.0
    g_d = hdisp > 2.0 ? 0.0003 : 0.0
    fvpd = exp(-g_d * vpd)
    slai = 0.5 * lai

    if slai * frad * fswp * fvpd == 0.0
        rs = 5000.0
    else
        veg_idx = clamp(round(Int, veg), 1, length(rl))
        rs = min(rl[veg_idx] / (slai * frad * fswp * fvpd), 5000.0)
    end

    pet = (delta * rad + dens * cp * vpd / ra) / (delta + γ * (1.0 + rs / ra))
    pet = 3.0 * 3600.0 * pet / λ
    return pet
end

function potevap_shutteworth_wallace(i, j, deltat, tempk, rad, rshort, press, qair, wind, lai, veg, hhveg,
                                     delta, gamma, lambda, ra_a, ra_c, rs_c, R_a, R_s,
                                     pet_s, pet_c, pet_w, pet_i, floodflag)
    cp = 1013.0
    vk = 0.41
    Rd = 287.0

    rl = [150.0, 150.0, 500.0, 500.0, 175.0, 240.0, 110.0, 100.0, 250.0, 150.0,
          80.0, 225.0, 225.0, 250.0, 180.0, 180.0, 240.0, 500.0, 240.0, 500.0,
          175.0, 250.0, 250.0, 175.0, 225.0, 150.0, 110.0, 180.0, 250.0, 250.0]
    bioparms = [
        0.001 0.0;
        0.001 0.0;
        0.001 0.0;
        0.02 0.001;
        0.02 0.001;
        0.02 0.08;
        0.02 0.05;
        0.01 0.01;
        0.01 0.01;
        0.001 0.01;
        0.01 0.01;
        0.01 0.01;
        0.02 0.01;
        0.02 0.01;
        0.02 0.04;
        0.005 0.01;
        0.005 0.01;
        0.01 0.01;
        0.01 0.001;
        0.02 0.05;
        0.02 0.001;
        0.02 0.08;
        0.01 0.01;
        0.02 0.04;
        0.02 0.01;
        0.02 0.01;
        0.01 0.01;
        0.005 0.01;
        0.001 0.01;
        0.02 0.0
    ]

    z0gr = bioparms[:, 1]
    wmax = bioparms[:, 2]

    tempc = tempk - 273.15
    pressesat = 610.8 * exp(17.27 * tempc / (tempc + 237.3))
    pressvap = qair * press / (0.622 + qair)
    delta = 4098.0 * pressesat / (tempc + 237.3)^2
    vpd = pressesat - pressvap
    lambda = (2.501 - 0.002361 * tempc) * 1.0e6
    gamma = (cp * press) / (0.622 * lambda)
    dens = press / (Rd * tempk * (1.0 + 0.608 * qair))

    hveg = max(hhveg, 0.1)

    if round(Int, veg) <= 1
        pet_w = (delta * rad + gamma * 6.43 * (1.0 + 0.536 * wind) * vpd / (24.0 * 3600.0)) / (delta + gamma)
        pet_w = max(deltat * pet_w / lambda, 0.0)
        pet_s = 0.0
        pet_c = 0.0
        pet_i = 0.0
    else
        pet_w = 0.0
        Rn_s = rad * exp(-0.5 * lai)
        za = hveg + 2.0

        if hveg <= 1.0
            z0c = 0.13 * hveg
        elseif hveg < 10.0
            z0c = 0.139 * hveg - 0.009 * hveg^2
        else
            z0c = 0.05 * hveg
        end

        c_d = hveg == 0.0 ? 1.4e-3 : (-1.0 + exp(0.909 - 3.03 * z0c / hveg))^4 / 4.0
        d0 = lai >= 4.0 ? max(hveg - z0c / 0.3, 0.0) : 1.1 * hveg * log(1.0 + (c_d * lai)^0.25)

        idx = clamp(round(Int, veg), 1, length(rl))
        z0g = floodflag == 0 ? z0gr[idx] : z0gr[1]
        z0 = min(0.3 * (hveg - d0), z0g + 0.3 * hveg * sqrt(c_d * lai))
        z0 = max(z0, z0g)
        ustar = vk * wind / log(10.0 / z0)
        K_h = vk * ustar * (hveg - d0)

        if hveg <= 1.0
            n = 2.5
        elseif hveg < 10.0
            n = 2.306 + 0.194 * hveg
        else
            n = 4.25
        end

        Z0 = 0.13 * hveg
        dp = 0.63 * hveg
        ra_a = log((za - d0) / (hveg - d0)) / (vk * ustar) + hveg * (exp(n * (1.0 - (Z0 + dp) / hveg)) - 1.0) / (n * K_h)
        ra_s = hveg * exp(n) * (exp(-n * z0g / hveg) - exp(-n * (Z0 + dp) / hveg)) / (n * K_h)

        uc = ustar * log((hveg - d0) / z0) / vk
        wleaf = idx in (4, 5, 13, 20, 21) ? wmax[idx] * (1.0 - exp(-0.6 * lai)) : wmax[idx]
        rb = 100.0 * sqrt(wleaf / uc) / ((1.0 - exp(-n / 2.0)) * n)
        ra_c = lai > 0.1 ? rb * 0.5 / lai : 0.0

        frad = min(1.0, (0.004 * rshort + 0.05) / (0.81 * (1.0 + 0.004 * rshort)))
        fswp = 1.0
        g_d = d0 > 2.0 ? 0.0003 : 0.0
        fvpd = exp(-g_d * vpd)
        slai = 0.5 * lai
        rs_c = slai * frad * fswp * fvpd == 0.0 ? 5000.0 : min(rl[idx] / (slai * frad * fswp * fvpd), 5000.0)

        R_a = (delta + gamma) * ra_a
        pet_c = delta * rad + (dens * cp * vpd - delta * ra_c * Rn_s) / (ra_a + ra_c)
        pet_s = delta * rad + (dens * cp * vpd - delta * ra_s * (rad - Rn_s)) / (ra_a + ra_s)

        R_c = (delta + gamma) * ra_c
        R_s = (delta + gamma) * ra_s
        C_c = 1.0 / (1.0 + R_a * R_c / (R_s * (R_c + R_a)))
        lai < 0.001 && (C_c = 0.0)

        pet_i = C_c * (delta * rad + (dens * cp * vpd - delta * ra_c * Rn_s) / (ra_a + ra_c)) / (delta + gamma)
        pet_i = max(deltat * pet_i / lambda, 0.0)
        pet_w = 0.0

        if isnan(pet_c) || isnan(pet_s) || isnan(pet_i)
            println("Warning: NaN detected in potential evapotranspiration: pet_c=", pet_c, " pet_s=", pet_s, " pet_i=", pet_i,
                    " ra_a=", ra_a, " ra_c=", ra_c, " rs_c=", rs_c, " ra_s=", ra_s)
            println("forcings: i=", i, " j=", j, " tempk=", tempk, " rad=", rad, " rshort=", rshort,
                    " press=", press, " qair=", qair, " wind=", wind, " lai=", lai, " veg=", veg, " hhveg=", hhveg)
        end
    end

    return delta, gamma, lambda, ra_a, ra_c, rs_c, R_a, R_s, pet_s, pet_c, pet_w, pet_i
end
