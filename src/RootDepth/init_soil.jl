# Translated from init_soil.f90

function init_soil_depth_clm(nzg, slz, dz)
    slz2 = zeros(eltype(slz), nzg + 1)
    dz2 = zeros(eltype(dz), nzg)
    vctr4 = zeros(eltype(dz), nzg + 1)

    for k in 1:(nzg + 1)
        vctr4[k] = 0.025 * (exp(0.5 * (float(k) - 0.5)) - 1.0)
    end

    for k in 2:(nzg - 1)
        dz2[k] = 0.5 * (vctr4[k + 1] - vctr4[k - 1])
    end
    dz2[1] = 0.5 * (vctr4[1] + vctr4[2])
    dz2[nzg] = vctr4[nzg] - vctr4[nzg - 1]

    for k in 1:nzg
        slz2[k] = 0.5 * (vctr4[k] + vctr4[k + 1])
    end
    slz2[nzg] = vctr4[nzg] + 0.5 * dz2[nzg]

    for k in 1:nzg
        kk = nzg - k + 1
        slz[k] = -slz2[kk]
        dz[k] = dz2[kk]
    end
    slz[nzg + 1] = 0.0
    return nothing
end

function init_soil_depth(nzg, slz, dz)
    dz2 = [0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.2,
           0.3, 0.3, 0.3, 0.3, 0.4, 0.4, 0.4, 0.5, 0.5, 0.6,
           0.7, 0.7, 0.8, 0.9, 1.0, 1.0, 1.2, 1.2, 1.5, 1.5,
           2.0, 2.0, 3.0, 6.0, 11.0, 20.0, 50.0, 100.0, 250.0, 540.0]

    slz[nzg + 1] = 0.0
    for k in nzg:-1:1
        dz[k] = dz2[nzg - k + 1]
        slz[k] = slz[k + 1] - dz[k]
    end
    return nothing
end

khyd(smoi, i) = Ksat(i) * (smoi / theta_sat(i))^(2.0 * slbs(i) + 3.0)

function init_soil_param(fieldcp, nzg)
    potwilt = -153.0
    for nsoil in 1:nstyp
        SLWILT[nsoil] = theta_sat(nsoil) * (slpots(nsoil) / potwilt)^(1.0 / slbs(nsoil))
    end
    return nothing
end
