# Translated from soilfluxes.f90

function soilfluxes(i, j, nzg, freedrain, dtll, slz, dz, soiltxt, smoiwtd, transp, transpdeep,
                    smoi, wtd, rech, deeprech, precip, pet_s, et_s, runoff, flux, fdepth, qlat,
                    qlatflux, qrf, qrfcorrect, flood, icefactor)
    vctr2 = similar(dz)
    vctr4 = similar(dz)
    vctr5 = similar(dz)
    vctr6 = similar(dz)
    kfmid = zeros(eltype(dz), nzg)
    diffmid = zeros(eltype(dz), nzg)
    aa = zeros(eltype(dz), nzg)
    bb = zeros(eltype(dz), nzg)
    cc = zeros(eltype(dz), nzg)
    rr = zeros(eltype(dz), nzg)
    vt3di = zeros(eltype(dz), nzg + 1)

    for k in 1:nzg
        vctr2[k] = 1.0 / dz[k]
        vctr4[k] = 0.5 * (slz[k] + slz[k + 1])
    end
    for k in 2:nzg
        vctr5[k] = vctr4[k] - vctr4[k - 1]
        vctr6[k] = 1.0 / vctr5[k]
    end

    rech = 0.0
    runoff = 0.0
    qgw = qlat - qrf
    vt3di[nzg + 1] = (-precip + pet_s) * 1.0e-3 - flood

    iwtd = if freedrain == 0
        k = 1
        for kk in 1:nzg
            k = kk
            wtd < slz[kk] && break
        end
        k
    else
        0
    end

    k = max(iwtd - 1, 1)
    qlatflux[k] += qgw

    for k in max(iwtd - 1, 2):nzg
        wgpmid = smoi[k] + (smoi[k] - smoi[k - 1]) * (slz[k] - vctr4[k]) * vctr6[k]
        nsoil = slz[k] < -0.30 ? soiltxt[1] : soiltxt[2]

        hydcon = Ksat(nsoil) * clamp(exp((slz[k] + 1.5) / fdepth), 0.1, 1.0)
        smoisat = theta_sat(nsoil) * clamp(exp((slz[k] + 1.5) / fdepth), 0.1, 1.0)
        psisat = slpots(nsoil) * min(max(exp(-(slz[k] + 1.5) / fdepth), 1.0), 10.0)
        wgpmid = min(wgpmid, smoisat)
        icefac = icefactor[k] == 0 ? 1.0 : 0.0

        kfmid[k] = icefac * hydcon * (wgpmid / smoisat)^(2.0 * slbs(nsoil) + 3.0)
        diffmid[k] = -icefac * (hydcon * psisat * slbs(nsoil) / smoisat) * (wgpmid / smoisat)^(slbs(nsoil) + 2.0)
    end

    for k in max(iwtd - 2, 2):nzg
        aa[k] = diffmid[k] * vctr6[k]
        cc[k] = diffmid[k + 1] * vctr6[k + 1]
        bb[k] = -(aa[k] + cc[k] + dz[k] / dtll)
        rr[k] = -smoi[k] * dz[k] / dtll - kfmid[k + 1] + kfmid[k] + transp[k] / dtll
        if k == iwtd - 1
            rr[k] -= qgw / dtll
        end
    end

    aa[nzg] = diffmid[nzg] * vctr6[nzg]
    bb[nzg] = -aa[nzg] - dz[nzg] / dtll
    rr[nzg] = vt3di[nzg + 1] / dtll - smoi[nzg] * dz[nzg] / dtll + kfmid[nzg] + transp[nzg] / dtll
    if iwtd - 1 == nzg
        rr[nzg] -= qgw / dtll
    end

    if freedrain != 1
        if iwtd <= 3
            aa[1] = 0.0
            cc[1] = diffmid[2] * vctr6[2]
            bb[1] = -(aa[1] + cc[1] + dz[1] / dtll)
            rr[1] = -smoi[1] * dz[1] / dtll - kfmid[2] + transp[1] / dtll
            if iwtd <= 2
                rr[1] -= qgw / dtll
            end
        else
            for k in 1:max(iwtd - 3, 1)
                aa[k] = 0.0
                cc[k] = 0.0
                bb[k] = 1.0
                rr[k] = smoi[k]
            end
        end
    else
        nsoil = soiltxt[1]
        hydcon = Ksat(nsoil) * max(min(exp((slz[1] + 1.5) / fdepth), 1.0), 0.1)
        smoisat = theta_sat(nsoil) * max(min(exp((slz[1] + 1.5) / fdepth), 1.0), 0.1)
        kfmid[1] = hydcon * (smoi[1] / smoisat)^(2.0 * slbs(nsoil) + 3.0)

        aa[1] = 0.0
        cc[1] = diffmid[2] * vctr6[2]
        bb[1] = -(cc[1] + dz[1] / dtll)
        rr[1] = -smoi[1] * dz[1] / dtll - kfmid[2] + kfmid[1] + transp[1] / dtll
    end

    tridag(aa, bb, cc, rr, smoi, nzg)

    for k in 2:nzg
        vt3di[k] = (-aa[k] * (smoi[k] - smoi[k - 1]) - kfmid[k]) * dtll
    end
    vt3di[1] = 0.0

    for k in 1:nzg
        nsoil = slz[k] < -0.30 ? soiltxt[1] : soiltxt[2]
        smoisat = theta_sat(nsoil) * max(min(exp((vctr4[k] + 1.5) / fdepth), 1.0), 0.1)

        if smoi[k] > smoisat
            dsmoi = max((smoi[k] - smoisat) * dz[k], 0.0)
            if k < nzg
                smoi[k + 1] += dsmoi * vctr2[k + 1]
                vt3di[k + 1] += dsmoi
            else
                runoff = dsmoi
            end
            smoi[k] = smoisat
        end
    end

    for k in 1:(nzg - 1)
        nsoil = slz[k] < -0.30 ? soiltxt[1] : soiltxt[2]
        smoicp = theta_cp(nsoil) * max(min(exp((vctr4[k] + 1.5) / fdepth), 1.0), 0.1)

        if smoi[k] < smoicp
            dsmoi = max((smoicp - smoi[k]) * dz[k], 0.0)
            smoi[k + 1] -= dsmoi * vctr2[k + 1]
            vt3di[k + 1] -= dsmoi
            smoi[k] = smoicp
        end
    end

    k = nzg
    nsoil = slz[k] < -0.30 ? soiltxt[1] : soiltxt[2]
    smoicp = theta_cp(nsoil) * max(min(exp((vctr4[k] + 1.5) / fdepth), 1.0), 0.1)

    if smoi[k] < smoicp
        dsmoi = max((smoicp - smoi[k]) * dz[k], 0.0)
        if vt3di[k + 1] > dsmoi
            et_s = max(0.0, pet_s - dsmoi * 1.0e3)
            smoi[k] = smoicp
        else
            et_s = max(0.0, pet_s - max(vt3di[k + 1], 0.0) * 1.0e3)
            smoi[k] += max(vt3di[k + 1], 0.0) / dz[k]
            dsmoi = max((smoicp - smoi[k]) * dz[k], 0.0)
            smoi[k - 1] -= dsmoi * vctr2[k - 1]
            vt3di[k] += dsmoi
            smoi[k] = smoicp

            for kk in (nzg - 1):-1:1
                nsoil = slz[kk] < -0.30 ? soiltxt[1] : soiltxt[2]
                smoicp = theta_cp(nsoil) * max(min(exp((vctr4[kk] + 1.5) / fdepth), 1.0), 0.1)
                if smoi[kk] < smoicp
                    dsmoi = max((smoicp - smoi[kk]) * dz[kk], 0.0)
                    if kk > 1
                        smoi[kk - 1] -= dsmoi * vctr2[kk - 1]
                    end
                    vt3di[kk] += dsmoi
                    smoi[kk] = smoicp
                end
            end
        end
    else
        et_s = pet_s
    end

    if vt3di[1] > 0.0
        qrfcorrect = -min(vt3di[1], max(qrf, 0.0))
    else
        qrfcorrect = 0.0
    end

    flux .+= vt3di

    if freedrain == 1
        rech = vt3di[1]
        smoiwtd -= vt3di[1]
    end

    return smoiwtd, wtd, rech, et_s, runoff, qrfcorrect
end
