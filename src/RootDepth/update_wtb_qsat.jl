# Translated from update_wtb_qsat.f90

function update_wtb_qlat(nzg, slz, dz, wtd, qspring, qlat, smoi, smoieq, soiltextures, smoiwtd, qlatflux, fdepth)
    vctr4 = zeros(eltype(slz), nzg)
    for k in 1:nzg
        vctr4[k] = 0.5 * (slz[k] + slz[k + 1])
    end

    soiltxt = [slz[k] < -0.3 ? soiltextures[1] : soiltextures[2] for k in 1:nzg]
    qspring = 0.0
    totwater = qlat

    if totwater > 0.0
        k = 2
        for kk in 2:nzg
            k = kk
            wtd < slz[kk] && break
        end
        iwtd = k
        kwtd = iwtd - 1
        nsoil = soiltxt[kwtd]
        smoisat = theta_sat(nsoil) * max(min(exp((vctr4[kwtd] + 1.5) / fdepth), 1.0), 0.1)
        maxwatup = dz[kwtd] * (smoisat - smoi[kwtd])

        if totwater <= maxwatup
            smoi[kwtd] += totwater / dz[kwtd]
            qlatflux[kwtd] += totwater
            smoi[kwtd] = min(smoi[kwtd], smoisat)
            if smoi[kwtd] > smoieq[kwtd]
                wtd = min((smoi[kwtd] * dz[kwtd] - smoieq[kwtd] * slz[iwtd] + smoisat * slz[kwtd]) /
                          (smoisat - smoieq[kwtd]), slz[iwtd])
            end
            totwater = 0.0
        else
            smoi[kwtd] = smoisat
            qlatflux[kwtd] += maxwatup
            totwater -= maxwatup
            for k in iwtd:nzg
                wtd = slz[k]
                iwtd = k + 1
                nsoil = soiltxt[k]
                smoisat = theta_sat(nsoil) * max(min(exp((vctr4[k] + 1.5) / fdepth), 1.0), 0.1)
                maxwatup = dz[k] * (smoisat - smoi[k])
                if totwater <= maxwatup
                    smoi[k] += totwater / dz[k]
                    qlatflux[k] += totwater
                    smoi[k] = min(smoi[k], smoisat)
                    if smoi[k] > smoieq[k]
                        wtd = min((smoi[k] * dz[k] - smoieq[k] * slz[iwtd] + smoisat * slz[k]) /
                                  (smoisat - smoieq[k]), slz[iwtd])
                    end
                    totwater = 0.0
                    break
                else
                    smoi[k] = smoisat
                    qlatflux[k] += maxwatup
                    totwater -= maxwatup
                end
            end
            if totwater > 0.0
                wtd = slz[nzg + 1]
            end
        end
        qspring = totwater

    elseif totwater < 0.0
        k = 2
        for kk in 2:nzg
            k = kk
            wtd < slz[kk] && break
        end
        iwtd = k

        for kwtd in (iwtd - 1):-1:1
            nsoil = soiltxt[kwtd]
            smoisat = theta_sat(nsoil) * max(min(exp((vctr4[kwtd] + 1.5) / fdepth), 1.0), 0.1)
            maxwatdw = dz[kwtd] * (smoi[kwtd] - smoieq[kwtd])

            if -totwater <= maxwatdw
                smoi[kwtd] += totwater / dz[kwtd]
                qlatflux[kwtd] += totwater
                if smoi[kwtd] > smoieq[kwtd]
                    wtd = (smoi[kwtd] * dz[kwtd] - smoieq[kwtd] * slz[iwtd] + smoisat * slz[kwtd]) /
                          (smoisat - smoieq[kwtd])
                else
                    wtd = slz[kwtd]
                    iwtd -= 1
                end
                totwater = 0.0
                break
            else
                wtd = slz[kwtd]
                iwtd -= 1
                if maxwatdw >= 0.0
                    smoi[kwtd] = smoieq[kwtd]
                    qlatflux[kwtd] += maxwatdw
                    totwater += maxwatdw
                end
            end
        end

        if iwtd == 1 && totwater < 0.0
            nsoil = soiltxt[1]
            smoisat = theta_sat(nsoil) * max(min(exp((vctr4[1] + 1.5) / fdepth), 1.0), 0.1)
            smoi[1] += totwater / dz[1]
            qlatflux[1] += totwater
            wtd = max((smoi[1] * dz[1] - smoieq[1] * slz[2] + smoisat * slz[1]) /
                      (smoisat - smoieq[1]), slz[1])
        end
        qspring = 0.0
    end

    return wtd, qspring
end
