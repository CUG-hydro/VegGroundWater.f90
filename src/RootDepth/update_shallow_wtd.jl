# Translated from update_shallow_wtd.f90

function update_shallow_wtd(i, j, nzg, freedrain, slz, dz, soiltxt, smoieq, smoiwtd, smoi, wtd, rech, fdepth)
    rech = 0.0
    flag = 0
    vctr4 = zeros(eltype(slz), nzg)

    for k in 1:nzg
        vctr4[k] = 0.5 * (slz[k] + slz[k + 1])
    end

    k = 1
    for kk in 1:nzg
        k = kk
        wtd < slz[kk] && break
    end
    iwtd = k
    kwtd = iwtd - 1

    if kwtd > 0
        wtdold = wtd
        nsoil = slz[kwtd] < -0.30 ? soiltxt[1] : soiltxt[2]

        if kwtd > 1
            smoisat = theta_sat(nsoil) * max(min(exp((vctr4[kwtd - 1] + 1.5) / fdepth), 1.0), 0.1)
            if wtd < slz[kwtd] + 0.01 && smoi[kwtd - 1] < smoisat
                flag = 1
            end
        end

        smoisat = theta_sat(nsoil) * max(min(exp((vctr4[kwtd] + 1.5) / fdepth), 1.0), 0.1)

        if smoi[kwtd] > smoieq[kwtd] && flag == 0
            if smoi[kwtd] == smoisat
                wtd = slz[iwtd]
                rech = (wtdold - wtd) * (smoisat - smoieq[kwtd])
                iwtd += 1
                kwtd += 1

                if kwtd <= nzg && smoi[kwtd] > smoieq[kwtd]
                    wtdold = wtd
                    nsoil = slz[kwtd] < -0.30 ? soiltxt[1] : soiltxt[2]
                    smoisat = theta_sat(nsoil) * max(min(exp((vctr4[kwtd] + 1.5) / fdepth), 1.0), 0.1)
                    wtd = min((smoi[kwtd] * dz[kwtd] - smoieq[kwtd] * slz[iwtd] + smoisat * slz[kwtd]) /
                              (smoisat - smoieq[kwtd]), slz[iwtd])
                    rech += (wtdold - wtd) * (smoisat - smoieq[kwtd])
                end
            else
                wtd = min((smoi[kwtd] * dz[kwtd] - smoieq[kwtd] * slz[iwtd] + smoisat * slz[kwtd]) /
                          (smoisat - smoieq[kwtd]), slz[iwtd])
                rech = (wtdold - wtd) * (smoisat - smoieq[kwtd])
            end
        else
            wtd = slz[kwtd]
            rech = (wtdold - wtd) * (smoisat - smoieq[kwtd])
            kwtd -= 1
            iwtd -= 1

            if kwtd >= 1
                wtdold = wtd
                nsoil = slz[kwtd] < -0.30 ? soiltxt[1] : soiltxt[2]
                smoisat = theta_sat(nsoil) * max(min(exp((vctr4[kwtd] + 1.5) / fdepth), 1.0), 0.1)

                if smoi[kwtd] > smoieq[kwtd]
                    wtd = min((smoi[kwtd] * dz[kwtd] - smoieq[kwtd] * slz[iwtd] + smoisat * slz[kwtd]) /
                              (smoisat - smoieq[kwtd]), slz[iwtd])
                else
                    wtd = slz[kwtd]
                end
                rech += (wtdold - wtd) * (smoisat - smoieq[kwtd])
            end
        end
    end

    if wtd < slz[1]
        println("Warning: water table depth below minimum layer: wtd=", wtd, " i=", i, " j=", j)
    end

    return wtd, rech
end
