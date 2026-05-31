# Translated from module_rootdepth.f90

const nvtyp = 30
const nstyp = 13

const THETA_SAT = [0.395, 0.410, 0.435, 0.485, 0.451, 0.420,
                   0.477, 0.476, 0.426, 0.492, 0.482, 0.863, 0.476]
const THETA_CP = [0.050, 0.052, 0.092, 0.170, 0.125, 0.148,
                  0.195, 0.235, 0.202, 0.257, 0.268, 0.195, 0.235]
const SLBS = [4.05, 4.38, 4.9, 5.3, 5.39, 7.12, 7.75, 8.52,
              10.4, 10.4, 11.4, 7.75, 8.52]
const KSAT = [0.000176, 0.0001563, 0.00003467,
              0.0000072, 0.00000695, 0.0000063,
              0.0000017, 0.00000245, 0.000002167,
              0.000001033, 0.000001283, 0.0000080, 0.000005787]
const SLPOTS = [-0.121, -0.090, -0.218, -0.786, -0.478, -0.299,
                -0.356, -0.630, -0.153, -0.490, -0.405, -0.356, -0.630]
const KLATFACTOR = [2.0, 3.0, 4.0, 10.0, 12.0, 14.0, 20.0, 24.0, 28.0, 40.0, 48.0, 48.0, 48.0]
const SLWILT = zeros(Float64, nstyp)

# helper accessors aligned with previous Julia translations
θ_sat(i) = THETA_SAT[i]
θ_cp(i) = THETA_CP[i]
B(i) = SLBS[i]
ψ_sat(i) = SLPOTS[i]
Ksat(i) = KSAT[i]
slbs(i) = SLBS[i]
theta_sat(i) = THETA_SAT[i]
theta_cp(i) = THETA_CP[i]
slpots(i) = SLPOTS[i]

include(joinpath(@__DIR__, "tridag.jl"))
include(joinpath(@__DIR__, "update_shallow_wtd.jl"))
include(joinpath(@__DIR__, "update_wtb_qsat.jl"))
include(joinpath(@__DIR__, "extraction.jl"))
include(joinpath(@__DIR__, "potevap.jl"))
include(joinpath(@__DIR__, "soilfluxes.jl"))
include(joinpath(@__DIR__, "interception.jl"))
include(joinpath(@__DIR__, "init_soil.jl"))

function rootdepth(freedrain, imax, js, je, nzg, slz, dz, deltat, landmask, veg, hveg, soiltxt, wind, temp, qair, press, netrad, rshort,
                   lai, precip, qsrun, smoi, smoieq, smoiwtd, wtd, waterdeficit, watext, watextdeep, rech, deeprech,
                   et_s, et_i, et_c, intercepstore, ppacum, pppendepth, pppendepthold,
                   qlat, qlatsum, qsprings, inactivedays, maxinactivedays, fieldcp, fdepth, steps, floodheight,
                   qrf, delsfcwat, icefactor, wtdflux, et_s_daily, et_c_daily, transptop, infilk)

    minpprate = 0.01
    icefac = zeros(Int8, max(nzg, 40))

    for j in (js + 1):(je - 1)
        for i in 1:imax
            landmask[i, j] == 0 && continue

            if nzg >= 40
                icefac[26:40] .= icefactor[i, j, 26:40]
            else
                icefac .= 0
            end
            floodflag = floodheight[i, j] > 0.05 ? 1 : 0

            delta, gamma, lambda, ra_a, ra_c, rs_c, R_a, R_s, petfactor_s, petfactor_c, petstep_w, petstep_i =
                potevap_shutteworth_wallace(i, j, deltat, temp[i, j], netrad[i, j], rshort[i, j], press[i, j], qair[i, j],
                                            wind[i, j], lai[i, j], veg[i, j], hveg[i, j],
                                            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                            0.0, 0.0, 0.0, 0.0, floodflag)

            et_s[i, j] += petstep_w
            if floodflag == 1 && round(Int, veg[i, j]) <= 1
                delsfcwat[i, j] -= petstep_w * 1.0e-3
            end

            round(Int, veg[i, j]) <= 1 && continue

            intercepstore[i, j], ppdrip, etstep_i = interception(minpprate, precip[i, j], lai[i, j], intercepstore[i, j], 0.0, petstep_i, 0.0)
            et_i[i, j] += etstep_i

            ppdrip /= steps
            floodstep = floodheight[i, j] / steps
            qlatstep = qlat[i, j] / steps
            qrfstep = qrf[i, j] / steps

            flux = zeros(Float64, nzg + 1)
            qlatflux = zeros(Float64, nzg + 2)
            wtdold = wtd[i, j]

            for _itime in 1:round(Int, steps)
                dsmoi = zeros(Float64, nzg)
                watdef = 0.0
                dsmoideep = 0.0
                petstep_s = 0.0
                petstep_c = 0.0

                extraction(i, j, nzg, slz, dz, deltat / steps, soiltxt[1, i, j], wtd[i, j], smoi[:, i, j],
                           delta, gamma, lambda, lai[i, j], ra_a, ra_c, rs_c, R_a, R_s, petfactor_s, petfactor_c,
                           petstep_s, petstep_c, watdef, dsmoi, inactivedays[:, i, j], maxinactivedays, fieldcp,
                           hveg[i, j], fdepth[i, j], icefac[1:nzg])

                et_c[i, j] += petstep_c - watdef * 1.0e3
                waterdeficit[i, j] += watdef * 1.0e3
                watext[:, i, j] .+= dsmoi .* 1.0e3
                transptop[i, j] += dsmoi[nzg] * 1.0e3
                et_c_daily[i, j] += petstep_c - watdef * 1.0e3

                smoiwtd[i, j], wtd[i, j], rechstep, etstep_s, runoff, qrfcorrect = soilfluxes(
                    i, j, nzg, freedrain, deltat / steps, slz, dz, soiltxt[:, i, j], smoiwtd[i, j], dsmoi, dsmoideep,
                    smoi[:, i, j], wtd[i, j], 0.0, deeprech[i, j], ppdrip, petstep_s, 0.0, 0.0, flux, fdepth[i, j],
                    qlatstep, qlatflux, qrfstep, 0.0, floodstep, icefac[1:nzg])

                delsfcwat[i, j] -= max(floodstep - runoff, 0.0)
                qsrun[i, j] += max(runoff - floodstep, 0.0)
                rech[i, j] += rechstep * 1.0e3
                et_s[i, j] += etstep_s
                et_s_daily[i, j] += etstep_s
                ppacum[i, j] += ppdrip
                qrf[i, j] += qrfcorrect

                wtd[i, j], rechstep = update_shallow_wtd(i, j, nzg, freedrain, slz, dz, soiltxt[:, i, j], smoieq[:, i, j],
                                                         smoiwtd[i, j], smoi[:, i, j], wtd[i, j], 0.0, fdepth[i, j])
                rech[i, j] += rechstep * 1.0e3
                qlatsum[i, j] += qlatstep
            end

            infilkstep = nzg + 1
            pppendepthstep = 0.0
            flux[nzg + 1] = -1.0
            for k in nzg:-1:0
                if k <= nzg - 2 && pppendepthold[i, j] >= k + 3
                    break
                end
                if flux[k + 1] < -0.333e-5
                    if k == 0
                        if -flux[1] > -qlatflux[k + 1] && pppendepthstep > slz[1]
                            pppendepthstep = slz[1]
                            infilkstep = 1
                        end
                    elseif -flux[k + 1] + flux[k] > -qlatflux[k + 1] + 0.0 && pppendepthstep > slz[k + 1]
                        pppendepthstep = slz[k + 1]
                        infilkstep = k + 1
                    end
                end
            end

            pppendepthold[i, j] = infilkstep
            if pppendepth[i, j] > pppendepthstep
                pppendepth[i, j] = pppendepthstep
            end
            if slz[max(infilkstep - 1, 1)] <= wtdold
                wtdflux[i, j] -= flux[infilkstep] * 1.0e3
            end
            if infilk[i, j] > infilkstep
                infilk[i, j] = infilkstep
            end
        end
    end
    return nothing
end
