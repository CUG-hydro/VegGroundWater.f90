# Translated from module_wtable.f90

module ModuleWTable

include(joinpath(@__DIR__, "..", "RootDepth", "module_rootdepth.jl"))

const pi4 = 4π

include("flowdir.jl")
include("flooding.jl")
include("GW2river.jl")
include("lateralflow.jl")
include("moveqrf.jl")
include("rivers_kw_flood.jl")
include("update_deep_wtb.jl")
include("update_wtd.jl")

function wtable(imax, jmax, js, je, nzg, slz, dz, area, soiltxt, wtd, bottomflux, rech, qslat, fdepth, topo, landmask, deltat,
                smoi, smoieq, smoiwtd, qsprings)

  # Calculate lateral flow
  qlat = zeros(eltype(wtd), size(wtd))
  klat = zeros(eltype(wtd), size(wtd))
  for j in js:je
    for i in 1:imax
      nsoil = soiltxt[1, i, j]
      klat[i, j] = Ksat(nsoil) * klatfactor(nsoil)
    end
  end

  lateralflow(imax, jmax, js, je, wtd, qlat, fdepth, topo, landmask, deltat, area, klat)

  qslat .+= qlat .* 1.0e3

  # now calculate deep recharge
  deeprech = zeros(eltype(wtd), size(wtd))

  for j in js+1:je-1
    for i in 1:imax
      if landmask[i, j] == 1
        if wtd[i, j] < slz[1] - dz[1]
          # calculate k for drainage
          nsoil = soiltxt[1, i, j]
          wgpmid = 0.5 * (smoiwtd[i, j] + theta_sat(nsoil))
          kfup = Ksat(nsoil) * (wgpmid / theta_sat(nsoil))^(2.0 * slbs(nsoil) + 3.0)

          # now calculate moisture potential
          vt3dbdw = slpots(nsoil) * (theta_sat(nsoil) / smoiwtd[i, j])^slbs(nsoil)

          # and now flux (=recharge)
          deeprech[i, j] = deltat * kfup * ((slpots(nsoil) - vt3dbdw) / (slz[1] - wtd[i, j]) - 1.0)

          # now update smoiwtd
          newwgp = smoiwtd[i, j] + (deeprech[i, j] - bottomflux[i, j]) / (slz[1] - wtd[i, j])
          if newwgp < theta_cp(nsoil)
            deeprech[i, j] += (theta_cp(nsoil) - newwgp) * (slz[1] - wtd[i, j])
            newwgp = theta_cp(nsoil)
          end
          if newwgp > theta_sat(nsoil)
            deeprech[i, j] -= (theta_sat(nsoil) - newwgp) * (slz[1] - wtd[i, j])
            newwgp = theta_sat(nsoil)
          end

          smoiwtd[i, j] = newwgp
          rech[i, j] += deeprech[i, j] * 1.0e3
        end
      end
    end
  end

  bottomflux .= 0.0

  # Now update water table and soil moisture
  for j in js+1:je-1
    for i in 1:imax
      if landmask[i, j] == 1
        # Total groundwater balance in the cell
        totwater = qlat[i, j] - deeprech[i, j]
        if isnan(qlat[i, j])
          println("gran problema!", wtd[i, j], qlat[i, j], i, j)
        end

        qspring = 0.0
        wtd[i, j], qspring, smoiwtd[i, j] = update_wtd(
          nzg, slz, dz, wtd[i, j], qspring, totwater, smoi[:, i, j], smoieq[:, i, j], soiltxt[:, i, j], smoiwtd[i, j]
        )

        qsprings[i, j] += qspring * 1.0e3
      end
    end
  end
end

end # module ModuleWTable
