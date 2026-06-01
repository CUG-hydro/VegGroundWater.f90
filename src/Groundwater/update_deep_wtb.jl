function update_deep_wtb(imax, jmax, nzg, slz, dz, soiltxt, wtd, bottomflux, rech, qslat, qlat, landmask, deltat, smoi, smoieq, smoiwtd, qsprings)
  deeprech = zeros(size(wtd))

  for j in 2:jmax-1
    for i in 1:imax
      if landmask[i, j] == 1
        if wtd[i, j] < slz[1] - dz[1]
          nsoil = soiltxt[1, i, j]
          wgpmid = 0.5 * (smoiwtd[i, j] + theta_sat(nsoil))
          kfup = Ksat(nsoil) * (wgpmid / theta_sat(nsoil))^(2 * slbs(nsoil) + 3)
          vt3dbdw = slpots(nsoil) * (theta_sat(nsoil) / smoiwtd[i, j])^slbs(nsoil)
          deeprech[i, j] = deltat * kfup * ((slpots(nsoil) - vt3dbdw) / (slz[1] - wtd[i, j]) - 1)
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
          rech[i, j] += deeprech[i, j] * 1e3
        end
      end
    end
  end

  bottomflux .= 0

  for j in 2:jmax-1
    for i in 1:imax
      if landmask[i, j] == 1
        totwater = qlat[i, j] - qslat[i, j] - deeprech[i, j]
        qspring = 0.0
        wtd[i, j], qspring, smoiwtd[i, j] = update_wtd(
          nzg, slz, dz, wtd[i, j], qspring, totwater, smoi[:, i, j], smoieq[:, i, j], soiltxt[:, i, j], smoiwtd[i, j]
        )
        qsprings[i, j] += qspring * 1e3
      end
    end
  end

end
