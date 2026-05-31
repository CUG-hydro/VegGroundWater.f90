# Translated from GW2river.f90

function gw2river(imax, js, je, nzg, slz, deltat, soiltxt, landmask, wtd, maxdepth, riverdepth, width, length, area, fdepth, qrf)
  soilwatercap = 0.0
  qrf .= 0.0

  for j in js+1:je-1
    for i in 2:imax-1
      (landmask[i, j] == 0 || width[i, j] == 0.0) && continue
      rdepth = max(riverdepth[i, j], 0.0)

      nsoil = soiltxt[2, i, j]
      riversurface = -(maxdepth[i, j] - rdepth)
      riversurface >= 0.0 && continue

      # Fan 2007, Eq. 6
      T2 = Ksat(nsoil) * clamp(exp((-maxdepth[i, j] + 1.5) / fdepth[i, j]), 0.1, 1.0)
      rcond = width[i, j] * length[i, j] * T2

      if wtd[i, j] > riversurface
        qrf[i, j] = rcond * (wtd[i, j] - riversurface) * deltat / area[i, j]
      elseif wtd[i, j] > -maxdepth[i, j]
        soilwatercap = -rcond * (wtd[i, j] - riversurface) * (deltat / area[i, j])
        qrf[i, j] = -max(min(soilwatercap, riverdepth[i, j]), 0.0) * min(width[i, j] * length[i, j] / area[i, j], 1.0)
      else
        # Water table below riverbed: disconnected from channel, infiltration at Ksat only
        frac = min(width[i, j] * length[i, j] / area[i, j], 1.0)
        qrf[i, j] = -max(min(Ksat(nsoil) * deltat, rdepth), 0.0) * frac
      end
    end
  end
end
