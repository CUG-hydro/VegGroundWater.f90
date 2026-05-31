# Translated from lateralflow.f90

function lateral(imax, jmax, js, je, soiltxt, wtd, qlat, fdepth, topo, landmask, deltat, area, lats, dxy)
  qlat .= 0.0
  klat = zeros(eltype(wtd), size(wtd))
  for j in js:je
    for i in 1:imax
      nsoil = soiltxt[1, i, j]
      klat[i, j] = Ksat(nsoil) * klatfactor(nsoil)
    end
  end
  lateralflow4(imax, jmax, js, je, wtd, qlat, fdepth, topo, landmask, deltat, area, klat, lats, dxy)
end

function lateralflow(imax, jmax, js, je, wtd, qlat, fdepth, topo, landmask, deltat, area, klat)
  fangle = sqrt(tan(pi4 / 32.0)) / (2.0 * sqrt(2.0))
  kcell = zeros(eltype(wtd), size(wtd))
  head = zeros(eltype(wtd), size(wtd))

  # lateral flow calculation
  for j in max(js, 1):min(je, jmax)
    for i in 1:imax
      if fdepth[i, j] < 1.0e-6
        kcell[i, j] = 0.0
      elseif wtd[i, j] < -1.5
        kcell[i, j] = fdepth[i, j] * klat[i, j] * exp((wtd[i, j] + 1.5) / fdepth[i, j])
      else
        kcell[i, j] = klat[i, j] * (wtd[i, j] + 1.5 + fdepth[i, j])
      end
      head[i, j] = topo[i, j] + wtd[i, j]
    end
  end

  for j in js+1:je-1
    for i in 2:imax-1
      if landmask[i, j] == 1
        q = 0.0

        q += (kcell[i-1, j+1] + kcell[i, j]) * (head[i-1, j+1] - head[i, j]) / sqrt(2.0)
        q += (kcell[i-1, j] + kcell[i, j]) * (head[i-1, j] - head[i, j])
        q += (kcell[i-1, j-1] + kcell[i, j]) * (head[i-1, j-1] - head[i, j]) / sqrt(2.0)
        q += (kcell[i, j+1] + kcell[i, j]) * (head[i, j+1] - head[i, j])
        q += (kcell[i, j-1] + kcell[i, j]) * (head[i, j-1] - head[i, j])
        q += (kcell[i+1, j+1] + kcell[i, j]) * (head[i+1, j+1] - head[i, j]) / sqrt(2.0)
        q += (kcell[i+1, j] + kcell[i, j]) * (head[i+1, j] - head[i, j])
        q += (kcell[i+1, j-1] + kcell[i, j]) * (head[i+1, j-1] - head[i, j]) / sqrt(2.0)

        qlat[i, j] = fangle * q * deltat / area[i, j]
      end
    end
  end
end

function lateralflow4(imax, jmax, js, je, wtd, qlat, fdepth, topo, landmask, deltat, area, klat, xlat, dxy)
  d2r = π / 180.0
  kcell = zeros(eltype(wtd), size(wtd))
  head = zeros(eltype(wtd), size(wtd))

  # lateral flow calculation
  for j in max(js, 1):min(je, jmax)
    for i in 1:imax
      if fdepth[i, j] < 1.0e-6
        kcell[i, j] = 0.0
      elseif wtd[i, j] < -1.5
        kcell[i, j] = fdepth[i, j] * klat[i, j] * exp((wtd[i, j] + 1.5) / fdepth[i, j])
      else
        kcell[i, j] = klat[i, j] * (wtd[i, j] + 1.5 + fdepth[i, j])
      end
      head[i, j] = topo[i, j] + wtd[i, j]
    end
  end

  for j in js+1:je-1
    for i in 2:imax-1
      if landmask[i, j] == 1
        q = 0.0
        # north
        q += (kcell[i, j+1] + kcell[i, j]) * (head[i, j+1] - head[i, j]) * cos(d2r * (xlat[i, j] + 0.5 * dxy))
        # south
        q += (kcell[i, j-1] + kcell[i, j]) * (head[i, j-1] - head[i, j]) * cos(d2r * (xlat[i, j] - 0.5 * dxy))
        # west
        q += (kcell[i-1, j] + kcell[i, j]) * (head[i-1, j] - head[i, j]) / cos(d2r * xlat[i, j])
        # east
        q += (kcell[i+1, j] + kcell[i, j]) * (head[i+1, j] - head[i, j]) / cos(d2r * xlat[i, j])

        qlat[i, j] = 0.5 * q * deltat / area[i, j]
      end
    end
  end
end
