# Translated from moveqrf.f90

function moveqrf(imax, js, je, fd, qrf, area, width)
  qrfextra = zeros(eltype(qrf), size(qrf))

  for j in js+1:je-1
    for i in 2:imax-1
      if fd[i, j] > 0 && width[i, j] < 1.0
        iout, jout = flowdir(imax, js, je, fd, i, j)
        qrfextra[iout, jout] += qrf[i, j] * area[i, j] / area[iout, jout]
        qrf[i, j] = 0.0
      end
    end
  end

  qrf .+= qrfextra
end
