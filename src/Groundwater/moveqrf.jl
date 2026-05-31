# Translated from moveqrf.f90

function moveqrf(imax, js, je, fd, qrf, area, width)
  qrfextra = zeros(eltype(qrf), size(qrf))

  for j in js+1:je-1
    for i in 2:imax-1
      if fd[i, j] > 0
        if width[i, j] < 1.0
          iout, jout = flowdir(imax, js, je, fd, i, j)
          qrfextra[iout, jout] += qrf[i, j] * area[i, j] / area[iout, jout]
          qrf[i, j] = 0.0
        end
      end
    end
  end

  if numtasks > 1
    qrfextra2 = copy(qrfextra)
    reqsu2, reqsd2, reqru2, reqrd2 = sendbordersflood(imax, js, je, qrfextra2)
  end

  # make sure that the borders are received before calculating anything
  if pid == 1
    MPI_wait(reqru2, status, ierr)
  elseif pid == numtasks - 2
    MPI_wait(reqrd2, status, ierr)
  elseif pid > 1 && pid < numtasks - 2
    MPI_wait(reqru2, status, ierr)
    MPI_wait(reqrd2, status, ierr)
  end

  if pid == 1
    qrfextra[1:imax, je-1] .+= qrfextra2[1:imax, je-1]
  elseif pid == numtasks - 2
    qrfextra[1:imax, js+1] .+= qrfextra2[1:imax, js+1]
  elseif pid > 1 && pid < numtasks - 2
    qrfextra[1:imax, js+1] .+= qrfextra2[1:imax, js+1]
    qrfextra[1:imax, je-1] .+= qrfextra2[1:imax, je-1]
  end

  # change qrf
  qrf .+= qrfextra

  # before changing qrfextra make sure that the borders have been received
  if pid == 1
    MPI_wait(reqsu2, status, ierr)
  elseif pid == numtasks - 2
    MPI_wait(reqsd2, status, ierr)
  elseif pid > 1 && pid < numtasks - 2
    MPI_wait(reqsu2, status, ierr)
    MPI_wait(reqsd2, status, ierr)
  end
end
