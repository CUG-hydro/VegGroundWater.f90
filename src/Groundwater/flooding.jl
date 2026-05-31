# Translated from flooding.f90

function flooding(imax, js, je, deltat, fd, bfd, topo, area, riverwidth, riverlength, riverdepth, floodheight, delsfcwat)
  ntsplit = 1

  dflood = zeros(eltype(floodheight), size(floodheight))
  dflood2 = zeros(eltype(floodheight), size(floodheight))

  for _ in 1:ntsplit
    # communicate flood water height to neighboring cells
    if numtasks > 1
      reqsu, reqsd, reqru, reqrd = sendborders(imax, js, je, floodheight)
    end

    # make sure that the borders are received before calculating anything
    if pid == 1
      MPI_wait(reqru, status, ierr)
    elseif pid == numtasks - 2
      MPI_wait(reqrd, status, ierr)
    elseif pid > 1 && pid < numtasks - 2
      MPI_wait(reqru, status, ierr)
      MPI_wait(reqrd, status, ierr)
    end

    for j in js+1:je-1
      for i in 2:imax-1
        fd[i, j] == 0 && continue
        if floodheight[i, j] > 0.05
          # find the lowest elevation neighbour that is not along the main river channel
          dhmax = 0.0
          ilow = i
          jlow = j
          for jj in j-1:j+1
            for ii in i-1:i+1
              (ii == i && jj == j) && continue

              dh = floodheight[i, j] + topo[i, j] - (floodheight[ii, jj] + topo[ii, jj])
              if ii != i && jj != j
                dh /= sqrt(2.0)
              end
              if dh > dhmax
                ilow = ii
                jlow = jj
                dhmax = dh
              end
            end
          end

          # now flood the lowest elevation neighbour
          if dhmax > 0.0
            i1, j1 = flowdir(imax, js, je, fd, i, j)

            dtotal = floodheight[i, j] + floodheight[ilow, jlow]
            dij = max(floodheight[i, j] - max(0.5 * (topo[ilow, jlow] - topo[i, j] + dtotal), 0.0), 0.0)
            if ilow == i1 && jlow == j1
              # the flow along the river channel is taken care of by the river routine
              dij = max(dij - (riverwidth[i, j] * floodheight[i, j] * riverlength[i, j]) / area[i, j], 0.0)
            end
            if delsfcwat[i, j] < 0.0
              dij = max(min(dij, floodheight[i, j] + delsfcwat[i, j]), 0.0)
            end
            dflood[i, j] -= dij
            dflood[ilow, jlow] += dij * area[i, j] / area[ilow, jlow]
          end
        end
      end
    end

    if numtasks > 1
      dflood2 .= dflood
      reqsu2, reqsd2, reqru2, reqrd2 = sendbordersflood(imax, js, je, dflood2)
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
      dflood[1:imax, je-1] .+= dflood2[1:imax, je-1]
    elseif pid == numtasks - 2
      dflood[1:imax, js+1] .+= dflood2[1:imax, js+1]
    elseif pid > 1 && pid < numtasks - 2
      dflood[1:imax, js+1] .+= dflood2[1:imax, js+1]
      dflood[1:imax, je-1] .+= dflood2[1:imax, je-1]
    end

    # before changing floodheight make sure that the borders have been received
    if pid == 1
      MPI_wait(reqsu, status, ierr)
    elseif pid == numtasks - 2
      MPI_wait(reqsd, status, ierr)
    elseif pid > 1 && pid < numtasks - 2
      MPI_wait(reqsu, status, ierr)
      MPI_wait(reqsd, status, ierr)
    end

    # update floodheight and riverdepth
    delsfcwat .+= dflood

    # before changing dflood make sure that the borders have been received
    if pid == 1
      MPI_wait(reqsu2, status, ierr)
    elseif pid == numtasks - 2
      MPI_wait(reqsd2, status, ierr)
    elseif pid > 1 && pid < numtasks - 2
      MPI_wait(reqsu2, status, ierr)
      MPI_wait(reqsd2, status, ierr)
    end
  end
end
