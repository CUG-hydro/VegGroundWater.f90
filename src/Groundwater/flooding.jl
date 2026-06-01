# Translated from flooding.f90

function flooding(imax, js, je, deltat, fd, bfd, topo, area, riverwidth, riverlength, riverdepth, floodheight, delsfcwat)
  dflood = zeros(eltype(floodheight), size(floodheight))

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

  delsfcwat .+= dflood
end
