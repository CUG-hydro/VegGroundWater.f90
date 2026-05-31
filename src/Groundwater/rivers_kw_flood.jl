# Translated from rivers_kw_flood.f90

function rivers_kw_flood(imax, js, je, deltat, dtlr, fd, bfd, qnew, qs, qrf, delsfcwat,
                         slope, depth, width, length, maxdepth, area, riverarea, floodarea, riverchannel,
                         qmean, floodheight, topo)

  q = similar(qnew)
  qin = zeros(eltype(qnew), size(qnew))
  qext = zeros(eltype(qnew), size(qnew))

  for j in js+1:je-1
    for i in 2:imax-1
      if fd[i, j] != 0
        qext[i, j] = (qrf[i, j] + qs[i, j] + delsfcwat[i, j]) / deltat * area[i, j]
      end
    end
  end

  q .= qnew
  qin .= 0.0

  for j in js:je
    for i in 1:imax
      if fd[i, j] > 0
        i1, j1 = flowdir(imax, js, je, fd, i, j)
        if i1 > 1 && i1 < imax && j1 > js && j1 < je
          qin[i1, j1] += q[i, j]
        end
      end
    end
  end

  for j in js+1:je-1
    for i in 2:imax-1
      if fd[i, j] != 0
        # calculate total inflow into cell i j
        dsnew = qin[i, j] - q[i, j]

        # Taquari
        if i == 4498 && j == 4535
          dsnew += q[4499, 4534] / 4.0
        end
        if i == 4498 && j == 4534
          dsnew -= q[4499, 4534] / 4.0
        end
        if i == 4464 && j == 4536
          dsnew += q[4465, 4535] / 2.0
        end
        if i == 4465 && j == 4534
          dsnew -= q[4465, 4535] / 2.0
        end
        if i == 4346 && j == 4560
          dsnew += q[4346, 4561] / 3.0
        end
        if i == 4345 && j == 4561
          dsnew -= q[4346, 4561] / 3.0
        end
        if i == 4444 && j == 4551
          dsnew += q[4444, 4552] / 3.0
        end
        if i == 4443 && j == 4553
          dsnew -= q[4444, 4552] / 3.0
        end
        if i == 4350 && j == 4497
          dsnew += q[4352, 4496] / 3.0
        end
        if i == 4351 && j == 4496
          dsnew -= q[4352, 4496] / 3.0
        end

        # Sao Lourenco
        if i == 4439 && j == 4772
          dsnew += q[4440, 4773] / 2.0
        end
        if i == 4440 && j == 4772
          dsnew -= q[4440, 4773] / 2.0
        end
        if i == 4400 && j == 4685
          dsnew += q[4401, 4685] / 2.0
        end
        if i == 4400 && j == 4684
          dsnew -= q[4401, 4685] / 2.0
        end
        if i == 4418 && j == 4688
          dsnew += q[4418, 4689] / 5.0
        end
        if i == 4417 && j == 4689
          dsnew -= q[4418, 4689] / 5.0
        end
        if i == 4367 && j == 4698
          dsnew += q[4368, 4699] / 6.0
        end
        if i == 4368 && j == 4698
          dsnew -= q[4368, 4699] / 6.0
        end
        if i == 4363 && j == 4667
          dsnew += q[4364, 4668] / 6.0
        end
        if i == 4364 && j == 4667
          dsnew -= q[4364, 4668] / 6.0
        end
        if i == 4475 && j == 4718
          dsnew += q[4475, 4717] / 6.0
        end
        if i == 4474 && j == 4717
          dsnew -= q[4475, 4717] / 6.0
        end

        snew = depth[i, j] * riverarea[i, j] + floodheight[i, j] * floodarea[i, j] + (dsnew + qext[i, j]) * dtlr

        # now redistribute water between river channel and floodplain and calculate new riverdepth and floodheight
        if isnan(snew)
          println("problem with snew ", i, " ", j, " ", dsnew, " ", floodheight[i, j], " ", depth[i, j], " ", qext[i, j])
        end

        if snew >= riverchannel[i, j]
          floodheight[i, j] = (snew - riverchannel[i, j]) / max(area[i, j], riverarea[i, j])
          depth[i, j] = floodheight[i, j] + maxdepth[i, j]
        else
          floodheight[i, j] = 0.0
          if riverarea[i, j] > 0.0
            depth[i, j] = snew / riverarea[i, j]
          else
            depth[i, j] = 0.0
          end
        end

        if isnan(depth[i, j])
          println("problem with depth ", i, " ", j, " ", qrf[i, j], " ", qs[i, j], " ", delsfcwat[i, j], " ", qnew[i, j], " ", floodheight[i, j])
        end
      end
    end
  end

  for j in js+1:je-1
    for i in 2:imax-1
      flowwidth = width[i, j]
      if fd[i, j] != 0
        if width[i, j] * depth[i, j] > 1.0e-9 && fd[i, j] > 0
          # calculate speed from Manning's formula
          aa = depth[i, j] * width[i, j] / (2.0 * depth[i, j] + width[i, j])
          if floodheight[i, j] > 0.05
            i1, j1 = flowdir(imax, js, je, fd, i, j)
            waterelevij = topo[i, j] - maxdepth[i, j] + depth[i, j]
            waterelevi1j1 = topo[i1, j1] - maxdepth[i1, j1] + max(depth[i1, j1], 0.0)
            slopefor = (waterelevij - waterelevi1j1) / (0.5 * (length[i, j] + length[i1, j1]))
            if bfd[i, j] > 0
              i2, j2 = flowdir(imax, js, je, bfd, i, j)
              waterelevi2j2 = topo[i2, j2] - maxdepth[i2, j2] + max(depth[i2, j2], 0.0)
              slopeback = (waterelevi2j2 - waterelevij) / (0.5 * (length[i2, j2] + length[i, j]))
              slopeinst = 0.5 * (slopefor + slopeback)
            else
              slopeinst = slopefor
            end
            slopeinst = 0.25 * slopeinst + 0.75 * slope[i, j]
            if slopeinst < 0.0
              slopeinst = slope[i, j]
            end
          else
            slopeinst = slope[i, j]
          end
          speed = (aa^(2.0 / 3.0)) * sqrt(slopeinst) / 0.03
          speed = max(min(speed, length[i, j] / dtlr), 0.01)
        else
          speed = 0.0
        end
        # now calculate the new q
        qnew[i, j] = speed * depth[i, j] * flowwidth
      else
        qnew[i, j] = 0.0
      end
    end
  end

  qmean .+= qnew .* dtlr
end
