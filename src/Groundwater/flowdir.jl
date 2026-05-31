# Translated from flowdir.f90

function flowdir(imax, js, je, fd, ii, jj)
  j = if fd[ii, jj] in (2, 4, 8)
    jj - 1
  elseif fd[ii, jj] in (1, 16)
    jj
  elseif fd[ii, jj] in (32, 64, 128)
    jj + 1
  else
    0
  end

  i = if fd[ii, jj] in (128, 1, 2)
    ii + 1
  elseif fd[ii, jj] in (4, 64)
    ii
  elseif fd[ii, jj] in (8, 16, 32)
    ii - 1
  else
    0
  end

  return i, j
end
