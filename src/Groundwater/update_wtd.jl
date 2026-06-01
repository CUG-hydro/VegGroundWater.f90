# j   : SM
# j-1 : GW & SM 混合层
function find_jwt(zwt, z₋ₕ)
  nzg = length(z₋ₕ)
  for k in 2:nzg
    zwt < z₋ₕ[k] && return k
  end
  return nzg
end


# qspring = 0.0
function update_wtd(nzg, z₋ₕ, dz, zwt, qspring, ∑, θ, θ_eq, soiltextures, θ_wtd)
  soiltxt = [z₋ₕ[i] < -0.3 ? soiltextures[1] : soiltextures[2] for i in 1:(nzg+1)]

  # Case 1: Water table going up (totwater > 0)
  if ∑ > 0.0
    if zwt >= z₋ₕ[1]
      iwtd = find_jwt(zwt, z₋ₕ)  # SM
      kwtd = iwtd - 1            # GW & SM 混合层
      nsoil = soiltxt[kwtd]

      # Maximum water that fits in the layer
      _layer = dz[kwtd] * (θ_sat(nsoil) - θ[kwtd])   # 最大补给量, dz应该修改，这里写的不好。

      if ∑ <= _layer # 补给量仅够补给一层
        θ[kwtd] += ∑ / dz[kwtd]
        θ[kwtd] = min(θ[kwtd], θ_sat(nsoil))

        if θ[kwtd] > θ_eq[kwtd]
          sy = θ_sat(nsoil) - θ_eq[kwtd]
          zwt = (θ[kwtd] * dz[kwtd] - θ_eq[kwtd] * z₋ₕ[iwtd] + θ_sat(nsoil) * z₋ₕ[kwtd]) / sy
          zwt = min(zwt, z₋ₕ[iwtd])
        end
        ∑ = 0.0
      else
        # Water enough to saturate the layer
        θ[kwtd] = θ_sat(nsoil)
        ∑ -= _layer
        k1 = iwtd
        for k in k1:(nzg+1)
          zwt = z₋ₕ[k]
          iwtd = k + 1
          k == nzg + 1 && break

          nsoil = soiltxt[k]
          _layer = dz[k] * (θ_sat(nsoil) - θ[k])
          if ∑ <= _layer
            θ[k] += ∑ / dz[k]
            θ[k] = min(θ[k], θ_sat(nsoil))

            # 这里比较合理的地方，每一层补给过程汇中，地下水的水位是动态上升的。
            if θ[k] > θ_eq[k]
              sy = θ_sat(nsoil) - θ_eq[k]
              zwt = (θ[k] * dz[k] - θ_eq[k] * z₋ₕ[iwtd] + θ_sat(nsoil) * z₋ₕ[k]) / sy
              zwt = min(zwt, z₋ₕ[iwtd])
            end
            ∑ = 0.0
            break
          else
            θ[k] = θ_sat(nsoil)
            ∑ -= _layer
          end
        end
      end
    elseif zwt >= z₋ₕ[1] - dz[1]
      # 这个情景是在干嘛？
      nsoil = soiltxt[1]
      _layer = (θ_sat(nsoil) - θ_wtd) * dz[1]

      if ∑ <= _layer
        θ_eq_wtd = θ_sat(nsoil) * (ψ_sat(nsoil) / (ψ_sat(nsoil) - dz[1]))^(1.0 / B(nsoil))
        θ_eq_wtd = max(θ_eq_wtd, θ_cp(nsoil))

        θ_wtd += ∑ / dz[1]
        θ_wtd = min(θ_wtd, θ_sat(nsoil))
        if θ_wtd > θ_eq_wtd
          zwt = min((θ_wtd * dz[1] - θ_eq_wtd * z₋ₕ[1] + θ_sat(nsoil) * (z₋ₕ[1] - dz[1])) / (θ_sat(nsoil) - θ_eq_wtd), z₋ₕ[1])
        end
        ∑ = 0.0
      else
        θ_wtd = θ_sat(nsoil)
        ∑ -= _layer
        for k in 1:(nzg+1)
          zwt = z₋ₕ[k]
          iwtd = k + 1
          if k == nzg + 1
            break
          end
          nsoil = soiltxt[k]
          _layer = dz[k] * (θ_sat(nsoil) - θ[k])
          if ∑ <= _layer
            θ[k] = min(θ[k] + ∑ / dz[k], θ_sat(nsoil))
            if θ[k] > θ_eq[k]
              zwt = min((θ[k] * dz[k] - θ_eq[k] * z₋ₕ[iwtd] + θ_sat(nsoil) * z₋ₕ[k]) / (θ_sat(nsoil) - θ_eq[k]), z₋ₕ[iwtd])
            end
            ∑ = 0.0
            break
          else
            θ[k] = θ_sat(nsoil)
            ∑ -= _layer
          end
        end
      end
    else
      nsoil = soiltxt[1]
      _layer = (θ_sat(nsoil) - θ_wtd) * (z₋ₕ[1] - dz[1] - zwt)
      if ∑ <= _layer
        zwt += ∑ / (θ_sat(nsoil) - θ_wtd)
        ∑ = 0.0
      else
        ∑ -= _layer
        zwt = z₋ₕ[1] - dz[1]
        _layer = (θ_sat(nsoil) - θ_wtd) * dz[1]
        if ∑ <= _layer
          θ_eq_wtd = θ_sat(nsoil) * (ψ_sat(nsoil) / (ψ_sat(nsoil) - dz[1]))^(1.0 / B(nsoil))
          θ_eq_wtd = max(θ_eq_wtd, θ_cp(nsoil))

          θ_wtd += ∑ / dz[1]
          θ_wtd = min(θ_wtd, θ_sat(nsoil))
          zwt = (θ_wtd * dz[1] - θ_eq_wtd * z₋ₕ[1] + θ_sat(nsoil) * (z₋ₕ[1] - dz[1])) / (θ_sat(nsoil) - θ_eq_wtd)
          ∑ = 0.0
        else
          θ_wtd = θ_sat(nsoil)
          ∑ -= _layer
          for k in 1:(nzg+1)
            zwt = z₋ₕ[k]
            iwtd = k + 1
            if k == nzg + 1
              break
            end
            nsoil = soiltxt[k]
            _layer = dz[k] * (θ_sat(nsoil) - θ[k])

            if ∑ <= _layer
              θ[k] += ∑ / dz[k]
              θ[k] = min(θ[k], θ_sat(nsoil))
              if θ[k] > θ_eq[k]
                zwt = (θ[k] * dz[k] - θ_eq[k] * z₋ₕ[iwtd] + θ_sat(nsoil) * z₋ₕ[k]) / (θ_sat(nsoil) - θ_eq[k])
              end
              ∑ = 0.0
              break
            else
              θ[k] = θ_sat(nsoil)
              ∑ -= _layer
            end
          end
        end
      end
    end

    # Water springing at the surface
    qspring = ∑

    # Case 2: Water table going down (totwater < 0)
  elseif ∑ < 0.0
    if zwt >= z₋ₕ[1]
      for k in 2:nzg
        if zwt < z₋ₕ[k]
          break
        end
      end
      iwtd = k

      k1 = iwtd - 1
      for kwtd in k1:-1:1
        nsoil = soiltxt[kwtd]

        # Max water that the layer can yield
        maxwatdw = dz[kwtd] * (θ[kwtd] - θ_eq[kwtd])

        if -∑ <= maxwatdw
          θ[kwtd] += ∑ / dz[kwtd]
          if θ[kwtd] > θ_eq[kwtd]
            zwt = (θ[kwtd] * dz[kwtd] - θ_eq[kwtd] * z₋ₕ[iwtd] + θ_sat(nsoil) * z₋ₕ[kwtd]) / (θ_sat(nsoil) - θ_eq[kwtd])
          else
            zwt = z₋ₕ[kwtd]
            iwtd -= 1
          end
          ∑ = 0.0
          break
        else
          zwt = z₋ₕ[kwtd]
          iwtd -= 1
          if maxwatdw >= 0.0
            θ[kwtd] = θ_eq[kwtd]
            ∑ += maxwatdw
          end
        end
      end

      if iwtd == 1 && ∑ < 0.0
        nsoil = soiltxt[1]
        θ_eq_wtd = θ_sat(nsoil) * (ψ_sat(nsoil) / (ψ_sat(nsoil) - dz[1]))^(1.0 / B(nsoil))
        θ_eq_wtd = max(θ_eq_wtd, θ_cp(nsoil))

        maxwatdw = dz[1] * (θ_wtd - θ_eq_wtd)

        if -∑ <= maxwatdw
          θ_wtd += ∑ / dz[1]
          zwt = max((θ_wtd * dz[1] - θ_eq_wtd * z₋ₕ[1] + θ_sat(nsoil) * (z₋ₕ[1] - dz[1])) / (θ_sat(nsoil) - θ_eq_wtd), z₋ₕ[1] - dz[1])
        else
          zwt = z₋ₕ[1] - dz[1]
          θ_wtd += ∑ / dz[1]
          dzup = (θ_eq_wtd - θ_wtd) * dz[1] / (θ_sat(nsoil) - θ_eq_wtd)
          zwt -= dzup
          θ_wtd = θ_eq_wtd
        end
      end
    elseif zwt >= z₋ₕ[1] - dz[1]
      nsoil = soiltxt[1]
      θ_eq_wtd = θ_sat(nsoil) * (ψ_sat(nsoil) / (ψ_sat(nsoil) - dz[1]))^(1.0 / B(nsoil))
      θ_eq_wtd = max(θ_eq_wtd, θ_cp(nsoil))

      maxwatdw = dz[1] * (θ_wtd - θ_eq_wtd)

      if -∑ <= maxwatdw
        θ_wtd += ∑ / dz[1]
        zwt = max((θ_wtd * dz[1] - θ_eq_wtd * z₋ₕ[1] + θ_sat(nsoil) * (z₋ₕ[1] - dz[1])) / (θ_sat(nsoil) - θ_eq_wtd), z₋ₕ[1] - dz[1])
      else
        zwt = z₋ₕ[1] - dz[1]
        θ_wtd += ∑ / dz[1]
        dzup = (θ_eq_wtd - θ_wtd) * dz[1] / (θ_sat(nsoil) - θ_eq_wtd)
        zwt -= dzup
        θ_wtd = θ_eq_wtd
      end
    else
      nsoil = soiltxt[1]
      wgpmid = θ_sat(nsoil) * (ψ_sat(nsoil) / (ψ_sat(nsoil) - (z₋ₕ[1] - zwt)))^(1.0 / B(nsoil))
      wgpmid = max(wgpmid, θ_cp(nsoil))
      s_yield_dw = θ_sat(nsoil) - wgpmid
      wtdold = zwt
      zwt = wtdold + ∑ / s_yield_dw
      θ_wtd = (θ_wtd * (z₋ₕ[1] - wtdold) + wgpmid * (wtdold - zwt)) / (z₋ₕ[1] - zwt)
    end
    qspring = 0.0
  end

  return zwt, qspring, θ_wtd
end
