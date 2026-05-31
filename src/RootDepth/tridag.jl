# Translated from tridag.f90

function tridag(a, b, c, r, u, n)
    gam = zeros(eltype(u), n)
    b[1] == 0.0 && error("tridag: rewrite equations")

    bet = b[1]
    u[1] = r[1] / bet
    for j in 2:n
        gam[j] = c[j - 1] / bet
        bet = b[j] - a[j] * gam[j]
        if bet == 0.0
            println("tridag failed at j=", j, " b=", b[j], " a=", a[j], " gam=", gam[j])
            error("tridag failed")
        end
        u[j] = (r[j] - a[j] * u[j - 1]) / bet
    end

    for j in (n - 1):-1:1
        u[j] = u[j] - gam[j + 1] * u[j + 1]
    end
    return nothing
end
