# Translated from interception.f90

function interception(minpprate, precip, lai, intercepstore, ppdrip, pet_i, et_i)
    intercepmax = 0.2 * lai
    deficit = intercepmax - intercepstore

    if precip > deficit
        if precip < minpprate
            et_i = min(intercepmax, pet_i)
        else
            et_i = 0.0
        end
        intercepstore = intercepmax - et_i
        ppdrip = precip - deficit
    else
        if precip < minpprate
            et_i = min(intercepstore + precip, pet_i)
        else
            et_i = 0.0
        end
        intercepstore = intercepstore + precip - et_i
        ppdrip = 0.0
    end

    return intercepstore, ppdrip, et_i
end
