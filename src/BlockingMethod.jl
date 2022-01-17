module BlockingMethod
# Author: David James, davidabraham@ucla.edu
# Reference: Flyvbjerg and Petersen, J. Chem. Phys., 91, 461, 1989

export estimate

"""
    estimate(x::AbstractArray)

Performs time series analysis on the `x` column of data outputting
the mean and standard deviation.

Reference Flyvbjerb and Petersen '89 for more info.
"""
function estimate(x::AbstractArray)

    n = length(x)

    xm = sum(x) / n
    x2m = sum(y -> y^2, x) / n
    σ = 0

    # if too small exit out
    if (n < 2) @goto fin end

    # number of bins
    nrbin = log(n) / log(2) - 1

    ci = (x2m - xm^2) / (n - 1)
    faci = 1 / sqrt(2*(n - 1))

    dc = faci * ci

    c_max = ci
    c_old = 0

    # calculating sigma
    for i in 1:nrbin

        c = 0
        n = floor(Int64, n / 2)

        for j in 1:n
            if n != 1
                x[j] = (x[2*j] + x[2*j-1]) / 2
                c += (x[j] - xm)^2 / (n * (n - 1))
            end
        end

        fac = 1 / sqrt(2*(n - 1))
        dc = fac*sqrt(c)
        diff = sqrt(c) - sqrt(c_old)

        if (abs(diff) < dc && c > c_max) c_max = c end

        c_old = c

    end

    σ = sqrt(c_max)

    @label fin
    return (xm, σ)

end

end # module
