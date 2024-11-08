module BlockingMethod
# Author: David James, davidabraham@ucla.edu
# Reference: Flyvbjerg and Petersen, J. Chem. Phys., 91, 461, 1989

export estimate

"""
    estimate(x::Vector{Float64})::NTuple{2, Float64}

Performs time series analysis on the `x` column of data outputting
the mean and standard error (variance).

# Flyvbjerb and Petersen 1989 - Abstract
We describe how the true statistical error on an average of correlated data
can be obtained with ease and efficiency by renormalization group method [...]
Reference article https://doi.org/10.1063/1.457480 for more info.
"""
function estimate(x::Vector{Float64})::NTuple{2, Float64}

    n = length(x)
    fn = Float64(n)

    xm = sum(x) / fn
    x2m = sum(y -> y^2, x) / fn
    σ = 0.0

    # if too small exit out
    if (n < 2) @goto fin end

    # number of bins
    nrbin = log(fn) / log(2) - 1.

    c_max = (x2m - xm^2) / (fn - 1.)
    fac = 1. / sqrt(2. * fn - 2.)

    dc = fac * c_max

    c_old = 0.0

    # calculating sigma
    for i in 1:nrbin

        c = 0.0
        n = floor(Int64, n / 2)
        fn = Float64(n)

        for j in 1:n
            if(n == 1)
                continue
            else
                fac = 1. / sqrt(2. * fn - 2.)
                x[j] = (x[2*j] + x[2*j-1]) / 2.
                c += (x[j] - xm)^2 / (fn^2 - fn)
            end
        end

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
