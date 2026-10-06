module BlockingMethod
# Author: David James, djames@epss.ucla.edu
# Reference: Flyvbjerg and Petersen, J. Chem. Phys., 91, 461, 1989

export block_mean

Base.@deprecate estimate(x) block_mean(x)

"""
    block_mean(x) -> (mean, σ)

Mean of a correlated series with its Flyvbjerg–Petersen blocking error.
Does not modify `x`.
"""
block_mean(x::AbstractVector{<:Real}) = block_mean!(collect(Float64, x))

"""
    block_mean(x::AbstractVector{<:Float64})::NTuple{2, Float64}

Same as `block_mean`, but uses `x` as workspace. On return, `x` holds block
averages instead of the original series. Use it only when you don't need `x`
afterwards and want to skip the copy.

# Flyvbjerb and Petersen 1989 - Abstract
We describe how the true statistical error on an average of correlated data
can be obtained with ease and efficiency by renormalization group method [...]
Reference article https://doi.org/10.1063/1.457480 for more info.
"""
function block_mean!(x::AbstractVector{<:Float64})::NTuple{2, Float64}

    n = length(x)
    fn = Float64(n)

    xm = sum(x) / fn
    x2m = sum(y -> y^2, x) / fn
    σ = 0.0

    # if too small exit out
    if n < 2
        # Handle edge case where there are too few points.
        return (xm, σ)
    end

    # number of bins
    nrbin = log(fn) / log(2) - 1.0

    c_max = (x2m - xm^2) / (fn - 1.0)
    fac = 1. / sqrt(2. * fn - 2.0)

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
                fac = 1. / sqrt(2.0 * fn - 2.0)
                x[j] = (x[2*j] + x[2*j-1]) / 2.0
                c += (x[j] - xm)^2 / (fn^2 - fn)
            end
        end

        dc = fac*sqrt(c)
        diff = sqrt(c) - sqrt(c_old)

        if (abs(diff) < dc && c > c_max) c_max = c end

        c_old = c

    end

    σ = sqrt(c_max)

    return (xm, σ)
end

end # module
