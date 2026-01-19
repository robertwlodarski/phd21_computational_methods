# Standard macros file used across projects
# Ignore if you don't see that it's particularly used. 

# Content 
# 1. macro insert!
# 2. displayall
# 3. ⪅ / ⪆

# 1. In-place version of @insert from Accessors.jl
macro insert!(ex)
    Accessors.insertmacro(identity, ex, overwrite=true)
end

# 2. Display setting
function displayall(x)
    show(stdout, "text/plain", x)
end

# 3. Inequality signs 
⪅(x, y) = x < y || x ≈ y
⪆(x, y) = x > y || x ≈ y

# 4. Discrete distributions sampling from StatsBase.jl
# Note: This is required for CUDA, which means it may not need to be used. 
function sample_gpu(rng, wv::AbstractVector)
    1 == firstindex(wv) ||
        throw(ArgumentError("non 1-based arrays are not supported"))
    t = rand(rng) * sum(wv)
    n = length(wv)
    i = 1
    cw = wv[1]
    while cw < t && i < n
        i += 1
        @inbounds cw += wv[i]
    end
    return i
end
sample_gpu(rng, a::AbstractArray, wv::AbstractVector) = a[sample_gpu(rng, wv)]