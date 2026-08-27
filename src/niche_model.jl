function niche(S::Int, C::AbstractFloat, rng=Random.default_rng())
    # variables definition
    adj_matrix   =  zeros(Int64, S, S)  # all species start having no interactions
    in_degree    = zeros(S)
    out_degree   = zeros(S)
    disconnected = 1 : S  # no interactions at the beginning
    n_disconnected = S > 1 ? S : 0  # if S = 1 there are no nodes to connect

    n_i = zeros(S)  # 'niche value' parameter
    r_i = zeros(S)  # feeding range
    c_i = zeros(S)  # range center

    #= 
    Species i consumes all species falling in a range `r_i`` that is placed by uniformly drawing the 
    centre of the range `c_i` from [r_i/2, n_i]. The size of r_i is assigned by using a beta function to randomly 
    draw values from [0, 1] whose expected value is 2C and then multiplying that value by `n_i`
    to obtain the desired C.
    =#
    β      = (1 / (2 * C)) - 1
    β_dist = Beta(1, β)

    while n_disconnected > 0
        connect_nodes!(adj_matrix, n_i, r_i, c_i, disconnected, β_dist, rng)

        in_degree[:]   = sum(adj_matrix, dims=1)  # preys
        out_degree[:]  = sum(adj_matrix, dims=2)  # predators
        disconnected   = findall((in_degree .== 0) .& (out_degree .== 0))  # no preys no predators
        n_disconnected = length(disconnected)
    end
    return adj_matrix
end

function connect_nodes!(adj_matrix, n_i, r_i, c_i, disconnected, β_dist, rng=Random.default_rng())
    n_disconnected = length(disconnected)

    tmp_n_i            = zeros(n_disconnected)
    fast_rand!(rng, tmp_n_i, Uniform(0, 1))
    sort!(tmp_n_i)
    n_i[disconnected] .= tmp_n_i
    tmp_r_i            = zeros(n_disconnected)
    fast_rand!(rng, tmp_r_i, β_dist)
    r_i[disconnected] .= tmp_r_i .* tmp_n_i

    r_i[argmin(n_i)] = 0  # The species with the smallest n_i has r_i = 0 so that every web has at least one basal species

    min_dist = r_i[disconnected] ./ 2  # minimum distance from the center c_i
    max_dist = n_i[disconnected]       # maximum distance from the center c_i
    c_i[disconnected] .= rand.(rng, Uniform.(min_dist, max_dist))
    
    min_x_i = c_i .- (r_i ./ 2)
    max_x_i = c_i .+ (r_i ./ 2)
    
    for i in disconnected  # axes(adj_matrix, 2)
        @. adj_matrix[i, :] = !isless(n_i, min_x_i[i]) * isless(n_i, max_x_i[i])  # i eats j (row -> column)
    end
    
    for i in disconnected  # axes(adj_matrix, 1)
        adj_matrix[i, i] = 0
    end
    nothing
end

"""
    fast_rand!([rng=default_rng()], A::AbstractArray, distribution::Sampleable)

Samples in-place from the sampler and stores the result in the provided array. Custom function. Faster than `Distributions.rand!`.
"""
function fast_rand!(rng, A::AbstractArray, distribution::Sampleable)
    for i in eachindex(A)
        @inbounds A[i] = rand(rng, distribution)
    end
    nothing
end

function fast_rand!(A::AbstractArray, distribution::Sampleable)
    for i in eachindex(A)
        @inbounds A[i] = rand(default_rng(), distribution)
    end
    nothing
end

# Populates array A with size(A) random elements. Each element is sampled from the corresponding element of an array of distributions.
"""
    fast_rand!([rng=default_rng()], A::AbstractArray, distributions::Array{T} where {T<:Sampleable})

Samples in-place from the array of samplers and stores the result in the provided array. Custom function.
"""
function fast_rand!(rng, A::AbstractArray, distributions::AbstractArray{T} where {T<:Sampleable})
    for i in eachindex(A)
        @inbounds A[i] = rand(rng, distributions[i])
    end
    nothing
end

function fast_rand!(A::AbstractArray, distributions::AbstractArray{T} where {T<:Sampleable})
    for i in eachindex(A)
        @inbounds A[i] = rand(default_rng(), distributions[i])
    end
    nothing
end