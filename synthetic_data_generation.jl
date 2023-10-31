using Random, LinearAlgebra

function generate_1spike_cov(p, k_true, β; seed=8735, binarize::Bool=false)
        Random.seed!(seed)
        x1 = zeros(p)
        x1[1:k_true] = rand(k_true)
        if binarize
            x1[1:k_true] .= sign.(x1[1:k_true] .- .5)
        end

        shufflecoords = randperm(p)
        x1 = x1[shufflecoords]        
        x1 /= norm(x1)
        @assert sum(abs.(x1) .> 0) == k_true

        S = β*x1*x1'+ Matrix(1.0*I, p, p)
        S = (S + S')/2
        return S, x1 
    end

    
function generate_2spike_cov(p, k_true, prop_overlap, β; seed=8735, binarize::Bool=false)
    Random.seed!(seed)

    k_overlap = floor(Int, prop_overlap*k_true)
    k_overlap += 1*(k_overlap % 2 == 1) # Make sure k_overlap is even
    if k_overlap == k_true
        k_overlap -= 2
    end
    k_overlap_half =Int(k_overlap/2)
    k_nonoverlap = k_true - k_overlap

    x1 = zeros(p); x2=zeros(p)
    if binarize
        #Non-overlapping support: no restrictions
        x1[1:k_nonoverlap] = sign.(rand(k_nonoverlap) .- .5)
        x2[(k_nonoverlap+1):(k_nonoverlap+k_nonoverlap)] = sign.(rand(k_nonoverlap) .- .5)
        #Overlapping support:
        x1[(2*k_nonoverlap + 1):(2*k_nonoverlap + k_overlap)] = sign.( rand(k_overlap) .-.5 )
        #x1 and x2 agree on the first half 
        x2[(2*k_nonoverlap + 1):(2*k_nonoverlap + k_overlap_half)] .= x1[(2*k_nonoverlap + 1):(2*k_nonoverlap + k_overlap_half)]
        #x1 and x2 disagree on the secon half 
        x2[(2*k_nonoverlap + k_overlap_half + 1):(2*k_nonoverlap + k_overlap)] .= -x1[(2*k_nonoverlap + k_overlap_half + 1):(2*k_nonoverlap + k_overlap)]
    else 
        #Non-overlapping support: no restrictions
        x1[1:k_nonoverlap] = (rand(k_nonoverlap) .- .5)
        x2[(k_nonoverlap+1):(k_nonoverlap+k_nonoverlap)] = (rand(k_nonoverlap) .- .5)
        #Overlapping support:
        ov_index = (2*k_nonoverlap + 1):(2*k_nonoverlap + k_overlap)
        x1[ov_index] = ( rand(k_overlap) .-.5 )
        #x1 and x2 agree on the first half 
        x2[ov_index] = ( rand(k_overlap) .-.5 )
        sc = dot(x1[ov_index], x2[ov_index])/dot(x1[ov_index],x1[ov_index])
        x2[ov_index] .-= sc*x1[ov_index]
        x2[ov_index] ./= maximum(abs.(x2[ov_index]))
    end    
    shufflecoords = randperm(p)
    x1 = x1[shufflecoords]; x2=x2[shufflecoords] 
    
    @assert sum(abs.(x1) .> 0) == k_true
    @assert sum(abs.(x2) .> 0) == k_true
    @assert abs(dot(x1,x2)) ≤ 1e-10
    
    x1 /= norm(x1); x2 /= norm(x2) 

    S = β*x1*x1'+β*x2*x2'+ Matrix(1.0*I, p, p)
    S = (S + S')/2

    return S, x1, x2 
end

# support_overlap(a,b) = length(intersect(a,b))