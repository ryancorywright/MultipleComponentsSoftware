include("core_julia1p7.jl")
include("synthetic_data_generation.jl")

using DelimitedFiles

mkpath("synthetic_data_experiments/data/")


## Parameters of the experiments
p = 50
k_true = 20
snr = 1.0

r = 2
for (p,k_true) in [(50,20)] 
    for q in [0.1, 0.5, 0.90] 
    # for q in [0.1] 
        #Generate true covariance model
        Σ, x1, x2 = generate_2spike_cov(p, k_true, q, snr, seed=Int(1532+q*20), binarize=true)

        DelimitedFiles.writedlm("synthetic_data_experiments/data/spikes_p_$(p)_k_$(k_true)_q_$(q).csv", [x1 x2], ',')
        DelimitedFiles.writedlm("synthetic_data_experiments/data/covariance_matrix_p_$(p)_k_$(k_true)_q_$(q).csv", Σ, ',')


        #Generate observations
        using Distributions
        d = MvNormal(zeros(p), Σ)
        Random.seed!(Int(1532+q*20))
        X = rand(d, 10000) #p by N matrix of observations

        for nexp in 1:20
            nseed = 1532+nexp*78
            Random.seed!(nseed)
            X = X[:,shuffle(1:10000)]

            nrange = union((10:30:310), (350:50:1200))
            for n = nrange
                @show q, nexp, n
                Sigma_hat = cov(X[:,1:n]') #Empirical covariance matrix

                DelimitedFiles.writedlm("synthetic_data_experiments/data/covariance_matrix_p_$(p)_k_$(k_true)_q_$(q)_n_$(n)_nexp_$(nexp).csv", Sigma_hat, ',')
            end
        end 
    end
end