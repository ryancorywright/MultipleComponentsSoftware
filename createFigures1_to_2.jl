include("core_julia1p7.jl")
include("synthetic_data_generation.jl")

include("sdp_bounds.jl")
include("iterative_deflation3.jl")
include("covariance_thresholding.jl")


include("custombb/branchAndBound.jl")
include("custombb/multiComponents.jl")
include("custombb/utilities.jl")

#### Support Recovery Metrics ####
support_overlap(a,b) = length(intersect(a,b))

# One spike
function support_recovery_metrics(Xhat, x1)
    support_truth_1 = findall(abs.(x1) .> 1e-8)

    support_recovered_1 = findall(abs.(Xhat)[:] .> 1e-8)

    A = 1.0*support_overlap(support_recovered_1,support_truth_1) / (length(support_truth_1))
    FDR = 1.0*length(setdiff(support_recovered_1,support_truth_1)) / (length(support_recovered_1))

    return A, FDR, A
end
# Two spikes
function support_recovery_metrics(Xhat, x1, x2)
    support_truth_1 = findall(abs.(x1) .> 1e-8)
    support_truth_2 = findall(abs.(x2) .> 1e-8)
    union_support_truth = unique(union(support_truth_1, support_truth_2))

    support_recovered_1 = findall(abs.(Xhat[:,1]) .> 1e-8)
    support_recovered_2 = findall(abs.(Xhat[:,2]) .> 1e-8)
    union_support_recovered = unique(union(support_recovered_1, support_recovered_2))

    option_1 = support_overlap(support_recovered_1,support_truth_1) + support_overlap(support_recovered_2,support_truth_2) 
    option_2 = support_overlap(support_recovered_1,support_truth_2) + support_overlap(support_recovered_2,support_truth_1) 

    A_perPC = 1.0*max(option_1, option_2) / (length(support_truth_1) + length(support_truth_2))
    
    A = 1.0*support_overlap(union_support_recovered,union_support_truth) / (length(union_support_truth))
    FDR = 1.0*length(setdiff(union_support_recovered,union_support_truth)) / (length(union_support_recovered))

    # @show (length(support_truth_1) + length(support_truth_2)), (length(union_support_truth)), (length(union_support_recovered))
    return A, FDR, A_perPC
end
#### #### #### ####


mkpath("synthetic_data_experiments/results/")

using DataFrames 
results = DataFrame(nexp=[], snr=[], overlap=[], nsamples=[], method=[], time=[], A_union=[], FDR_union=[], A_perPC=[], var_explained=[], inner_prod=[])


## Parameters of the experiments
snr = 2.0
r = 2

array_num = parse(Int, ARGS[1]) #Input parameter: Controls overlap (q) + id of the experiments (nexp)

q = [0.1, 0.5, 0.90][(array_num % 3)+1]
nexp = div(array_num, 3) + 1


p = 50
k_true = 20 
fixing = true #Decide whether to run the methods with/without variable fixing

using DelimitedFiles, LinearAlgebra

#Generate true covariance model
# Σ, x1, x2 = generate_2spike_cov(p, k_true, q, snr, seed=Int(1532+q*20), binarize=true)
Σ = DelimitedFiles.readdlm("synthetic_data_experiments/data/covariance_matrix_p_$(p)_k_$(k_true)_q_$(q).csv", ',')
X = DelimitedFiles.readdlm("synthetic_data_experiments/data/spikes_p_$(p)_k_$(k_true)_q_$(q).csv",  ',')
x1 = X[:,1]; x2 = X[:,2]

# #Generate observations
# using Distributions
# d = MvNormal(zeros(p), Σ)
# Random.seed!(Int(1532+q*20))
# X = rand(d, 10000) #p by N matrix of observations

# for nexp in 1:20
    nseed = 1532+nexp*78
    Random.seed!(nseed)
    # X = X[:,shuffle(1:10000)]

    nrange = union((10:30:310),(350:50:1200))
    for n = nrange
        @show q, nexp, n
        Sigma_hat = DelimitedFiles.readdlm("synthetic_data_experiments/data/covariance_matrix_p_$(p)_k_$(k_true)_q_$(q)_n_$(n)_nexp_$(nexp).csv", ',')
        D = inv(sqrt(Diagonal(Sigma_hat)))
        Cor_hat = D*Sigma_hat*D
        # Sigma_hat = cov(X[:,1:n]') #Empirical covariance matrix
        # Cor_hat = cor(X[:,1:n]') #Empirical correlation matrix

        indices_reduced = collect(1:p)
        
        if fixing 
            ### Variable Fixing: Solve PSD relaxation and then round it, to obtain the feasible solution
            println("-----Solving relaxation for variable fixing-----")
            run_relax = @elapsed _, _, _, _, z_relax = getSDPUpperBound_gd_multiple_permutation(Sigma_hat, r, repeat([k_true], r),  
                        usePSDs=true, useSOCs=true, useCuts=false, addDisjointConstraint=false, generateDisjointFeasible=false, verbose=false)

            A, FDR, A_perPC = support_recovery_metrics(z_relax, x1, x2)
            
            push!(results, [nexp, snr, q, n, "Relaxation",  run_relax, A, FDR, A_perPC, 0., Inf])

            indices_reduced = findall(z_relax*ones(r) .>= 1e-3)
            @show length(indices_reduced)
            # indices_reduced = 1:p
            Sigma_reduced = Sigma_hat[indices_reduced,indices_reduced]
            Cor_reduced = Cor_hat[indices_reduced,indices_reduced]


            ### Algorithm 1: Solve PSD relaxation and then round it, to obtain the feasible solution
            println("-----Algorithm 1-----")
            #Covariance matrix

            #OLD Method
            run_sdp_disj = @elapsed UB_sdp_disj, x_rounded, LB_gd_disj, violation_gd_disj, z_relax = getSDPUpperBound_gd_multiple_permutation(Sigma_reduced, r, repeat([k_true], r),  
                    usePSDs = true, useSOCs=false, useCuts=false, addDisjointConstraint=false, generateDisjointFeasible=true, verbose = false,
                    breakSymmetry = false, UBpenalty_rounding = 0.0)

            U = zeros(p,r); U[indices_reduced,:] = x_rounded
            A, FDR, A_perPC = support_recovery_metrics(U, x1, x2)
            var_explained = tr(U'*Σ*U)/tr(Σ)
            push!(results, [nexp, snr, q, n, "Algorithm 1 - OLD (fixed) - COV",  run_sdp_disj, A, FDR, A_perPC, var_explained, dot(U[:,1],U[:,2])])


            #New Rounding
            run_sdp_disj = @elapsed UB_sdp_disj, x_rounded, LB_gd_disj, violation_gd_disj, z_relax = getSDPUpperBound_gd_multiple_permutation(Sigma_reduced, r, repeat([k_true], r),  
                    usePSDs = true, useSOCs=false, useCuts=false, addDisjointConstraint=false, generateDisjointFeasible=true, verbose = false,
                    breakSymmetry = false, UBpenalty_rounding = 1.0)

            U = zeros(p,r); U[indices_reduced,:] = x_rounded
            A, FDR, A_perPC = support_recovery_metrics(U, x1, x2)
            var_explained = tr(U'*Σ*U)/tr(Σ)
            push!(results, [nexp, snr, q, n, "Algorithm 1 - New rounding (fixed) - COV",  run_sdp_disj, A, FDR, A_perPC, var_explained, dot(U[:,1],U[:,2])])

            # #New Rounding + disjoint constraint
            # run_sdp_disj = @elapsed UB_sdp_disj, x_rounded, LB_gd_disj, violation_gd_disj, z_relax = getSDPUpperBound_gd_multiple_permutation(Sigma_reduced, r, repeat([k_true], r),  
            #         usePSDs = true, useSOCs=false, useCuts=false, addDisjointConstraint=true, generateDisjointFeasible=true, verbose = false,
            #         breakSymmetry = false, UBpenalty_rounding = 1.0)

            # U = zeros(p,r); U[indices_reduced,:] = x_rounded
            # A, FDR, A_perPC = support_recovery_metrics(U, x1, x2)
            # var_explained = tr(U'*Σ*U)/tr(Σ)
            # push!(results, [nexp, snr, q, n, "Algorithm 1 - New rounding + Disj constraint (fixed) - COV",  run_sdp_disj, A, FDR, A_perPC, var_explained, dot(U[:,1],U[:,2])])


            #Correlation matrix

            #OLD Method
            run_sdp_disj = @elapsed UB_sdp_disj, x_rounded, LB_gd_disj, violation_gd_disj, z_relax = getSDPUpperBound_gd_multiple_permutation(Cor_reduced, r, repeat([k_true], r),  
                    usePSDs = true, useSOCs=false, useCuts=false, addDisjointConstraint=false, generateDisjointFeasible=true, verbose = false,
                    breakSymmetry = false, UBpenalty_rounding = 0.0)

            U = zeros(p,r); U[indices_reduced,:] = x_rounded
            A, FDR, A_perPC = support_recovery_metrics(U, x1, x2)
            var_explained = tr(U'*Σ*U)/tr(Σ)
            push!(results, [nexp, snr, q, n, "Algorithm 1 - OLD (fixed) - COR",  run_sdp_disj, A, FDR, A_perPC, var_explained, dot(U[:,1],U[:,2])])


            #New Rounding
            run_sdp_disj = @elapsed UB_sdp_disj, x_rounded, LB_gd_disj, violation_gd_disj, z_relax = getSDPUpperBound_gd_multiple_permutation(Cor_reduced, r, repeat([k_true], r),  
                    usePSDs = true, useSOCs=false, useCuts=false, addDisjointConstraint=false, generateDisjointFeasible=true, verbose = false,
                    breakSymmetry = false, UBpenalty_rounding = 1.0)

            U = zeros(p,r); U[indices_reduced,:] = x_rounded
            A, FDR, A_perPC = support_recovery_metrics(U, x1, x2)
            var_explained = tr(U'*Σ*U)/tr(Σ)
            push!(results, [nexp, snr, q, n, "Algorithm 1 - New rounding (fixed) - COR",  run_sdp_disj, A, FDR, A_perPC, var_explained, dot(U[:,1],U[:,2])])


            # #New Rounding + disjoint constraint
            # run_sdp_disj = @elapsed UB_sdp_disj, x_rounded, LB_gd_disj, violation_gd_disj, z_relax = getSDPUpperBound_gd_multiple_permutation(Cor_reduced, r, repeat([k_true], r),  
            #         usePSDs = true, useSOCs=false, useCuts=false, addDisjointConstraint=true, generateDisjointFeasible=true, verbose = false,
            #         breakSymmetry = false, UBpenalty_rounding = 1.0)

            # U = zeros(p,r); U[indices_reduced,:] = x_rounded
            # A, FDR, A_perPC = support_recovery_metrics(U, x1, x2)
            # var_explained = tr(U'*Σ*U)/tr(Σ)
            # push!(results, [nexp, snr, q, n, "Algorithm 1 - New rounding + Disj constraint (fixed) - COR",  run_sdp_disj, A, FDR, A_perPC, var_explained, dot(U[:,1],U[:,2])])

            
            ### Algorithm 2: Iterative deflation
            println("-----Algorithm 2-----")
            numIters = 200
            @time _, _, runtime_AM, x_AM_reduced = findmultPCs_deflation(Sigma_reduced, r, repeat([k_true], r), numIters=numIters, verbose=false)

            U = zeros(p,r); U[indices_reduced,:] = x_AM_reduced
            A, FDR, A_perPC = support_recovery_metrics(U, x1, x2)
            var_explained = tr(U'*Σ*U)/tr(Σ)

            push!(results, [nexp, snr, q, n, "Algorithm 2 (fixed) - COV",  runtime_AM, A, FDR, A_perPC, var_explained, dot(U[:,1],U[:,2])])

        
            @time _, _, runtime_AM, x_AM_reduced = findmultPCs_deflation(Cor_reduced, r, repeat([k_true], r), numIters=numIters, verbose=false)

            U = zeros(p,r); U[indices_reduced,:] = x_AM_reduced
            A, FDR, A_perPC = support_recovery_metrics(U, x1, x2)
            var_explained = tr(U'*Σ*U)/tr(Σ)

            push!(results, [nexp, snr, q, n, "Algorithm 2 (fixed) - COR",  runtime_AM, A, FDR, A_perPC, var_explained, dot(U[:,1],U[:,2])])
        

            ### Bertsimas and Berk (2019)
            println("-----Berk and Bertsimas----")
            ## With variable fixing -- not huge change
            B = sqrt(Sigma_reduced)
            theProb = problem(B, Sigma_reduced)
            theU, theRuntime = multiOptimalSPCA(theProb, k_true, r)

            U = zeros(p,r); U[indices_reduced,:] = theU'
            A, FDR, A_perPC = support_recovery_metrics(U, x1, x2)
            var_explained = tr(U'*Σ*U)/tr(Σ)
            
            push!(results, [nexp, snr, q, n, "Berk and Bertsimas (2019) (fixed) - COV",  theRuntime, A, FDR, A_perPC, var_explained, dot(U[:,1],U[:,2])])

            try
                B = sqrt(Cor_reduced)
                theProb = problem(B, Cor_reduced)
                theU, theRuntime = multiOptimalSPCA(theProb, k_true, r)

                U = zeros(p,r); U[indices_reduced,:] = theU'
                A, FDR, A_perPC = support_recovery_metrics(U, x1, x2)
                var_explained = tr(U'*Σ*U)/tr(Σ)
                
                push!(results, [nexp, snr, q, n, "Berk and Bertsimas (2019) (fixed) - COR",  theRuntime, A, FDR, A_perPC, var_explained, dot(U[:,1],U[:,2])])
            catch
                push!(results, [nexp, snr, q, n, "Berk and Bertsimas (2019) (fixed) - COR",  NaN, 0., 1., 0., 0., 0.])
            end
        end

        if !fixing
            ### Algorithm 2: Iterative deflation
            println("-----Algorithm 2-----")
            numIters = 200
            @time _, _, runtime_AM, x_AM_reduced = findmultPCs_deflation(Sigma_hat, r, repeat([k_true], r), numIters=numIters, verbose=false)

            U = x_AM_reduced
            A, FDR, A_perPC = support_recovery_metrics(U, x1, x2)
            var_explained = tr(U'*Σ*U)/tr(Σ)

            push!(results, [nexp, snr, q, n, "Algorithm 2 - COV",  runtime_AM, A, FDR, A_perPC, var_explained, dot(U[:,1],U[:,2])])

        
            @time _, _, runtime_AM, x_AM_reduced = findmultPCs_deflation(Cor_hat, r, repeat([k_true], r), numIters=numIters, verbose=false)

            U = x_AM_reduced
            A, FDR, A_perPC = support_recovery_metrics(U, x1, x2)
            var_explained = tr(U'*Σ*U)/tr(Σ)

            push!(results, [nexp, snr, q, n, "Algorithm 2 - COR",  runtime_AM, A, FDR, A_perPC, var_explained, dot(U[:,1],U[:,2])])


            ### Covariance Thresholding
            println("-----Covariance Thresholding----")
            run_covthresh = @elapsed V = covarianceThresholding_adaptive(Sigma_hat, r, repeat([k_true], r)) 
            A, FDR, A_perPC = support_recovery_metrics(V, x1, x2)
            var_explained = tr(V'*Σ*V)/tr(Σ) 
            
            push!(results, [nexp, snr, q, n, "Covariance Thresholding - COV",  run_covthresh, A, FDR, A_perPC, var_explained, dot(V[:,1],V[:,2])])

            run_covthresh = @elapsed V = covarianceThresholding_adaptive(Cor_hat, r, repeat([k_true], r)) 
            A, FDR, A_perPC = support_recovery_metrics(V, x1, x2)
            var_explained = tr(V'*Σ*V)/tr(Σ) 
            
            push!(results, [nexp, snr, q, n, "Covariance Thresholding - COR",  run_covthresh, A, FDR, A_perPC, var_explained, dot(V[:,1],V[:,2])])


            ### Bertsimas and Berk (2019)
            println("-----Berk and Bertsimas----")
            B = sqrt(Sigma_hat)
            theProb = problem(B, Sigma_hat)
            theU, theRuntime = multiOptimalSPCA(theProb, k_true, r)

            U = theU'
            A, FDR, A_perPC = support_recovery_metrics(U, x1, x2)
            var_explained = tr(U'*Σ*U)/tr(Σ)
            
            push!(results, [nexp, snr, q, n, "Berk and Bertsimas (2019) - COV",  theRuntime, A, FDR, A_perPC, var_explained, dot(U[:,1],U[:,2])])
            
            try
                B = sqrt(Cor_hat)
                theProb = problem(B, Cor_hat)
                theU, theRuntime = multiOptimalSPCA(theProb, k_true, r)

                U = theU'
                A, FDR, A_perPC = support_recovery_metrics(U, x1, x2)
                var_explained = tr(U'*Σ*U)/tr(Σ)
                
                push!(results, [nexp, snr, q, n, "Berk and Bertsimas (2019) - COR",  theRuntime, A, FDR, A_perPC, var_explained, dot(U[:,1],U[:,2])])
            catch
                push!(results, [nexp, snr, q, n, "Berk and Bertsimas (2019) - COR",  NaN, 0., 1., 0., 0., 0.])
            end
        end

        # CSV.write("synthetic_data_experiments/results/results_2spike_p_$(p)_k_$(k_true)"*(fixing ? "_fixing" : "_notfixing")*"_q_$(q)_nexp_$(nexp).csv", results)
    end
# end