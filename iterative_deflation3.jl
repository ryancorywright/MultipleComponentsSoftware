include("sdp_bounds.jl")
include("covariance_thresholding.jl")

#results_template2 = DataFrame(Iter=Int[], ofv=Real[], viol=Real[])
using Printf
function findmultPCs_deflation(Sigma::Array{Float64, 2}, r::Int64, ks::Array{Int64,1};
                    numIters::Int64 = 30, 
                    # violationPenalty::Float64 = -10.0, 
                    innerSolveMethod::Symbol=:BB_warmstart,
                    x_init::Array{Float64, 2} = zeros(size(Sigma,1), r),
                    verbose::Bool = true,
                    violation_tolerance::Real = 1e-4 # Pick the best solution within the solution path that satisfies this violation tolerance, or the final one if none do
                    )
    n = size(Sigma, 1)

    ofv_best = -1e10
    violation_best = n
    x_best = zeros(n,r)
    x_current = x_init[:,:] # Trying this to see if it improves our results, which are sometimes not so great when we start at 0
    violation_init = max(sum(abs.(x_init'*x_init .- Diagonal(ones(r)))),1e-7)
    ofv_prev = 0.0
    ofv_overall = 0.0
    
    # Step size tuning
    weights = zeros(r)

    theLambda = sum(weights) > 0. ? 0.1*violation_init : 0. #Vector of penalty parameters
    @show weights, theLambda
    stepSize = 0.0
    slowPeriod = ceil(Int, 0.15*numIters)
    fastPeriod = ceil(Int, 0.75*numIters)
    # Parameters for solving the relaxation of the rank-1 SPCA problem 
    useSOC= (n > 50)
    usePSD= (n <= 50)

    if verbose 
        println("---- Iterative deflation algorithm for sparse PCA with multiple PCs ---")
        println("Dimension: $(n)")
        println("Number of PCs: $(r)")
        println("Sparsity pattern: $(ks)")
        println()
        Printf.@printf(" %10s | %20s | %25s | %10s ", "Iteration", "Objective value", "Orthogonality Violation", "Time")
        println()
    end

    start_time = time()
    for theIter = 1:numIters
        theLambda += stepSize

        for t = (1:r)
            sigma_current = Sigma-theLambda*sum(weights[s]*x_current[:,s]*x_current[:,s]' for s=1:r if s != t);
            sigma_current = (sigma_current + sigma_current')/2
            λ0 = -eigmin(sigma_current)+1e-4
            
            sigma_current += λ0*I
            sigma_current = real.(sigma_current)
            
            lambda_partial = 0; x_output=zeros(n)
            if innerSolveMethod == :BB_warmstart
                lambda_partial, x_output = subset(problem(sqrt(sigma_current), sigma_current), ks[t], timeLimit = 20)
            elseif innerSolveMethod == :SOC_relaxation
                _, lambda_partial, x_output, = getSDPUpperBound_gd_permutation(sigma_current, ks[t],
                            useSOCS=true, usePSDs=false, verbose= false)
            elseif innerSolveMethod == :SDP_relaxation
                _, lambda_partial, x_output, = getSDPUpperBound_gd_permutation(sigma_current, ks[t],
                            useSOCS=false, usePSDs=true, verbose= false)
            elseif innerSolveMethod == :Cov_Thresholding
                x_output = covarianceThresholding_adaptive(sigma_current, 1, [ks[t]])
                x_output = vec(x_output)
                lambda_partial = tr(x_output'*sigma_current*x_output)
            else 
                error("Unknown innerSolveMethod (method used for rank-1 problems at each step in the deflation)")
            end
            x_current[:,t] .= x_output

            if theIter == 1
                weights[t] = lambda_partial
            end
        end

        ofv_prev = ofv_overall
        ofv_overall = tr(x_current'*Sigma*x_current) 
        if theIter == 1
            ofv_prev = ofv_overall
        end

        violation = sum(abs.(x_current'*x_current .- Diagonal(ones(r))))
        violation = max(violation, 1e-7)
        stepSize = (theIter < fastPeriod ? 0.01 : 0.05)*(theIter < slowPeriod ? violation : ofv_overall/violation) #OLD RULE

        if verbose 
            if numIters <= 25 || theIter % 10 == 0
                Printf.@printf(" %10d | %20.3f | %25.2e | %10.3f \n", theIter, ofv_overall / tr(Sigma), violation, (time() - start_time) )
            end
        end

        if violation < violation_tolerance || (theIter==numIters && ofv_best < 0) # Find the best solution which satisfies the feasibility tolerance, or report the last solution if nothing found so far
            ofv_current=tr(x_current'*Sigma*x_current) 
            if ofv_best < ofv_current
                x_best = x_current
                ofv_best= ofv_current
            end
        end
       
        if abs(ofv_prev - ofv_overall) < 1e-8 && violation < violation_tolerance #theIter > fastPeriod
            if ofv_best < 0 #If no solution found, return the last one
                x_best = x_current
                ofv_best= tr(x_current'*Sigma*x_current) 
            end
            break
        end
    end
    runtime = (time() - start_time)
    violation_best = sum(abs.(x_best'*x_best .- Diagonal(ones(r))))
    ofv_best = tr(x_best'*Sigma*x_best) 

return ofv_best, violation_best, runtime, x_best
end



