using JuMP, MosekTools, MathOptInterface, TSVD

include("roundings.jl")

function getSDPUpperBound_gd(Sigma::Array{Float64, 2}, k::Int64, useSOCS::Bool=false, usePSDs::Bool=true; 
        verbose::Bool = false)
# Note: greedy in terms of the rounding mechanism.

    n = size(Sigma, 1)

    sdp_gd = Model(Mosek.Optimizer)
    if !verbose
        set_silent(sdp_gd)
    end
    if usePSDs
        @variable(sdp_gd, X[1:n, 1:n], PSD)
    else
        @variable(sdp_gd, X[1:n, 1:n], Symmetric)
        @constraint(sdp_gd, diagnonneg[i=1:n], X[i,i]>=0.0)
        @constraint(sdp_gd, twobytwominor[i=1:n, j=1:n], [X[i,i]+X[j,j]; X[i,i]-X[j,j]; 2.0*X[i,j]] in SecondOrderCone()) #SOC approx of SDP conde
    end
    @variable(sdp_gd, z[1:n] >= 0.0)
    @constraint(sdp_gd, z .<= 1.0)

    @constraint(sdp_gd, sum(diag(X)) == 1.0)
    @constraint(sdp_gd, sum(z) <= k)
    @constraint(sdp_gd, diagX[i=1:n], X[i,i] <= z[i])
    @constraint(sdp_gd, offDiagX[i=1:n, j=1:n], (i!=j)*(X[i,j]-0.5*z[i])<=0.0)
    @constraint(sdp_gd, offDiagX2[i=1:n, j=1:n], (i!=j)*(-X[i,j]-0.5*z[i])<=0.0)
    if useSOCS
        # Enforcing constraint on L1 norm of X
        @variable(sdp_gd, U[1:n, 1:n])
        @constraint(sdp_gd, sum(U) <= k)
        @constraint(sdp_gd, X.-U .<=0.0)
        @constraint(sdp_gd, X.+U .>=0.0)
        # \sum_{j} X_{ij}^2 ≤ X_ii z_ii
        @constraint(sdp_gd, perspectiveRelaxation[i=1:n], [z[i]+X[i,i]; 2.0*X[:,i];z[i]-X[i,i]] in SecondOrderCone())
    end

    @objective(sdp_gd, Max, LinearAlgebra.dot(Sigma, X))

    @suppress optimize!(sdp_gd)

    #Get greedy solution and evaluate objective
    indices_gd = sortperm(value.(z), rev=true)[1:k]
    Sigma_reduced = Sigma[indices_gd, indices_gd]
    λs,xs = eigs(Sigma_reduced, nev=1)
    λ = λs[1]
    x_full = zeros(n)
    x_full[indices_gd] .= xs[:,1]

    return (JuMP.objective_value(sdp_gd)), λ, x_full
end

#T-relaxation from Kim et al. (2022)
function getSDPUpperBound_gd_permutation(Sigma::Array{Float64, 2}, k::Int64;
                        useSOCS::Bool=false, usePSDs::Bool=true,
                        verbose::Bool = true, alwaysSOCforY::Bool = true)
# Note: greedy in terms of the rounding mechanism.
    n = size(Sigma, 1)
    
    sdp_gd = Model(Mosek.Optimizer)
    if !verbose
        set_silent(sdp_gd)
    end
    # set_attribute(sdp_gd, "MSK_IPAR_NUM_THREADS", 0)

    #Note: X encodes  xx^⊤; Y encodes |X|=; U = uu^⊤, with u that majorizes |x| 
    if usePSDs
        @variable(sdp_gd, X[1:n, 1:n], PSD)
        @variable(sdp_gd, U[1:n, 1:n], PSD) # Doubly non-negative actually
    else
        @variable(sdp_gd, X[1:n, 1:n], Symmetric)
        @variable(sdp_gd, U[1:n, 1:n], Symmetric)
        @constraint(sdp_gd, twobytwominorX[i=1:n, j=1:n], [X[i,i]+X[j,j]; X[i,i]-X[j,j]; 2.0*X[i,j]] in SecondOrderCone())
        @constraint(sdp_gd, twobytwominorU[i=1:n, j=1:n], [U[i,i]+U[j,j]; U[i,i]-U[j,j]; 2.0*U[i,j]] in SecondOrderCone())
    end
    if alwaysSOCforY || !usePSDs
        @variable(sdp_gd, Y[1:n, 1:n], Symmetric)
        @constraint(sdp_gd, twobytwominorY[i=1:n, j=1:n], [Y[i,i]+Y[j,j]; Y[i,i]-Y[j,j]; 2.0*Y[i,j]] in SecondOrderCone())
    else 
        @variable(sdp_gd, Y[1:n, 1:n], PSD) # Kim et al require that this has 2x2 minors, but it's stronger to require that it's PSD
    end
    @constraint(sdp_gd, U .>= 0.0)


    # More compact formulation for (37)
    @constraint(sdp_gd, X .<= Y)
    @constraint(sdp_gd, -X .<= Y)

    # Impose permutation constraints (38)
    @constraint(sdp_gd, perm1[i=1:(k-1), j=1:k], U[i,j] >= U[i+1,j])
    @constraint(sdp_gd, perm2[i=(k+1):n, j=1:n], U[i,j] == 0.0)
    @constraint(sdp_gd, perm3[i=1:n, j=(k+1):n], U[i,j] == 0.0)

    # Impose absolute sum of U, Y matching (39b)
    @constraint(sdp_gd, sum(U .- Y)==0.0)

    # Impose trace constraints (39a) and (41)
    @constraint(sdp_gd, sum(X[i,i] for i=1:n) == 1.0)
    @constraint(sdp_gd, sum(Y[i,i] for i=1:n) == 1.0)
    @constraint(sdp_gd, sum(U[i,i] for i=1:n) == 1.0)

    # Constraint (44)
    @variable(sdp_gd, rd[1:(n-1)])
    @variable(sdp_gd, td[1:n, 1:(n-1)] >=0.0)
    @constraint(sdp_gd, imposemajorization1[j=1:(n-1)], sum(U[i,i] for i=1:j) >= j*rd[j]+sum(td[i,j] for i=1:n))
    # Note: typo in their paper: they indexed the last term in the constraint by j in the sum, should have been i
    @constraint(sdp_gd, imposemajorization2[i=1:n, j=1:(n-1)], X[i,i] <= td[i,j]+rd[j])

    # Constraint (49d-f)
    @variable(sdp_gd, z[1:n] >=0.0)
    @constraint(sdp_gd, z .<= 1.0) #49f
    @constraint(sdp_gd, sum(z) <= k) #49e
    # @constraint(sdp_gd, [i=1:n], X[i,i] <= z[i]) #49d -- implied

    # Constraint (50)
    @variable(sdp_gd, T[1:n, 1:n] >= 0.0) # Todo: impose (50), test formulation on pitprops, edit paper to eliminate V, W (and see if we can further simplify)
    @constraint(sdp_gd, twobytwominor2[i=1:n, j=(i+1):n], [T[i,j]+T[j,i]; T[i,j]-T[j,i]; 2.0*Y[i,j]] in SecondOrderCone()) # 50(a)
    @constraint(sdp_gd, matchdiag[i=1:n], T[i,i] .== Y[i,i]) # 50(b)
    @constraint(sdp_gd, matchz[i=1:n], z[i] .== sum(T[i,j] for j=1:n))  # 50(c)
    @constraint(sdp_gd, matchY[j=1:n], sum(T[i,j] for i=1:n) <= k*Y[j,j]) # 50(d) -- with inequality instead of equality
    @constraint(sdp_gd, leqY[i=1:n, j=1:n], T[i,j] <= Y[j,j]) # 50(e)
   
    @objective(sdp_gd, Max, LinearAlgebra.dot(Sigma, X))
    optimize!(sdp_gd)

    #Get greedy solution and evaluate objective
    x_full = zeros(n)
    λ = 0.0
    if k > 1
        indices_gd = sortperm(value.(z), rev=true)[1:k]
        Sigma_reduced = Sigma[indices_gd, indices_gd]
        λs, xs = eigs(Sigma_reduced, nev=1)
        λ = λs[1]
        x_full[indices_gd] .= xs[:,1]
    else
        indices_gd = sortperm(value.(z), rev=true)[1:k]
        x_full[1] = 1.0
        λ = Sigma[indices_gd[1], indices_gd[1]]
    end

    return (JuMP.objective_value(sdp_gd)), λ, x_full
end


function getSDPUpperBound_gd_multiple_extended(Sigma::Array{Float64, 2}, r::Int64, targetSparsity::Array{Int,1}; 
        useSOCs::Bool=false, usePSDs::Bool=true, 
        addDisjointConstraint::Bool=false, #Add constraint for disjoint support directly on the relaxation
        generateDisjointFeasible::Bool=false, #Generates rounded solution with disjoint support
        useCuts::Bool=false,
        verbose::Bool = true, maxCuts::Int=50,
        k_tot::Int= sum(targetSparsity)) #nonconvexRounding::Bool=false,
# Note: greedy in terms of the rounding mechanism.

    n = size(Sigma, 1)

    sdp_gd = Model(Mosek.Optimizer)
    if !verbose
        set_silent(sdp_gd)
    end
    # set_optimizer_attribute(sdp_gd, "MSK_IPAR_PRESOLVE_USE", 0)
    # set_optimizer_attribute(sdp_gd, "MSK_DPAR_INTPNT_CO_TOL_DFEAS", 1e-12)
    # set_optimizer_attribute(sdp_gd, "MSK_DPAR_INTPNT_CO_TOL_PFEAS", 1e-12)
    # set_optimizer_attribute(sdp_gd, "MSK_DPAR_INTPNT_CO_TOL_NEAR_REL", 10)
    # set_optimizer_attribute(sdp_gd, "MSK_DPAR_INTPNT_QO_TOL_MU_RED", 1e-12)
    # set_optimizer_attribute(sdp_gd, "MSK_DPAR_INTPNT_CO_TOL_MU_RED", 1e-12)

    @variable(sdp_gd, w[1:n]>=0.0)
    @constraint(sdp_gd, w.<=1.0)

    if usePSDs
        @variable(sdp_gd, X[1:n, 1:n], Symmetric)
        if r>1
            @constraint(sdp_gd, Diagonal(w).-X in PSDCone())
        # @constraint(sdp_gd, Diagonal(sum(z[:,i] for i=1:t)).-sum(X_t[:,:,i] for i=1:t)  in PSDCone()) # Could also consider imposing this?

        end
        @variable(sdp_gd, X_t[1:n, 1:n, 1:r])
        @constraint(sdp_gd, imposeSymm[t=1:r], X_t[:,:,t] .== X_t[:,:,t]')
        @constraint(sdp_gd, imposePSD[t=1:r], X_t[:,:,t] in PSDCone())

    else
        @variable(sdp_gd, X[1:n, 1:n], Symmetric)
        @variable(sdp_gd, X_t[1:n, 1:n, 1:r])
        # if sdp_light
        #     @constraint(sdp_gd, Diagonal(w).-X in PSDCone())
        #
        # else
        @constraint(sdp_gd, twobytwoOrth[i=1:n, j=(i+1):n], [w[i]+w[j]-X[i,i]-X[j,j]; w[j]-w[i]+X[i,i]-X[j,j]; 2.0*X[i,j]] in SecondOrderCone())
        # end

        @constraint(sdp_gd, imposeSymm[t=1:r], X_t[:,:,t] .== X_t[:,:,t]')
        @constraint(sdp_gd, diagnonneg[i=1:n, t=1:r], X_t[i,i,t]>=0.0)

        @constraint(sdp_gd, twobytwominor[i=1:n, j=1:n, t=1:r], [X_t[i,i,t]+X_t[j,j,t]; X_t[i,i,t]-X_t[j,j,t]; 2.0*X_t[i,j,t]] in SecondOrderCone())
    end

    @variable(sdp_gd, z[1:n,1:r]>=0.0)
    @constraint(sdp_gd, z .<=1.0)
    @constraint(sdp_gd, z'*ones(n) .<= targetSparsity) # If we only impose an overall sparsity budget, we set targetSparsity=k for each PC
    @constraint(sdp_gd, w .<= z*ones(r))

    if addDisjointConstraint #Imposes disjoint support constraints directly on the relaxation. 
        # print("---NEW (R1) Disjoint Support Constraints on the Relaxation---")
        @constraint(sdp_gd, [i=1:n], sum(z[i,:]) <= 1.0)
    end

    @constraint(sdp_gd, X.==sum(X_t[:,:,t] for t=1:r))

    @constraint(sdp_gd, imposetraceone[t=1:r], sum(diag(X_t[:,:,t]))==1.0)
    @constraint(sdp_gd, sum(z)<=k_tot)

    @constraint(sdp_gd, diagX[i=1:n, t=1:r], X_t[i,i,t]<=z[i,t])
    @constraint(sdp_gd, offDiagX[i=1:n, j=1:n, t=1:r], (i!=j)*(X_t[i,j,t]-0.5*z[i,t])<=0.0)
    @constraint(sdp_gd, offDiagX2[i=1:n, j=1:n, t=1:r], (i!=j)*(-X_t[i,j,t]-0.5*z[i,t])<=0.0)

    if useSOCs
        # Define U, U_t
        @variable(sdp_gd, U_t[1:n, 1:n, 1:r])
        @constraint(sdp_gd, X_t.-U_t.<=0.0)
        @constraint(sdp_gd, X_t.+U_t.>=0.0)

        @variable(sdp_gd,  U[1:n, 1:n])
        @constraint(sdp_gd, X.-U.<=0.0)
        @constraint(sdp_gd, X.+U.>=0.0)

        # Impose inequalities from original problem
        @constraint(sdp_gd, perspectiverankr[i=1:n], [r*X[i,i]+w[i]; r*X[i,i]-w[i]; 2.0*X[:,i]] in SecondOrderCone())
        @constraint(sdp_gd, perspectiverankone[i=1:n, t=1:r], [z[i,t]+X_t[i,i,t]; 2.0*X_t[:,i,t];z[i,t]-X_t[i,i,t]] in SecondOrderCone())

        @constraint(sdp_gd, persp2rankr[i=1:n], [k_tot*X[i,i]+w[i]; k_tot*X[i,i]-w[i]; 2.0*sum(U[i,:])] in SecondOrderCone())

        @constraint(sdp_gd, soc_orth[j=1:n], [(k_tot-r+1.0)*w[j]+w[j]-X[j,j]; (k_tot-r+1.0)*w[j]-w[j]+X[j,j]; 2.0*((1:n).!=j).*X[:,j] ] in SecondOrderCone())


        # Impose additional inequalities
        @constraint(sdp_gd, persp2_plus[i=1:n, t=1:r], [targetSparsity[t]*X_t[i,i,t]+z[i,t]; targetSparsity[t]*X_t[i,i,t]-z[i,t]; 2.0*sum(U_t[i,:,t])] in SecondOrderCone())

        @constraint(sdp_gd, perspectiveRelaxation3[j=1:n, t=1:r], [targetSparsity[t]*z[j,t]-X_t[j,j,t]; (targetSparsity[t]-2.0)*z[j,t]+X_t[j,j,t]; 2.0*((1:n).!=j).*X_t[:,j,t]] in SecondOrderCone())
        @constraint(sdp_gd, soc_orth_rankr[j=1:n, t=1:r], [(targetSparsity[t]-1.0)*z[j,t]+z[j,t]-X_t[j,j,t]; (targetSparsity[t]-1.0)*z[j,t]-z[j,t]+X_t[j,j,t]; 2.0*((1:n).!=j).*X_t[:,j,t] ] in SecondOrderCone())

        # @constraint(sdp_gd, imposeCyclic[i=1:(n-2), j=(i+1):(n-1), l=(j+1):n, t=1:r], U_t[i,j, t]+U_t[j,l,t]+U_t[l,i,t]<=(z[i,t]+z[j,t]+z[l,t])/3)
    end

    @objective(sdp_gd, Max, LinearAlgebra.dot(Sigma, X))

    optimize!(sdp_gd)
    @show (JuMP.objective_value(sdp_gd))

    theX = value.(X)
    theXt = value.(X_t)
    theW = value.(w)

    # Impose PSD cuts (on each X_t seperately)
    # useCuts=true
    if useCuts
        j=0
        while j < maxCuts
            numCuts=0
            for t=1:r
                λ, ϕ=eigs(Matrix(theXt[:,:,t]),  nev=1, which=:SR, maxiter=10000, tol=1e-12)
                    if real(λ[1])<=-1e-4
                         cutV=real(ϕ[:, 1])
                         @constraint(sdp_gd, Compat.dot(cutV*cutV', X_t[:,:,t])>=0.0)
                         numCuts+=1
                    end
            end
            λ, ϕ=eigs(Diagonal(theW).-theX,  nev=1, which=:SR, maxiter=10000)
            if real(λ[1])<=-1e-4
                cutV=real(ϕ[:, 1])
                @constraint(sdp_gd, Compat.dot(cutV*cutV', Diagonal(w).-X)>=0.0)
                numCuts+=1
            end

            if numCuts==0
                break
            end
            @suppress optimize!(sdp_gd)
            theX=value.(X)
            theXt=value.(X_t)
            theW=value.(w)
            @show (JuMP.objective_value(sdp_gd))
            j+=1
        end
    end

    # Do disjoint rounding

    z_rounded = zeros(n,r); x_rounded=zeros(n,r)
    z_relax=value.(z)
    if generateDisjointFeasible
        x_rounded, z_rounded = disjointRounding(Sigma, r, targetSparsity, z_relax)
    end
    rounded_ofv = dot(Sigma*x_rounded,x_rounded)
    rounded_orth=sum(abs.(x_rounded'*x_rounded.-Diagonal(ones(r))))


    return (JuMP.objective_value(sdp_gd)), z_rounded, rounded_ofv, rounded_orth, z_relax#, nodesExpanded
end



function getSDPUpperBound_gd_multiple_permutation(Sigma::Array{Float64, 2}, r::Int64, targetSparsity::Array{Int,1};
                usePSDs::Bool=true, #If true, use PSD matrices, otherwise use SOC approximation of SDP cone
                useSOCs::Bool=true, 
                useCuts::Bool=false, 
                breakSymmetry::Bool=false,
                UBpenalty_rounding::Real=1.,
                withFrob::Bool=false, #If true, add - λ*||Y_t||_F^2 to the objective in the relaxation
                frobPenalty::Float64=1e-4, #Penalty for the Frobenius norm term
                addDisjointConstraint::Bool=false, #Add constraint for disjoint support directly on the relaxation
                generateDisjointFeasible::Bool=false, #Generates rounded solution with disjoint support
                verbose::Bool=true, maxCuts::Int=50,
                highDim::Bool=false # If true, do not impose SOC constraint systems with O(n^2) constraints, to reduce the memory burden of the relaxation. Assume SOC=true if this parameter is active
                )
# Assumption: vector targetSparsity of appropriate dimensionality. Could fix via an assert test.

    n = size(Sigma, 1)

    sdp_gd = Model(Mosek.Optimizer)
    if !verbose 
        set_silent(sdp_gd)
    end

    @variable(sdp_gd, X[1:n, 1:n, 1:r]) #X[:,:,t] should be the t-th component U_t U_t^⊤
    @variable(sdp_gd, Y[1:n, 1:n, 1:r] >= 0) #Y[:,:,t] should encode the absolute value of X[:,:,t]
    
    @constraint(sdp_gd, [t=1:r], X[:,:,t] .== X[:,:,t]') #Symmetric matrix
    @constraint(sdp_gd, [t=1:r], Y[:,:,t] .== Y[:,:,t]') #Symmetric matrix

    # For each t, U should model u u^⊤ with u an orderered version of |x| = y
    @variable(sdp_gd, U[1:n, 1:n, 1:r] >= 0) # Doubly non-negative in full formulation
    @constraint(sdp_gd, [t=1:r], U[:,:,t] .== U[:,:,t]') #Symmetric matrix

    # Define variable ranges
    if usePSDs    
        for t=1:r
            @constraint(sdp_gd, U[:,:,t] in PSDCone())

            @constraint(sdp_gd, X[:,:,t] in PSDCone())
            @constraint(sdp_gd, Y[:,:,t] in PSDCone())
        end
    else
        @constraint(sdp_gd, twobytwominorU[i=1:n, j=1:n, t=1:r], [U[i,i,t]+U[j,j,t]; U[i,i,t]-U[j,j,t]; 2.0*U[i,j,t]] in SecondOrderCone())
        
        if useSOCs
            if !(highDim)
                @constraint(sdp_gd, twobytwominorX[i=1:n, j=1:n, t=1:r], [X[i,i,t]+X[j,j,t]; X[i,i,t]-X[j,j,t]; 2.0*X[i,j,t]] in SecondOrderCone())
                @constraint(sdp_gd, twobytwominorY[i=1:n, j=1:n, t=1:r], [Y[i,i,t]+Y[j,j,t]; Y[i,i,t]-Y[j,j,t]; 2.0*Y[i,j,t]] in SecondOrderCone())
            end
        end
    end
    
    @variable(sdp_gd, X_overall[1:n, 1:n], Symmetric)
    @constraint(sdp_gd, X_overall .==sum(X[:,:,t] for t=1:r))

    @variable(sdp_gd, w[1:n] >= 0.0)    
    @variable(sdp_gd, z[1:n, 1:r] >= 0.0)

    # w_i = min(1, ∑_t Z_{i,t} )
    @constraint(sdp_gd, w .<= 1.0)
    @constraint(sdp_gd, w .<= z*ones(r))
    # X_overall = ∑_t Xt ⪯ Diag(w)
    if usePSDs
        @constraint(sdp_gd, Diagonal(w) .- X_overall in PSDCone())
    elseif useSOCs
        if !(highDim)
            @constraint(sdp_gd, twobytwoOrth[i=1:n, j=(i+1):n], [w[i]+w[j]-X_overall[i,i]-X_overall[j,j]; w[j]-w[i]+X_overall[i,i]-X_overall[j,j]; 2.0*X_overall[i,j]] in SecondOrderCone())
        end
    end

    if addDisjointConstraint #Imposes disjoint support constraints directly on the relaxation. 
        # print("---NEW (R1) Disjoint Support Constraints on the Relaxation---")
        @constraint(sdp_gd, [i=1:n], sum(z[i,:]) <= 1.0)
    end

    # More compact formulation for (37)
    @constraint(sdp_gd, X .<= Y)
    @constraint(sdp_gd, -X .<= Y)


    # Impose permutation constraints (38 in Kim et al.)
    @constraint(sdp_gd, perm1[t=1:r, i=1:(targetSparsity[t]-1), j=1:(targetSparsity[t])], U[i,j,t] >= U[i+1,j,t])
    @constraint(sdp_gd, perm2[t=1:r, i=(targetSparsity[t]+1):n, j=1:n], U[i,j,t] == 0.0)
    @constraint(sdp_gd, perm3[t=1:r, i=1:n, j=(targetSparsity[t]+1):n], U[i,j,t] == 0.0)

    # Impose trace constraints (39a) and (41)
    @constraint(sdp_gd, traceOne1[t=1:r], sum(X[i,i,t] for i=1:n) == 1.0)
    @constraint(sdp_gd, traceOne2[t=1:r], sum(Y[i,i,t] for i=1:n) == 1.0)
    @constraint(sdp_gd, traceOne3[t=1:r], sum(U[i,i,t] for i=1:n) == 1.0)

    # Impose full sum of U, Y matching (39b in  Kim et al.)
    @constraint(sdp_gd, [t=1:r], sum(U[:,:,t].-Y[:,:,t]) == 0.0)

    # Constraint (44)
    @variable(sdp_gd, rd[1:(n-1), 1:r])
    @variable(sdp_gd, td[1:n, 1:(n-1), 1:r] >=0.0)
    @constraint(sdp_gd, imposemajorization1[j=1:(n-1), t=1:r], sum(U[i,i,t] for i=1:j) >= j*rd[j,t] + sum(td[i,j,t] for i=1:n))
    # Note: typo in their paper: they indexed the last term in the constraint by j in the sum, should have been i
    @constraint(sdp_gd, imposemajorization2[i=1:n, j=1:(n-1), t=1:r], X[i,i,t] <= td[i,j,t]+rd[j,t])

    # Constraint 49(e)-49(f) in Kim et. al.
    @constraint(sdp_gd, z .<= 1.0)
    @constraint(sdp_gd, imposetargetsparsity[t=1:r], sum(z[:,t]) <= targetSparsity[t])

    # Constraints 50(a-e)
    @variable(sdp_gd, T[1:n, 1:n, 1:r]>=0.0)
    # if useSOCs | usePSDs # Not implied by the PSD constraints, need to impose explicitly
    if !(highDim)
        @constraint(sdp_gd, twobytwominor2[i=1:n, j=(i+1):n, t=1:r], [T[i,j,t]+T[j,i,t]; T[i,j,t]-T[j,i,t]; 2.0*Y[i,j,t]] in SecondOrderCone()) #50a
    end
    # end
    @constraint(sdp_gd, matchdiag[i=1:n, t=1:r], T[i,i,t].==Y[i,i,t]) #50b
    @constraint(sdp_gd, matchz[i=1:n, t=1:r], z[i,t] .== sum(T[i,j,t] for j=1:n)) #50c
    @constraint(sdp_gd, matchY[j=1:n, t=1:r], sum(T[i,j,t] for i=1:n) .<= targetSparsity[t]*Y[j,j,t]) #50d - with inequality in case zt not exactly kt sparse
    @constraint(sdp_gd, leqY[i=1:n, j=1:n, t=1:r], T[i,j,t] <= Y[j,j,t]) #50e


    # Impose more SOCs
    @variable(sdp_gd,  U_overall[1:n, 1:n])
    @constraint(sdp_gd, X_overall.-U_overall.<=0.0)
    @constraint(sdp_gd, X_overall.+U_overall.>=0.0)
    
    # Impose inequalities from original problem
    @constraint(sdp_gd, perspectiverankr[i=1:n], [r*X_overall[i,i]+w[i]; r*X_overall[i,i]-w[i]; 2.0*X_overall[:,i]] in SecondOrderCone())
    @constraint(sdp_gd, persp2rankr[i=1:n], [sum(targetSparsity)*X_overall[i,i]+w[i]; sum(targetSparsity)*X_overall[i,i]-w[i]; 2.0*sum(U_overall[i,:])] in SecondOrderCone())
    @constraint(sdp_gd, soc_orth[j=1:n], [(sum(targetSparsity)-r+1.0)*w[j]+w[j]-X_overall[j,j]; (sum(targetSparsity)-r+1.0)*w[j]-w[j]+X_overall[j,j]; 2.0*((1:n).!=j).*X_overall[:,j] ] in SecondOrderCone())

    # Impose additional inequalities
    @constraint(sdp_gd, soc_orth_rankr[j=1:n, t=1:r], [(targetSparsity[t]-1.0)*z[j,t]+z[j,t]-X[j,j,t]; (targetSparsity[t]-1.0)*z[j,t]-z[j,t]+X[j,j,t]; 2.0*((1:n).!=j).*X[:,j,t] ] in SecondOrderCone())

    if breakSymmetry
        _, λ, = tsvd(Sigma, r+1)
        @show λgap = minimum(targetSparsity)/n*minimum([(λ[i] - λ[i+1]) for i in 1:(length(λ)-1)])
        λgap=max(1e-8, min(1e-4, λgap))
        @constraint(sdp_gd, [t=1:(r-1)], LinearAlgebra.dot(Sigma, X[:,:,t]) >= λgap + LinearAlgebra.dot(Sigma, X[:,:,t+1])) #X[:,:,t] should be the t-th component U_t U_t^⊤
    end 

    if !withFrob
        @objective(sdp_gd, Max, LinearAlgebra.dot(Sigma, X_overall))
    else 
        λ0 = max(eigmin(Sigma)*frobPenalty, 1e-6)

        @variable(sdp_gd, ρ[1:n,1:n,1:r] >= 0)
        @constraint(sdp_gd, [i=1:n, j=1:n, t=1:r], [z[i,t] + ρ[i,j,t]; z[i,t] - ρ[i,j,t]; 2*Y[i,j,t]] in SecondOrderCone())
        @constraint(sdp_gd, [i=1:n, j=1:n, t=1:r], [z[j,t] + ρ[i,j,t]; z[j,t] - ρ[i,j,t]; 2*Y[i,j,t]] in SecondOrderCone())
        @objective(sdp_gd, Max, LinearAlgebra.dot(Sigma, X_overall) - λ0*sum(ρ))
    end
    optimize!(sdp_gd)

    # @assert all(value.(ρ) .>= 0)
    theX=value.(X_overall)
    theXt=value.(X)
    theYt=value.(Y)
    theUt=value.(U)
    theW=value.(w)

    # Impose PSD cuts (on each X_t seperately)
    # useCuts=true
    if useCuts
        j=0
        while j < maxCuts
            numCuts=0
            for t=1:r
                try
                    λ, ϕ=eigs(Matrix(theXt[:,:,t]),  nev=1, which=:SR, maxiter=50000, tol=1e-12) #By default this is machine precision, which is too strict
                    # U, s, _ = tsvd(Matrix(theXt[:,:,t]),  1)
                    # λ = s[1]; ϕ = U[:,1]
                    # @show λ 
                    if real(λ[1])<=-1e-4
                        cutV=real(ϕ[:, 1])
                        @constraint(sdp_gd, LinearAlgebra.dot(cutV*cutV', X[:,:,t])>=0.0)
                        numCuts+=1
                    end

                λ, ϕ=eigs(Matrix(theUt[:,:,t]),  nev=1, which=:SR, maxiter=50000, tol=1e-12)
                # U, s, _ = tsvd(Matrix(theUt[:,:,t]),  1)
                # λ = s[1]; ϕ = U[:,1]
                if real(λ[1])<=-1e-4
                    cutV=real(ϕ[:, 1])
                    @constraint(sdp_gd, LinearAlgebra.dot(cutV*cutV', U[:,:,t])>=0.0)
                    numCuts+=1
                end

                λ, ϕ=eigs(Matrix(theYt[:,:,t]),  nev=1, which=:SR, maxiter=50000, tol=1e-12)
                # U, s, _ = tsvd(Matrix(theYt[:,:,t]),  1)
                # λ = s[1]; ϕ = U[:,1]
                if real(λ[1])<=-1e-4
                    cutV=real(ϕ[:, 1])
                    @constraint(sdp_gd, LinearAlgebra.dot(cutV*cutV', Y[:,:,t])>=0.0)
                    numCuts+=1
                end
            catch 
                println("Exception in computing EV for cut")
            end
            
            end

            try
                λ, ϕ=eigs(Diagonal(theW).-theX,  nev=1, which=:SR, maxiter=50000, tol=1e-12)
                # U, s, _ = tsvd(Diagonal(theW).-theX,  1)
                # λ = s[1]; ϕ = U[:,1]
                if real(λ[1])<=-1e-4
                    cutV=real(ϕ[:, 1])
                    @constraint(sdp_gd, LinearAlgebra.dot(cutV*cutV', Diagonal(w).-X_overall)>=0.0)
                    numCuts+=1
                end
            catch 
                println("Exception in computing EV for cut")
            end


            if numCuts==0
                break
            end
            @suppress optimize!(sdp_gd)
            theX=value.(X_overall)
            theXt=value.(X)
            theUt=value.(U)
            theYt=value.(Y)
            theW=value.(w)
            @show (JuMP.objective_value(sdp_gd))
            j+=1
        end
    end


    # Do disjoint rounding
    z_rounded = zeros(n,r); x_rounded=zeros(n,r)
    z_relax=value.(z)
    if generateDisjointFeasible
        x_rounded, z_rounded = disjointRounding(Sigma, r, targetSparsity, z_relax, λUB=UBpenalty_rounding)
    end
    rounded_ofv = dot(Sigma*x_rounded,x_rounded)
    rounded_orth = sum(abs.(x_rounded'*x_rounded.-Diagonal(ones(r))))



    return (JuMP.objective_value(sdp_gd)), x_rounded, rounded_ofv, rounded_orth, z_relax
end