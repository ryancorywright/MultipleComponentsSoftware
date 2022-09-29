function getSDPUpperBound_gd(Sigma::Array{Float64, 2}, k::Int64, useSOCS::Bool=false, usePSDs::Bool=true)
# Note: greedy in terms of the rounding mechanism.

    n = size(Sigma, 1)

    sdp_gd = Model(Mosek.Optimizer)
    if usePSDs
        @variable(sdp_gd, X[1:n, 1:n], PSD)
    else
        @variable(sdp_gd, X[1:n, 1:n], Symmetric)
        @constraint(sdp_gd, diagnonneg[i=1:n], X[i,i]>=0.0)
        @constraint(sdp_gd, twobytwominor[i=1:n, j=1:n], [X[i,i]+X[j,j]; X[i,i]-X[j,j]; 2.0*X[i,j]] in SecondOrderCone())
    end
    @variable(sdp_gd, z[1:n]>=0.0)
    @constraint(sdp_gd, z.<=1.0)

    @constraint(sdp_gd, sum(diag(X))==1.0)
    @constraint(sdp_gd, sum(z)<=k)
    @constraint(sdp_gd, diagX[i=1:n], X[i,i]<=z[i])
    @constraint(sdp_gd, offDiagX[i=1:n, j=1:n], (i!=j)*(X[i,j]-0.5*z[i])<=0.0)
    @constraint(sdp_gd, offDiagX2[i=1:n, j=1:n], (i!=j)*(-X[i,j]-0.5*z[i])<=0.0)
    if useSOCS
        @variable(sdp_gd, U[1:n, 1:n])
        @constraint(sdp_gd, sum(U)<=k)
        @constraint(sdp_gd, X.-U.<=0.0)
        @constraint(sdp_gd, X.+U.>=0.0)
        @constraint(sdp_gd, perspectiveRelaxation[i=1:n], [z[i]+X[i,i]; 2.0*X[:,i];z[i]-X[i,i]] in SecondOrderCone())
    end

    @objective(sdp_gd, Max, LinearAlgebra.dot(Sigma, X))

    @suppress optimize!(sdp_gd)

    #Get greedy solution and evaluate objective
    indices_gd=sortperm(value.(z), rev=true)[1:k]
    Sigma_reduced=Sigma[indices_gd, indices_gd]
    λs,xs=eigs(Sigma_reduced, nev=1)
    λ=λs[1]
    x_full=zeros(n)
    x_full[indices_gd].=xs[:,1]

    return (JuMP.objective_value(sdp_gd)), λ, x_full
end

function getSDPUpperBound_gd_permutation(Sigma::Array{Float64, 2}, k::Int64, useSOCS::Bool=false, usePSDs::Bool=true)
# Note: greedy in terms of the rounding mechanism.

    n = size(Sigma, 1)

    sdp_gd = Model(Mosek.Optimizer)
    if usePSDs
        @variable(sdp_gd, X[1:n, 1:n], PSD)
        @variable(sdp_gd, U[1:n, 1:n], PSD) # Doubly non-negative actually
        @variable(sdp_gd, Y[1:n, 1:n], PSD) # Kim et al require that this has 2x2 minors, but it's stronger to require that it's PSD

    else
        @variable(sdp_gd, X[1:n, 1:n], Symmetric)
        @variable(sdp_gd, U[1:n, 1:n], Symmetric)
        @variable(sdp_gd, Y[1:n, 1:n], Symmetric)
        @constraint(sdp_gd, twobytwominorY[i=1:n, j=1:n], [Y[i,i]+Y[j,j]; Y[i,i]-Y[j,j]; 2.0*Y[i,j]] in SecondOrderCone())
        @constraint(sdp_gd, twobytwominorX[i=1:n, j=1:n], [X[i,i]+X[j,j]; X[i,i]-X[j,j]; 2.0*X[i,j]] in SecondOrderCone())
        @constraint(sdp_gd, twobytwominorU[i=1:n, j=1:n], [U[i,i]+U[j,j]; U[i,i]-U[j,j]; 2.0*U[i,j]] in SecondOrderCone())
    end

    @variable(sdp_gd, T[1:n, 1:n]>=0.0)
    @variable(sdp_gd, rd[1:(n-1)])
    @variable(sdp_gd, td[1:n, 1:(n-1)]>=0.0)
    @variable(sdp_gd, z[1:n]>=0.0)


    @constraint(sdp_gd, X.<=Y)
    @constraint(sdp_gd, -X.<=Y)

    # Impose permutation constraints (38 in Kim et al.)
    @constraint(sdp_gd, perm1[i=1:(k-1), j=1:k], U[i,j]>=U[i+1,j])
    @constraint(sdp_gd, perm2[i=(k+1):n, j=1:n], U[i,j]==0.0)
    @constraint(sdp_gd, perm3[i=1:n, j=(k+1):n], U[i,j]==0.0)


    # Impose absolute sum of U, Y matching (39b in  Kim et al.)
    @constraint(sdp_gd, sum(U.-Y)==0.0)

    # Constraint (44) in Kim et al.
    @constraint(sdp_gd, imposemajorization1[j=1:(n-1)], sum(U[i,i] for i=1:j)>=j*rd[j]+sum(td[i,j] for i=1:n))
    # Note: typo in their paper: they indexed the last term in the constraint by j in the sum, should have been i
    @constraint(sdp_gd, imposemajorization2[i=1:n, j=1:(n-1)], X[i,i]<=td[i,j]+rd[j])


    # Constraint 49(e)-49(f) in Kim et. al.
    @constraint(sdp_gd, z.<=1.0)
    @constraint(sdp_gd, sum(z)<=k)

    # Todo: impose (50), test formulation on pitprops, edit paper to eliminate V, W (and see if we can further simplify)
    # 50(a)
    @constraint(sdp_gd, twobytwominor2[i=1:n, j=(i+1):n], [T[i,j]+T[j,i]; T[i,j]-T[j,i]; 2.0*Y[i,j]] in SecondOrderCone())
    # 50(b)
    @constraint(sdp_gd, matchdiag[i=1:n], T[i,i].==Y[i,i])
    # 50(c)
    @constraint(sdp_gd, matchz[i=1:n], z[i].==sum(T[i,j] for j=1:n))
    # 50(d)
    @constraint(sdp_gd, matchY[j=1:n], sum(T[i,j] for i=1:n).==k*Y[j,j])
    # 50(e)
    @constraint(sdp_gd, leqY[i=1:n, j=1:n], T[i,j]<=Y[j,j])


    # Impose trace constraints
    @constraint(sdp_gd, sum(X[i,i] for i=1:n)==1.0)
    @constraint(sdp_gd, sum(Y[i,i] for i=1:n)==1.0)
    @constraint(sdp_gd, sum(U[i,i] for i=1:n)==1.0)



    @objective(sdp_gd, Max, LinearAlgebra.dot(Sigma, X))

    #@suppress optimize!(sdp_gd)
    @suppress optimize!(sdp_gd)
    #Get greedy solution and evaluate objective
    x_full=zeros(n)
    λ=0.0
    if k>1
        indices_gd=sortperm(value.(z), rev=true)[1:k]
        Sigma_reduced=Sigma[indices_gd, indices_gd]
        λs,xs=eigs(Sigma_reduced, nev=1)
        λ=λs[1]
        x_full[indices_gd].=xs[:,1]
    else
        indices_gd=sortperm(value.(z), rev=true)[1:k]
        x_full[1]=1.0
        λ=Sigma[indices_gd[1], indices_gd[1]]
    end



    return (JuMP.objective_value(sdp_gd)), λ, x_full
end


function getSDPUpperBound_gd_multiple_extended(Sigma::Array{Float64, 2}, k::Int64, r::Int64, useSOCS::Bool=false, usePSDs::Bool=true, doDisjoint::Bool=false, useCuts::Bool=false, targetSparsity::Array{Int,1}=k*ones(r)) #nonconvexRounding::Bool=false,
# Note: greedy in terms of the rounding mechanism.

    n = size(Sigma, 1)

    sdp_gd = Model(Mosek.Optimizer)
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
    @constraint(sdp_gd, z.<=1.0)
    @constraint(sdp_gd, z'*ones(n).<=targetSparsity) # If we only impose an overall sparsity budget, we set targetSparsity=k for each PC
    @constraint(sdp_gd, w.<=z*ones(r))



    @constraint(sdp_gd, X.==sum(X_t[:,:,t] for t=1:r))

    @constraint(sdp_gd, imposetraceone[t=1:r], sum(diag(X_t[:,:,t]))==1.0)
    @constraint(sdp_gd, sum(z)<=k)

    @constraint(sdp_gd, diagX[i=1:n, t=1:r], X_t[i,i,t]<=z[i,t])
    @constraint(sdp_gd, offDiagX[i=1:n, j=1:n, t=1:r], (i!=j)*(X_t[i,j,t]-0.5*z[i,t])<=0.0)
    @constraint(sdp_gd, offDiagX2[i=1:n, j=1:n, t=1:r], (i!=j)*(-X_t[i,j,t]-0.5*z[i,t])<=0.0)

    if useSOCS
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

            @constraint(sdp_gd, persp2rankr[i=1:n], [k*X[i,i]+w[i]; k*X[i,i]-w[i]; 2.0*sum(U[i,:])] in SecondOrderCone())

            @constraint(sdp_gd, soc_orth[j=1:n], [(k-r+1.0)*w[j]+w[j]-X[j,j]; (k-r+1.0)*w[j]-w[j]+X[j,j]; 2.0*((1:n).!=j).*X[:,j] ] in SecondOrderCone())


            # Impose additional inequalities
            @constraint(sdp_gd, persp2_plus[i=1:n, t=1:r], [targetSparsity[t]*X_t[i,i,t]+z[i,t]; targetSparsity[t]*X_t[i,i,t]-z[i,t]; 2.0*sum(U_t[i,:,t])] in SecondOrderCone())

            @constraint(sdp_gd, perspectiveRelaxation3[j=1:n, t=1:r], [targetSparsity[t]*z[j,t]-X_t[j,j,t]; (targetSparsity[t]-2.0)*z[j,t]+X_t[j,j,t]; 2.0*((1:n).!=j).*X_t[:,j,t]] in SecondOrderCone())
            @constraint(sdp_gd, soc_orth_rankr[j=1:n, t=1:r], [(targetSparsity[t]-1.0)*z[j,t]+z[j,t]-X_t[j,j,t]; (targetSparsity[t]-1.0)*z[j,t]-z[j,t]+X_t[j,j,t]; 2.0*((1:n).!=j).*X_t[:,j,t] ] in SecondOrderCone())

            # @constraint(sdp_gd, imposeCyclic[i=1:(n-2), j=(i+1):(n-1), l=(j+1):n, t=1:r], U_t[i,j, t]+U_t[j,l,t]+U_t[l,i,t]<=(z[i,t]+z[j,t]+z[l,t])/3)


        end


    @objective(sdp_gd, Max, LinearAlgebra.dot(Sigma, X))

    optimize!(sdp_gd)
    #@suppress optimize!(sdp_gd)
    @show (JuMP.objective_value(sdp_gd))

    theX=value.(X)
    theXt=value.(X_t)
    theW=value.(w)
    # Impose PSD cuts (on each X_t seperately)
    # useCuts=true
    if useCuts
     j=0
     maxCuts=50
         while j<maxCuts
             numCuts=0
             for t=1:r
                 λ, ϕ=eigs(Matrix(theXt[:,:,t]),  nev=1, which=:SR, maxiter=10000)
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
    z_rounded=zeros(n,r)
    x_rounded=zeros(n,r)
    nodesExpanded=0.0
    z_relax=value.(z)
    rounded_ofv=0.0
    rounded_orth=0.0
    if doDisjoint
        # Solve MIO to find closest disjoint solution
        m_disj=Model(Gurobi.Optimizer)
        @variable(m_disj, z_round[1:n, 1:r], Bin)
        @constraint(m_disj, sum(z_round)<=k)
        @constraint(m_disj, imposeorth2[t=1:r], sum(z_round[:,t])>=1.0)
        # if equalSparsity #&& k<n
        #     @constraint(m_disj, imposeorth3[t=1:r], sum(z_round[:,t])<=(k/r))
        # end
        @constraint(m_disj, imposeorth[i=1:n], sum(z_round[i,:])<=1.0)
        @objective(m_disj, Max, sum(z_round[i,t]*z_relax[i,t] for i=1:n for t=1:r))
        @suppress optimize!(m_disj)
        z_rounded=value.(z_round)
        x_rounded=zeros(n,r)
        #@show findall(z_rounded.>0.0)
        for t=1:r
            indices_disj=findall(z_rounded[:,t].>1e-4)
            if size(indices_disj,1)>1
                Sigma_reduced=Sigma[indices_disj, indices_disj]
                λs,xs=eigs(Sigma_reduced, nev=1)
                x_rounded[indices_disj, t].=xs
                rounded_ofv=rounded_ofv+λs[1]
            elseif size(indices_disj,1)>0
                x_rounded[indices_disj[1],t]=1.0
                rounded_ofv=rounded_ofv+Sigma[indices_disj[1], indices_disj[1]]
            end
        end
        rounded_orth=sum(abs.(x_rounded'*x_rounded.-Diagonal(ones(r))))

    # elseif nonconvexRounding
    #     # Solve MIO to find closest binary solution
    #     m_disj=Model(Gurobi.Optimizer)
    #     @variable(m_disj, z_round[1:n, 1:r], Bin)
    #     @constraint(m_disj, sum(z_round)<=k)
    #     @constraint(m_disj, imposeorth2[t=1:r], sum(z_round[:,t])>=1.0)
    #     if equalSparsity
    #         @constraint(m_disj, imposeorth3[t=1:r], sum(z_round[:,t])>=(k/r)-1e-4)
    #     end
    #     @objective(m_disj, Max, sum(z_round[i,t]*z_relax[i,t] for i=1:n for t=1:r))
    #     @suppress optimize!(m_disj)
    #     z_rounded=value.(z_round)
    #
    #     nonconvex_round_model=Model(Gurobi.Optimizer)
    #     @variable(nonconvex_round_model, U[1:n, 1:r])
    #
    #     set_optimizer_attribute(nonconvex_round_model, "NonConvex", 2)
    #     set_optimizer_attribute(nonconvex_round_model, "TimeLimit", 600.0)
    #     set_optimizer_attribute(nonconvex_round_model, "MIPGap", 1e-2)
    #     # set_optimizer_attribute(nonconvex_round_model, "FuncPieceLength", 1e-5)
    #     # set_optimizer_attribute(nonconvex_round_model, "FuncPieceError", 1e-5)
    #     # Orthogonality
    #     for i=1:r
    #         for j=i:r
    #             ind=1.0*(i==j)
    #             @constraint(nonconvex_round_model, U[:,i]'*U[:,j]<=1.0*ind+1e-4)
    #             @constraint(nonconvex_round_model, U[:,i]'*U[:,j]>=1.0*ind-1e-4)
    #         end
    #     end
    #
    #
    #     for i=1:n
    #         for t=1:r
    #             if z_rounded[i,t]<1e-4
    #                 @constraint(nonconvex_round_model , U[i,t]==0.0)
    #             end
    #         end
    #     end
    #
    #     @objective(nonconvex_round_model, Max, LinearAlgebra.dot(Sigma, U*U'))
    #
    #     @suppress optimize!(nonconvex_round_model)
    #     x_rounded=value.(U)
    #     rounded_ofv=tr(x_rounded'*Sigma*x_rounded)
    #     rounded_orth=sum(abs.(x_rounded'*x_rounded.-Diagonal(ones(r))))
    #     nodesExpanded=MOI.get(nonconvex_round_model, MOI.NodeCount())

    end



    return (JuMP.objective_value(sdp_gd)), z_rounded, rounded_ofv, rounded_orth#, nodesExpanded
end



function getSDPUpperBound_gd_multiple_permutation(Sigma::Array{Float64, 2}, r::Int64, targetSparsity::Array{Int,1}, usePSDs::Bool=true, useSOCs::Bool=true, useCuts::Bool=false, doDisjoint::Bool=false)
# Assumption: vector targetSparsity of appropriate dimensionality. Could fix via an assert test.

    n = size(Sigma, 1)

    sdp_gd = Model(Mosek.Optimizer)
    @variable(sdp_gd, X[1:n, 1:n, 1:r])
    @variable(sdp_gd, Y[1:n, 1:n, 1:r])
    @variable(sdp_gd, w[1:n]>=0.0)
    if usePSDs
        @variable(sdp_gd, X_overall[1:n, 1:n], Symmetric)
        if r>1
            @constraint(sdp_gd, Diagonal(w).-X_overall in PSDCone())
        end
    else
        @variable(sdp_gd, X_overall[1:n, 1:n], Symmetric)
        if useSOCs
            @constraint(sdp_gd, twobytwoOrth[i=1:n, j=(i+1):n], [w[i]+w[j]-X_overall[i,i]-X_overall[j,j]; w[j]-w[i]+X_overall[i,i]-X_overall[j,j]; 2.0*X_overall[i,j]] in SecondOrderCone())
        end
    end
    @variable(sdp_gd, U[1:n, 1:n, 1:r]) # Doubly non-negative in full formulation
    # Define variable ranges
    for t=1:r
        @constraint(sdp_gd, Y[:,:,t].==Y[:,:,t]')
        if usePSDs
            @constraint(sdp_gd, X[:,:,t] in PSDCone())
            @constraint(sdp_gd, U[:,:,t] in PSDCone())
            @constraint(sdp_gd, Y[:,:,t] in PSDCone())
        end
    end

    if !(usePSDs) # O(r^3), not O(p^2 r) of these, so retain even without using the more expressive family of SOCs
        @constraint(sdp_gd, twobytwominorU[i=1:n, j=1:n, t=1:r], [U[i,i,t]+U[j,j,t]; U[i,i,t]-U[j,j,t]; 2.0*U[i,j,t]] in SecondOrderCone())
    end
    if useSOCs
        @constraint(sdp_gd, twobytwominorX[i=1:n, j=1:n, t=1:r], [X[i,i,t]+X[j,j,t]; X[i,i,t]-X[j,j,t]; 2.0*X[i,j,t]] in SecondOrderCone())
        @constraint(sdp_gd, twobytwominorY[i=1:n, j=1:n, t=1:r], [Y[i,i,t]+Y[j,j,t]; Y[i,i,t]-Y[j,j,t]; 2.0*Y[i,j,t]] in SecondOrderCone())

    end

    @constraint(sdp_gd, X_overall.==sum(X[:,:,t] for t=1:r))





    @variable(sdp_gd, T[1:n, 1:n, 1:r]>=0.0)
    @variable(sdp_gd, rd[1:(n-1), 1:r])
    @variable(sdp_gd, td[1:n, 1:(n-1), 1:r]>=0.0)
    @variable(sdp_gd, z[1:n, 1:r]>=0.0)
    @constraint(sdp_gd, w.<=1.0)
    @constraint(sdp_gd, w.<=z*ones(r))


    @constraint(sdp_gd, X.<=Y)
    @constraint(sdp_gd, -X.<=Y)

    # Impose permutation constraints (38 in Kim et al.)
    @constraint(sdp_gd, perm1[t=1:r, i=1:(targetSparsity[t]-1), j=1:(targetSparsity[t])], U[i,j,t]>=U[i+1,j,t])
    @constraint(sdp_gd, perm2[t=1:r, i=(targetSparsity[t]+1):n, j=1:n], U[i,j,t]==0.0)
    @constraint(sdp_gd, perm3[t=1:r, i=1:n, j=(targetSparsity[t]+1):n], U[i,j,t]==0.0)


    # Impose absolute sum of U, Y matching (39b in  Kim et al.)
    for t=1:r
        @constraint(sdp_gd, sum(U[:,:,t].-Y[:,:,t])==0.0)
    end

    # Constraint (44) in Kim et al.
    @constraint(sdp_gd, imposemajorization1[j=1:(n-1), t=1:r], sum(U[i,i,t] for i=1:j)>=j*rd[j,t]+sum(td[i,j,t] for i=1:n))
    # Note: typo in their paper: they indexed the last term in the constraint by j in the sum, should have been i
    @constraint(sdp_gd, imposemajorization2[i=1:n, j=1:(n-1), t=1:r], X[i,i,t]<=td[i,j,t]+rd[j,t])


    # Constraint 49(e)-49(f) in Kim et. al.
    @constraint(sdp_gd, z.<=1.0)
    @constraint(sdp_gd, imposetargetsparsity[t=1:r], sum(z[:,t])<=targetSparsity[t])

    # 50(a) in Kim et al.
    if useSOCs | usePSDs # Not implied by the PSD constraints, need to impose explicitly
        @constraint(sdp_gd, twobytwominor2[i=1:n, j=(i+1):n, t=1:r], [T[i,j,t]+T[j,i,t]; T[i,j,t]-T[j,i,t]; 2.0*Y[i,j,t]] in SecondOrderCone())
    end
    # 50(b) in Kim et al.
    @constraint(sdp_gd, matchdiag[i=1:n, t=1:r], T[i,i,t].==Y[i,i,t])
    # 50(c)
    @constraint(sdp_gd, matchz[i=1:n, t=1:r], z[i,t].==sum(T[i,j,t] for j=1:n))
    # 50(d)
    @constraint(sdp_gd, matchY[j=1:n, t=1:r], sum(T[i,j,t] for i=1:n).==targetSparsity[t]*Y[j,j,t])
    # 50(e)
    @constraint(sdp_gd, leqY[i=1:n, j=1:n, t=1:r], T[i,j,t]<=Y[j,j,t])


    # Impose trace constraints
    @constraint(sdp_gd, traceOne1[t=1:r], sum(X[i,i,t] for i=1:n)==1.0)
    @constraint(sdp_gd, traceOne2[t=1:r], sum(Y[i,i,t] for i=1:n)==1.0)
    @constraint(sdp_gd, traceOne3[t=1:r], sum(U[i,i,t] for i=1:n)==1.0)

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



    @objective(sdp_gd, Max, LinearAlgebra.dot(Sigma, X_overall))

    #@suppress optimize!(sdp_gd)
    optimize!(sdp_gd)

    theX=value.(X_overall)
    theXt=value.(X)
    theYt=value.(Y)
    theUt=value.(U)
    theW=value.(w)
    # Impose PSD cuts (on each X_t seperately)
    # useCuts=true
    if useCuts
     j=0
     maxCuts=50
         while j<maxCuts
             numCuts=0
             for t=1:r
                 λ, ϕ=eigs(Matrix(theXt[:,:,t]),  nev=1, which=:SR, maxiter=10000)
                    if real(λ[1])<=-1e-4
                         cutV=real(ϕ[:, 1])
                         @constraint(sdp_gd, Compat.dot(cutV*cutV', X[:,:,t])>=0.0)
                         numCuts+=1
                    end
                λ, ϕ=eigs(Matrix(theUt[:,:,t]),  nev=1, which=:SR, maxiter=10000)
                   if real(λ[1])<=-1e-4
                        cutV=real(ϕ[:, 1])
                        @constraint(sdp_gd, Compat.dot(cutV*cutV', U[:,:,t])>=0.0)
                        numCuts+=1
                   end
               λ, ϕ=eigs(Matrix(theYt[:,:,t]),  nev=1, which=:SR, maxiter=10000)
                  if real(λ[1])<=-1e-4
                       cutV=real(ϕ[:, 1])
                       @constraint(sdp_gd, Compat.dot(cutV*cutV', Y[:,:,t])>=0.0)
                       numCuts+=1
                  end
            end
            λ, ϕ=eigs(Diagonal(theW).-theX,  nev=1, which=:SR, maxiter=10000)
            if real(λ[1])<=-1e-4
                cutV=real(ϕ[:, 1])
                @constraint(sdp_gd, Compat.dot(cutV*cutV', Diagonal(w).-X_overall)>=0.0)
                numCuts+=1
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
    z_rounded=zeros(n,r)
    x_rounded=zeros(n,r)
    nodesExpanded=0.0
    z_relax=value.(z)
    rounded_ofv=0.0
    rounded_orth=0.0
    if doDisjoint
        # Solve MIO to find closest disjoint solution
        m_disj=Model(Gurobi.Optimizer)
        @variable(m_disj, z_round[1:n, 1:r], Bin)
        @constraint(m_disj, target[t=1:r], sum(z_round[i,t] for i=1:n)<=targetSparsity[t])
        @constraint(m_disj, imposeorth2[t=1:r], sum(z_round[:,t])>=1.0)
        @constraint(m_disj, imposeorth[i=1:n], sum(z_round[i,:])<=1.0)
        @objective(m_disj, Max, sum(z_round[i,t]*z_relax[i,t] for i=1:n for t=1:r))
        @suppress optimize!(m_disj)
        z_rounded=value.(z_round)
        x_rounded=zeros(n,r)
        #@show findall(z_rounded.>0.0)
        for t=1:r
            indices_disj=findall(z_rounded[:,t].>1e-4)
            if size(indices_disj,1)>1
                Sigma_reduced=Sigma[indices_disj, indices_disj]
                λs,xs=eigs(Sigma_reduced, nev=1)
                x_rounded[indices_disj, t].=xs[:,1]
                rounded_ofv=rounded_ofv+λs[1]
            elseif size(indices_disj,1)>0
                x_rounded[indices_disj[1],t]=1.0
                rounded_ofv=rounded_ofv+Sigma[indices_disj[1], indices_disj[1]]
            end
        end
        rounded_orth=sum(abs.(x_rounded'*x_rounded.-Diagonal(ones(r))))
    end


    return (JuMP.objective_value(sdp_gd)), x_rounded, rounded_ofv, rounded_orth, z_relax
end


# function getSDPUpperBound_gd_multiple_compact(Sigma::Array{Float64, 2}, k::Int64, r::Int64, useSOCS::Bool=false, usePSDs::Bool=true)
# # Note: greedy in terms of the rounding mechanism.
#
#     n = size(Sigma, 1)
#
#     sdp_gd = Model(Mosek.Optimizer)
#     @variable(sdp_gd, z[1:n]>=0.0) # Only working with the sum anyway
#
#     if usePSDs
#         @variable(sdp_gd, X[1:n, 1:n], PSD)
#         if r>1
#             @constraint(sdp_gd, I.-X in PSDCone())
#         end
#
#     else
#         @variable(sdp_gd, X[1:n, 1:n], Symmetric)
#         @constraint(sdp_gd, diagnonneg[i=1:n], X[i,i]>=0.0)
#         @constraint(sdp_gd, twobytwominor[i=1:n, j=1:n], [X[i,i]+X[j,j]; X[i,i]-X[j,j]; 2.0*X[i,j]] in SecondOrderCone())
#         @constraint(sdp_gd, twobytwominor2[i=1:n, j=1:n], [z[i]+z[j]-X[i,i]-X[j,j]; (X[i,i]-X[j,j])-(z[i]-z[j]); 2.0*X[i,j]] in SecondOrderCone())
#         @constraint(sdp_gd, twobytwoOrth[i=1:n, j=1:n], [2.0-X[i,i]-X[j,j]; X[i,i]-X[j,j]; 2.0*X[i,j]] in SecondOrderCone())
#
#     end
#
#
#
#     @constraint(sdp_gd, z.<=r)
#
#     @constraint(sdp_gd, imposetracer, sum(diag(X[:,:]))==r)
#     @constraint(sdp_gd, sum(z)<=k)
#     @constraint(sdp_gd, diagX[i=1:n], X[i,i]<=z[i])
#     @constraint(sdp_gd, offDiagX[i=1:n, j=1:n], (i!=j)*(X[i,j]-0.5*z[i])<=0.0)
#     @constraint(sdp_gd, offDiagX2[i=1:n, j=1:n], (i!=j)*(-X[i,j]-0.5*z[i])<=0.0)
#
#     if useSOCS # Some redundancy with extended formulation case
#             @variable(sdp_gd,  U[1:n, 1:n])
#             @constraint(sdp_gd, X.-U.<=0.0)
#             @constraint(sdp_gd, X.+U.>=0.0)
#             @constraint(sdp_gd, persp2[i=1:n], [k*X[i,i]+z[i]; k*X[i,i]-z[i]; 2.0*sum(U[i,:])] in SecondOrderCone())
#             @constraint(sdp_gd, perspectiveRelaxation3[i=1:n], [r*z[i]+X[i,i]; 2.0*X[:,i];r*z[i]-X[i,i]] in SecondOrderCone())
#     end
#
#
#
#     @objective(sdp_gd, Max, LinearAlgebra.dot(Sigma, X))
#
#
#     #@suppress optimize!(sdp_gd)
#     @suppress optimize!(sdp_gd)
#
#
#     return (JuMP.objective_value(sdp_gd)), value.(z), value.(X)
# end
#
# function getSDPUpperBound_multiple_ell1(Sigma::Array{Float64, 2}, k::Int64, r::Int64, targetSparsity::Array{Int64,1}=k*ones(r))
# # Note: greedy in terms of the rounding mechanism.
#
#     n = size(Sigma, 1)
#
#     grb_mult = Model(Gurobi.Optimizer)
#     set_optimizer_attribute(grb_mult, "NonConvex", 2) # Non-convex in objective, convex constraints
#     set_optimizer_attribute(grb_mult, "TimeLimit", 600.0) # Non-convex in objective, convex constraints
#
#     @variable(grb_mult, U[1:n, 1:r])
#
#     @constraint(grb_mult, imposeell1[t=1:r], [sqrt(targetSparsity[t]); U[:,t]] in MOI.NormOneCone(n+1))
#     @constraint(grb_mult, [sqrt(k); vec(U)] in MOI.NormOneCone(n*r+1))
#     @constraint(grb_mult, imposeell2[t=1:r], [1.0; U[:,t]] in MOI.SecondOrderCone(n+1))
#     @constraint(grb_mult, imposeparb_1[t=1:r, s=1:(t-1)], [1.0; 1.0; U[:,t]-U[:,s]] in MOI.RotatedSecondOrderCone(n+2))
#     @constraint(grb_mult, imposeparb_2[t=1:r, s=1:(t-1)], [1.0; 1.0; U[:,t]+U[:,s]] in MOI.RotatedSecondOrderCone(n+2))
#     @objective(grb_mult, Max, LinearAlgebra.dot(Sigma, U*U'))
#
#
#     # # Alternative approach: use SOS-2 formulation
#     # numKnots=40
#     # thetas=zeros(n,r)
#     # sigma_mat=eigen(Sigma)
#     # gammas=zeros(n, r, 2*numKnots+1)
#     # for i=1:n
#     #     for j=1:r
#     #         thetas[i,j]=sqrt(sum(sort((sigma_mat.vectors[:,i]).^2, rev=true)[1:targetSparsity[j]]))
#     #         for l=1:(2*numKnots+1)
#     #             gammas[i,j,l]=(-1.0+2.0*(l/(2.0*numKnots+1)))*thetas[i,j]
#     #         end
#     #     end
#     # end
#
#
#     # @variable(grb_mult, eta[i=1:n,j=1:r,l=1:(2*numKnots+1)], Bin)
#     # @variable(grb_mult, g[i=1:n, j=1:r])
#     # @constraint(grb_mult, defineg[i=1:n, j=1:r], g[i,j]==sigma_mat.vectors[:,i]'*U[:,j])
#     # @constraint(grb_mult, defineg2[i=1:n, j=1:r], g[i,j]==sum(gammas[i,j,l]*eta[i,j,l] for l=1:(2*numKnots+1)))
#     # @constraint(grb_mult, boundg1[i=1:n, j=1:r], g[i,j]<=thetas[i,j])
#     # @constraint(grb_mult, boundg2[i=1:n, j=1:r], -g[i,j]<=thetas[i,j])
#     #
#     #
#     # @constraint(grb_mult, imposeSOS[i=1:n, j=1:r], eta[i,j,:] in MOI.SOS2(collect(1:(2*numKnots+1))))
#     # @variable(grb_mult, psi[i=1:n, j=1:r])
#     # @constraint(grb_mult, definePsi[i=1:n, j=1:r], psi[i,j]==sum(gammas[i,j,l]^2*eta[i,j,l] for l=1:(2*numKnots+1)))
#     # @objective(grb_mult, Max, sum(sigma_mat.values[i]*sum(psi[i,t] for t=1:r) for i=1:n))
#
#
#     # Also try imposing orthogonality here explicitly, because why not
#     # if r>1 # Don't need this in the rank-1 case, we already have the bound then
#     #     for i=2:r
#     #         for j=1:(i-1)
#     #             @constraint(grb_mult, U[:,i]'*U[:,j]<=1e-4)
#     #             @constraint(grb_mult, U[:,i]'*U[:,j]>=-1e-4)
#     #         end
#     #     end
#     # end
#
#     @objective(grb_mult, Max, LinearAlgebra.dot(Sigma, U*U'))
#     @suppress optimize!(grb_mult)
#
#
#     return (JuMP.objective_bound(grb_mult)), MOI.get(grb_mult, MOI.NodeCount()) ,value.(U) # Need bound rather than value, since may not actually solve the problem
# end
