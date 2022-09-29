function solveExact_direct(Sigma::Array{Float64, 2}, k::Int64, r::Int64, targetSparsity::Array{Int64,1}=k*ones(r), timeLimit::Float64=600.0, theGap::Float64=1e-2, use_ell1::Bool=false, warmStart::Array{Float64, 2}=zeros(n,r))
# Note: only use a target sparsity for each PC if specified by the user, otherwise make it a redundant constraint

    n = size(Sigma, 1)
    exact_model=Model(Gurobi.Optimizer)
    @variable(exact_model, z[1:n, 1:r], Bin)
    @variable(exact_model, U[1:n, 1:r])
    set_optimizer_attribute(exact_model, "NonConvex", 2)
    set_optimizer_attribute(exact_model, "TimeLimit", timeLimit)
    set_optimizer_attribute(exact_model, "MIPGap", theGap)
    set_optimizer_attribute(exact_model, "FuncPieceError", 1e-6)
    set_optimizer_attribute(exact_model, "FuncPieceLength", 1e-5)

    # Orthogonality

    for i=1:r
        for j=i:r
            ind=1.0*(i==j)
            @constraint(exact_model, U[:,i]'*U[:,j]<=1.0*ind+1e-8)
            @constraint(exact_model, U[:,i]'*U[:,j]>=1.0*ind-1e-8)
        end
    end
    @constraint(exact_model, sum(z)<=k)
    @constraint(exact_model, imposePCTarget[t=1:r], sum(z[:,t])<=targetSparsity[t])
    @constraint(exact_model, imposeLogical[i=1:n, t=1:r], !z[i,t] => {U[i,t]==0.0})

    #Impose l_1 target sparsity, not implementing SOC formulation for sake of tractability
    if use_ell1
        @variable(exact_model, U_abs[1:n, 1:r])
        @constraint(exact_model, U.<=U_abs)
        @constraint(exact_model, -U.<=U_abs)
        #@constraint(exact_model, imposePCTarget2[t=1:r], sum(U_abs[:,t])<=sqrt(targetSparsity[t]))
        @constraint(exact_model, imposePCTargetSOC[t=1:r], [sum(z[:,t])+1.0; sum(z[:,t])-1.0; 2.0*sum(U_abs[:,t])] in SecondOrderCone())

    end
    @objective(exact_model, Max, LinearAlgebra.dot(Sigma, U*U'))
    set_start_value.(z, warmStart)

    optimize!(exact_model)

    return JuMP.objective_value(exact_model), JuMP.objective_bound(exact_model), MOI.get(exact_model, MOI.NodeCount()), MOI.get(exact_model, MOI.RelativeGap()), value.(U), value.(z)
end
