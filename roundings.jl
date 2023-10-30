using JuMP, Gurobi, TSVD, LinearAlgebra

#Perform rounding of the solution of the relaxation z_relax 
#Constrained to return PCs with disjoint support 
#Objective is a mixture of minimizing distance to z_relax and maximizing objective estimate (via UB)
function disjointRounding(Sigma, r, targetSparsity, z_relax; 
        λUB::Real = 1, 
        verbose::Bool = false,
        timeLimit::Real=600.0)

    n = size(Sigma, 1)
    m_disj = Model(()->Gurobi.Optimizer(GRB_ENV))
    set_optimizer_attribute(m_disj, "TimeLimit", timeLimit)
    set_optimizer_attribute(m_disj, "MemLimit", 300) # Don't use >300 GB RAM; can change this parameter depending on machine spwecs


    if !verbose 
        set_silent(m_disj)
    end
    @variable(m_disj, z_round[1:n, 1:r], Bin)
    @constraint(m_disj, target[t=1:r], sum(z_round[i,t] for i=1:n) <= targetSparsity[t])

    @constraint(m_disj, [t=1:r], sum(z_round[:,t]) >=1.0)
    @constraint(m_disj, [i=1:n], sum(z_round[i,:]) <=1.0)

    @variable(m_disj, μ[1:n, 1:r], Bin)
    @variable(m_disj, ρ[1:n, 1:r]>=0.0)

    @constraint(m_disj, [i=1:n, t=1:r], !(μ[i,t]) => {ρ[i,t]==0.0})
    @constraint(m_disj, [i=1:n, t=1:r], μ[i,t] => {ρ[i,t] == sum(z_round[j,t]*abs.(Sigma[j,i]) for j=1:n)})
    
    @constraint(m_disj, [t=1:r], sum(μ[:,t]) == 1.0)
    @constraint(m_disj, [i=1:n], sum(μ[i,:]) <= 1.0)
    @constraint(m_disj, [i=1:n, t=1:r], !(z_round[i,t]) => {μ[i,t]==0.0}) # Constraint which makes upper bound non-trivial
    @constraint(m_disj, [i=1:n, t=1:r], μ[i,t]=> {z_round[i,t]==1.0})

    # @objective(m_disj, Max, sum(z_round[i,t]*z_relax[i,t] for i=1:n for t=1:r))
    @objective(m_disj, Max, sum(z_round[i,t]*z_relax[i,t] for i=1:n for t=1:r) + λUB*sum(ρ) / tr(Sigma))

    @suppress optimize!(m_disj)

    z_rounded=value.(z_round)
    x_rounded=zeros(n,r)
    for t=1:r
        indices_disj = findall(z_rounded[:,t] .> 1e-4)
        if length(indices_disj) > 1
            try
                xs, λs, _ = tsvd(Sigma[indices_disj, indices_disj],  1)
                x_rounded[indices_disj, t].=xs[:,1]
            catch # In case of LAPACK exception, use eigs here
                λs, xs =eigs(Sigma[indices_disj, indices_disj],  nev=1, which=:LR, maxiter=50000, tol=1e-12)
                x_rounded[indices_disj, t].=xs[:,1]
            end
        elseif length(indices_disj) > 0
            x_rounded[indices_disj[1],t] = 1.0
        end
    end

    return x_rounded, z_rounded
end
