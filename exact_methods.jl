function solveExact_direct(Sigma::Array{Float64, 2}, k::Int64, r::Int64; targetSparsity::Array{Int64,1}=k*ones(r), timeLimit::Float64=600.0, theGap::Float64=1e-2, use_ell1::Bool=false, warmStart::Array{Float64, 2}=zeros(n,r), useUB::Bool=false,
    withFrob::Bool=false, frobPenalty::Float64=1e-3, verbose::Bool=true, totalViol::Float64=1e-4)
# Note: only use a target sparsity for each PC if specified by the user, otherwise make it a redundant constraint

    n = size(Sigma, 1)
    exact_model=Model(()->Gurobi.Optimizer(GRB_ENV))

    if !verbose
        set_silent(exact_model)
    end

    @variable(exact_model, z[1:n, 1:r], Bin)
    @variable(exact_model, U[1:n, 1:r])
    @variable(exact_model, theta)
    set_optimizer_attribute(exact_model, "NonConvex", 2)
    set_optimizer_attribute(exact_model, "TimeLimit", timeLimit)
    set_optimizer_attribute(exact_model, "MIPGap", theGap)
    set_optimizer_attribute(exact_model, "FuncPieceError", 1e-6)
    set_optimizer_attribute(exact_model, "FuncPieceLength", 1e-5)
    set_optimizer_attribute(exact_model, "MemLimit", 300) # Don't use >300 GB RAM; can change this parameter depending on machine specs


    # Orthogonality
    for i=1:r
        for j=i:r
            ind=1.0*(i==j)
            @constraint(exact_model, U[:,i]'*U[:,j]<=1.0*ind+totalViol/r^2)
            @constraint(exact_model, U[:,i]'*U[:,j]>=1.0*ind-totalViol/r^2) # So that the overall constraint violation is at most 10^-4
        end
    end
    @constraint(exact_model, sum(z)<=k)
    @constraint(exact_model, imposePCTarget[t=1:r], sum(z[:,t]) <= targetSparsity[t])
    @constraint(exact_model, imposeLogical[i=1:n, t=1:r], !z[i,t] => {U[i,t]==0.0})

    if useUB # Use the upper bound developed in section 2.5 of the paper to accelerate the method further
        M_its=zeros(n,r)
        for i=1:n
            for t=1:r
                M_its[i,t]=sum(sort(abs.(Sigma[i,:]), rev=true)[1:targetSparsity[t]])
            end
        end

        @variable(exact_model, mu[1:n, 1:r], Bin)
        @variable(exact_model, rho[1:n, 1:r]>=0.0)
        @constraint(exact_model, imposeGershgorin[i=1:n, t=1:r], rho[i,t]<=sum(z[j,t]*abs.(Sigma[j,i]) for j=1:n)+M_its[i,t]*(1-mu[i,t]))
        @constraint(exact_model, imposeGershgorin2[i=1:n, t=1:r], rho[i,t]>=sum(z[j,t]*abs.(Sigma[j,i]) for j=1:n)-M_its[i,t]*(1-mu[i,t]))
        @constraint(exact_model, imposeLogical3[i=1:n, t=1:r], rho[i,t]<=M_its[i,t]*mu[i,t])
        @constraint(exact_model, definemu1[t=1:r], sum(mu[:,t])==1.0)
        @constraint(exact_model, definemu2[i=1:n], sum(mu[i,:])<=1.0)
        @constraint(exact_model, linkmuandz[i=1:n, t=1:r], !(z[i,t]) => {mu[i,t]==0.0}) # Constraint which makes upper bound non-trivial
        @constraint(exact_model, linkmuandz2[i=1:n, t=1:r], mu[i,t]=> {z[i,t]==1.0})

        @constraint(exact_model, theta <= sum(rho)) # Impose upper bound on epigraph version of objective

        # Set branching priority of mu to be lower than z
        for i=1:n
            for t=1:r
            MOI.set(exact_model, Gurobi.VariableAttribute("BranchPriority"), z[i,t], 1)
            MOI.set(exact_model, Gurobi.VariableAttribute("BranchPriority"), mu[i,t], 0)
            end
        end
    end
    
    incumbent_value=0.0 
    # Do branching callback here if useUB is true, to avoid expanding so many nodes—using a generalization of the bound from the rank-one case, but treating each PC seperately to avoid solving an optimization problem

    if useUB
        cb_calls = Cint[]
        function my_callback_function(cb_data, cb_where::Cint)
            # You can reference variables outside the function as normal
            push!(cb_calls, cb_where)

            if cb_where == GRB_CB_MIP 
                inc_val= Ref{Cdouble}()
                GRBcbget(cb_data, cb_where, GRB_CB_MIP_OBJBST, inc_val)
                
               incumbent_value=max(inc_val[], incumbent_value) # Keep track of incumbent value
            end

            # You can select where the callback is run
            if cb_where != GRB_CB_MIPSOL && cb_where != GRB_CB_MIPNODE 
                return
            end
            # You can query a callback attribute using GRBcbget


              if cb_where == GRB_CB_MIPNODE
                  resultP = Ref{Cint}()
                  GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_STATUS, resultP)
                  if resultP[] != GRB_OPTIMAL
                      return  # Solution is something other than optimal.
                  end
              end
            

            Gurobi.load_callback_variable_primal(cb_data, cb_where)
            z0=zeros(n, r)
            zZero=zeros(n,r)
            zOne=zeros(n,r)
            for i in 1:n
                for t in 1:r
                    z0[i,t]=callback_value(cb_data, z[i,t])
                end
            end
            bestbound_total=0.0

            

            for t=1:r
                bestbound=0.0
                zPositive=(z0[:,t].>1-1e-6)
                zPositiveInd=findall(z0[:,t].>1-1e-6)
                zFloating=(z0[:,t].<1-1e-4).*(z0[:,t].>1e-6)
                zFloatingInd=findall((z0[:,t].<1-1e-6).*(z0[:,t].>1e-6))
                zZero[:,t].=(z0[:,t].<1e-6)
                zOne[:,t].=(z0[:,t].>1-1e-6)
                toGo=targetSparsity[t]-sum(zPositive)
                if sum(zFloating)<toGo
                    toGo=round(sum(zFloating)) #All variables fixed, and we are below sparsity budget
                end
                
  
                    for i in 1:n
                        if z0[i,t]>1e-6
                            rowbound=sum(abs.(Sigma[i,zPositive]))+sum(sort(abs.(Sigma[i, zFloating]), rev=true)[1:toGo])
                            if bestbound<rowbound
                                bestbound=rowbound
                            end
                        end
                    end
                    bestbound_total+=bestbound
            end

            if bestbound_total<incumbent_value
                con=@build_constraint(sum(z[i,t]*(zZero[i,t]) for i in 1:n for t in 1:r)+sum((1.0-z[i,t])*zOne[i,t] for i in 1:n for t in 1:r)>=1.0)
                MOI.submit(exact_model, MOI.LazyConstraint(cb_data), con)
            end
      
            return
        end
        # You _must_ set this parameter if using lazy constraints.
        MOI.set(exact_model, MOI.RawOptimizerAttribute("LazyConstraints"), 1)
        MOI.set(exact_model, Gurobi.CallbackFunction(), my_callback_function)

    end


    if use_ell1
        @variable(exact_model, U_abs[1:n, 1:r])
        @constraint(exact_model, U.<=U_abs)
        @constraint(exact_model, -U.<=U_abs)
        @constraint(exact_model, imposePCTarget2[t=1:r], sum(U_abs[:,t])<=sqrt(targetSparsity[t]))
        if (minimum(targetSparsity)-k)<=1e-4 # This corresponds to the case where k_t is not specific, in which case we default to k_t=
            @constraint(exact_model, imposePCTargetSOC[t=1:r], [sum(z[:,t])+1.0; sum(z[:,t])-1.0; 2.0*sum(U_abs[:,t])] in SecondOrderCone())
        end

    end

    @constraint(exact_model, theta >= LinearAlgebra.dot(Sigma, U*U') ) # Objective in epigraph form
    
    @objective(exact_model, Max, LinearAlgebra.dot(Sigma, U*U'))
    set_start_value.(z, warmStart)

    optimize!(exact_model)

    @show value.(U)'*value.(U)
    @show LinearAlgebra.dot(Sigma, value.(U)*value.(U)')

    return JuMP.objective_value(exact_model), JuMP.objective_bound(exact_model), MOI.get(exact_model, MOI.NodeCount()), MOI.get(exact_model, MOI.RelativeGap()), value.(U), value.(z)
end

