#results_template2 = DataFrame(Iter=Int[], ofv=Real[], viol=Real[])

function findmultPCs_deflation(Sigma::Array{Float64, 2}, r::Int64, ks::Array{Int64,1}, numIters::Int64=30, violationPenalty::Float64=-10.0, x_init::Array{Float64, 2}=zeros(size(Sigma,1), r))
    n = size(Sigma, 1)

    # find best EVs using an alternating minimization/deflation method


    # Todo: also consider allowing warm-starts (e.g. from the SDP relaxation)

    ofv_best=-1e10
    violation_best=1e10
    x_best=zeros(n,r)
    x_current=x_init # Trying this to see if it improves our results, which are sometimes not so great when we start at 0
    ofv_prev=0.0
    ofv_overall=0.0
    theLambda=0.0
    stepSize=0.0

    runtime=@elapsed for theIter=1:numIters
        theLambda=theLambda+stepSize
        if theIter==50
            stepSize=stepSize*5
        end
        # if ofv_overall<ofv_prev-1.0
        #     theLambda=theLambda*3.0/4.0
        #     stepSize=stepSize/2.0
        # end
        @show theLambda
        for t=1:r


            theLambda2=-eigmin(Sigma-theLambda*sum(x_current[:,s]*x_current[:,s]'*(s!=t) for s=1:r))+1e-4

            # Trying averaging for stability [may not work due to sharp thresholds]
            sigma_current=real.(Sigma-theLambda*sum(x_current[:,s]*x_current[:,s]'*(s!=t) for s=1:r)+theLambda2*I)

            useSOC=false
            usePSD=true
            if n>50
                useSOC=true
                usePSD=false
            end
            ofv_partial, lambda_partial, x_output=getSDPUpperBound_gd_permutation(sigma_current, ks[t], useSOC, usePSD)
            x_current[:,t].=x_output
        end
        ofv_prev=ofv_overall
        ofv_overall=tr(x_current'*Sigma*x_current)
        if theIter==1
            ofv_prev=ofv_overall
        end

        violation=sum(abs.(x_current'*x_current.-Diagonal(ones(r))))
        if theIter==1
            stepSize=0.02*ofv_overall/violation
        end
        @show ofv_overall, violation, ofv_overall+violationPenalty*violation
        if ofv_overall+violationPenalty*violation > ofv_best+violationPenalty*violation_best
           ofv_best=ofv_overall
           violation_best=violation
           x_best=x_current

        end
       # results_run=similar(results_template2, 0)
       # push!(results_run, [theIter, ofv_overall, violation])
       # CSV.write("alg2_trace_results.csv", results_run, append=true)

        if abs(ofv_prev-ofv_overall)<1e-6 && theIter>20
            @show x_current
            break
        end
    end

    @show runtime
return ofv_best, violation_best, runtime, x_best
end
