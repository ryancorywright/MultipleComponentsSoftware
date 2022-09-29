include("core_julia1p7.jl")
include("sdp_bounds.jl")
include("iterative_deflation3.jl")
include("exact_methods.jl")
include("custombb/branchAndBound.jl")
include("custombb/multiComponents.jl")
include("custombb/utilities.jl")
using DelimitedFiles

results_template = DataFrame(SeedNum=Int[], prop_overlap=Real[], kTrue=Int[], rTrue=Int[], n=Int[], k=Int[],
upperbound_sdp_perm=Real[], lowerbound_sdp_disj=Real[], violation_gd_disj=Real[], TPR_disj=Real[], FPR_disj=Real[], runtime_disj=Real[],
upperbound_grb_fixed=Real[], lowerbound_grb_fixed=Real[], nodes_grb_fixed=Real[], viol_grb_fixed=Real[], TPR_GRB=Real[], FPR_GRB=Real[], runtime_grb=Real[],
lb_am_fixed=Real[], viol_am=Real[], TPR_am=Real[], FPR_am=Real[], runtime_am=Real[],
lb_bb=Real[], viol_bb=Real[], TPR_bb=Real[], FPR_bb=Real[], runtime_bb=Real[]
)


n=50
seeds=1:20
ks=1:50
k_true=20
prop_overlap=0.5 # Note: you will need to also set prop_overlap=0.05, 0.5 to fully reproduce results.
                  # Also note that prop_overlap is really a misnomer, to compute the real prop overlap you need to take 1-prop.
                  # Note: to fully reproduce results, you will also need to run method of Hein and Buhler and Zou et al yourself. You can write the test matrices to csv using a commented out line as a first step towards doing that.
numIters=200
runtime_grb=300.0 # We also tried running Gurobi, but didn't include it since it barely improved upon disjoint rounding

    # Remark: we generated the test instances using Julia 1.3. If you use Julia 1.7 or 1.8 you get different random matrices, and therefore slightly different solutions.

violationPenalty=-1e1
r=2
#
    for ARG in ARGS # Note: this script assumes you are running on a cluster which can enumerate arguments from 1-1000, as otherwise it will take a long time to run
                     # To run one instance you can change ARGS to e.g. ["21"]
        array_num = parse(Int, ARG)
        theSeed=seeds[(array_num-1)%20+1] # 20 is the number of seeds considered
        theK=ks[(array_num-1)÷20+1]
        Random.seed!(theSeed)
        U=rand(n, n)
        snr=2.0
        # Construct x_1 and x_2
        x_1=zeros(n)
        x_1[1:k_true].=1.0
        x_2=zeros(n)
        x_2[Int(k_true*prop_overlap+1):Int(k_true*prop_overlap+k_true)].=1.0

        Sigma_test=U'*U+snr*x_1*x_1'+snr*x_2*x_2'
        # write matrix to csv
        #writedlm("sigma_test_0p95_"*string(array_num)*"_.csv", Sigma_test, ", ")
    #end
    # Create test matrix from random seed
    #    Now run everything, and then write to csv
        results_run=similar(results_template, 0)

        # Solve PSD relaxation and then round it, to obtain the feasible solution
        usePSD=true
        useSOC=false
        useCuts=false
        doDisjoint=true
        run_sdp_disj=@elapsed UB_sdp_disj, x_rounded, LB_gd_disj, violation_gd_disj, z_relax = getSDPUpperBound_gd_multiple_permutation(Sigma_test, r, repeat([theK], r),  usePSD, useSOC, useCuts, doDisjoint)
        z_rounded=1.0*(abs.(x_rounded).>1e-4)
        @show UB_sdp_disj, run_sdp_disj
        @show x_rounded
        # Measure disjoint TPR and FPR
        TPR_disj=max(sum(abs.(x_rounded[:,1].*x_1).>1e-8)+sum(abs.(x_rounded[:,2].*x_2).>1e-8), sum(abs.(x_rounded[:,1].*x_2).>1e-8)+sum(abs.(x_rounded[:,2].*x_1).>1e-8))/(k_true*r)
        FPR_disj=min(sum(abs.(x_rounded[:,1].*(ones(n).-x_1)).>1e-8)+sum(abs.(x_rounded[:,2].*(ones(n).-x_2)).>1e-8), sum(abs.(x_rounded[:,1].*(ones(n).-x_2)).>1e-8)+sum(abs.(x_rounded[:,2].*(ones(n).-x_1)).>1e-8))/(n*r-k_true*r)


        # Next, do the same thing for GRB_exact
        indices_reduced=findall(z_relax*ones(r).>=1e-3)
        Sigma_test_reduced=Sigma_test[indices_reduced, indices_reduced]

        run_exact_reduced=@elapsed ofv_exact_reduced, bound_exact_reduced, nodes_exact_reduced, gap_exact_reduced, U_exact_reduced, z_exact_reduced=solveExact_direct(Sigma_test_reduced, theK*r, r, repeat([theK], r), runtime_grb, 1e-4, true, z_rounded[indices_reduced,:])
        ofv_exact_reduced=tr(U_exact_reduced'*Sigma_test_reduced*U_exact_reduced)
        U_full_grb=zeros(n,r)
        U_full_grb[indices_reduced,:].=U_exact_reduced
        violation_exact_reduced=sum(abs.(U_exact_reduced'*U_exact_reduced.-Diagonal(ones(r))))
        # Measure GRB TPR and FPR
        TPR_GRB=max(sum(abs.(U_full_grb[:,1].*x_1).>1e-8)+sum(abs.(U_full_grb[:,2].*x_2).>1e-8), sum(abs.(U_full_grb[:,1].*x_2).>1e-8)+sum(abs.(U_full_grb[:,2].*x_1).>1e-8))/(k_true*r)
        FPR_GRB=min(sum(abs.(U_full_grb[:,1].*(ones(n).-x_1)).>1e-8)+sum(abs.(U_full_grb[:,2].*(ones(n).-x_2)).>1e-8), sum(abs.(U_full_grb[:,1].*(ones(n).-x_2)).>1e-8)+sum(abs.(U_full_grb[:,2].*(ones(n).-x_1)).>1e-8))/(n*r-k_true*r)

        # Next, do the same thing for alternating minimization
        ofv_AM, violation_AM, runtime_AM, x_AM_reduced=findmultPCs_deflation(Sigma_test_reduced, r, repeat([theK], r), numIters, -20.0*n)
        U_full_am=zeros(n,r)
        U_full_am[indices_reduced,:].=x_AM_reduced

        TPR_AM=max(sum(abs.(U_full_am[:,1].*x_1).>1e-8)+sum(abs.(U_full_am[:,2].*x_2).>1e-8), sum(abs.(U_full_am[:,1].*x_2).>1e-8)+sum(abs.(U_full_am[:,2].*x_1).>1e-8))/(k_true*r)
        FPR_AM=min(sum(abs.(U_full_am[:,1].*(ones(n).-x_1)).>1e-8)+sum(abs.(U_full_am[:,2].*(ones(n).-x_2)).>1e-8), sum(abs.(U_full_am[:,1].*(ones(n).-x_2)).>1e-8)+sum(abs.(U_full_am[:,2].*(ones(n).-x_1)).>1e-8))/(n*r-k_true*r)

        # Finally, run the method of Bertsimas+Berk on the same dataset [greedy rounding with Global support too?]
        B=sqrt(Sigma_test)
        theProb=problem(B, Sigma_test)
        # Solve using Lauren's code and then write to csv
        run_BB=@elapsed theU, theRuntime=multiOptimalSPCA(theProb, theK, r)
        viol_BB=sum(abs.(theU*theU'.-Diagonal(ones(r))))
        ofv_BB=tr(theU*Sigma_test*theU')
        TPR_BB=max(sum(abs.(theU[1,:].*x_1).>1e-8)+sum(abs.(theU[2,:].*x_2).>1e-8), sum(abs.(theU[1,:].*x_2).>1e-8)+sum(abs.(theU[2,:].*x_1).>1e-8))/(k_true*r)
        FPR_BB=min(sum(abs.(theU[1,:].*(ones(n).-x_1)).>1e-8)+sum(abs.(theU[2,:].*(ones(n).-x_2)).>1e-8), sum(abs.(theU[1,:].*(ones(n).-x_2)).>1e-8)+sum(abs.(theU[2,:].*(ones(n).-x_1)).>1e-8))/(n*r-k_true*r)


        push!(results_run, [theSeed, prop_overlap, k_true, r, n, theK,
        UB_sdp_disj, LB_gd_disj, violation_gd_disj, TPR_disj, FPR_disj, run_sdp_disj,
        bound_exact_reduced, ofv_exact_reduced, nodes_exact_reduced, violation_exact_reduced, TPR_GRB, FPR_GRB, run_exact_reduced,
        ofv_AM, violation_AM, TPR_AM,  FPR_AM, runtime_AM,
        ofv_BB, viol_BB, TPR_BB, FPR_BB, run_BB
        ])


        CSV.write("figure2_raw.csv", results_run, append=true)

    end
