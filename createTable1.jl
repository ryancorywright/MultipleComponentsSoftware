include("core_julia1p7.jl")
include("sdp_bounds.jl")
include("exact_methods.jl")

results_template = DataFrame(n=Int[], r=Int[], k=Int[], k_split=Array{Int}[],
    Bound_strong=Real[], time_strong=Real[], bound_perm=Real[], time_perm=Real[],
    gurobi_bound=Real[], gurobi_nodes=Real[], gurobi_time=Real[]
)


pitprops=[[1,0.954,0.364,0.342,-0.129,0.313,0.496,0.424,0.592,0.545,0.084,-0.019,0.134];
       [0.954,1,0.297,0.284,-0.118,0.291,0.503,0.419,0.648,0.569,0.076,-0.036,0.144];
       [0.364,0.297,1,0.882,-0.148,0.153,-0.029,-0.054,0.125,-0.081,0.162,0.22,0.126];
       [0.342,0.284,0.882,1,0.22,0.381,0.174,-0.059,0.137,-0.014,0.097,0.169,0.015];
       [-0.129,-0.118,-0.148,0.22,1,0.364,0.296,0.004,-0.039,0.037,-0.091,-0.145,-0.208];
       [0.313,0.291,0.153,0.381,0.364,1,0.813,0.09,0.211,0.274,-0.036,0.024,-0.329];
       [0.496,0.503,-0.029,0.174,0.296,0.813,1,0.372,0.465,0.679,-0.113,-0.232,-0.424];
       [0.424,0.419,-0.054,-0.059,0.004,0.09,0.372,1,0.482,0.557,0.061,-0.357,-0.202];
       [0.592,0.648,0.125,0.137,-0.039,0.211,0.465,0.482,1,0.526,0.085,-0.127,-0.076];
       [0.545,0.569,-0.081,-0.014,0.037,0.274,0.679,0.557,0.526,1,-0.319,-0.368,-0.291];
       [0.084,0.076,0.162,0.097,-0.091,-0.036,-0.113,0.061,0.085,-0.319,1,0.029,0.007];
       [-0.019,-0.036,0.22,0.169,-0.145,0.024,-0.232,-0.357,-0.127,-0.368,0.029,1,0.184];
       [0.134,0.144,0.126,0.015,-0.208,-0.329,-0.424,-0.202,-0.076,-0.291,0.007,0.184,1];]
pitprops=reshape(pitprops, (13,13));
n=size(pitprops,1)
B=sqrt(pitprops);

r=2
for k in 4:2:10 
    results_run=similar(results_template, 0)



    run_extended=@elapsed UB_extended, = getSDPUpperBound_gd_multiple_extended(pitprops, r, repeat([k], r), useSOCs=true, usePSDs=true, generateDisjointFeasible=true, useCuts=false, verbose=false, k_tot=k)


    
    run_exact=@elapsed ofv_exact, bound_exact, nodes_exact, gap_exact, U_exact, z_exact=solveExact_direct(pitprops, k, r, targetSparsity=repeat([k], r), timeLimit=7200.0, theGap=1e-4, use_ell1=true, warmStart=zeros(n,r))


        push!(results_run, [n, r, k, [-1],
            UB_extended, run_extended, -1.0, -1.0,
            ofv_exact, nodes_exact, run_exact
        ])


        CSV.write("createTable1.csv", results_run, append=true)
        results_run=similar(results_template, 0)



end

r=2
for k_split in [[1,3], [2,2], [1,5], [2,4], [3,3], [1,7], [2,6], [3,5], [4,4], [1,9], [2,8], [3,7], [4,6], [5,5]]
    results_run=similar(results_template, 0)



    run_extended=@elapsed UB_extended, = getSDPUpperBound_gd_multiple_extended(pitprops, r, k_split, useSOCs=true, usePSDs=true, generateDisjointFeasible=true, useCuts=false, verbose=false, k_tot=sum(k_split))

    run_perm=@elapsed UB_perm, = getSDPUpperBound_gd_multiple_permutation(pitprops, r, k_split, usePSDs=true, useSOCs=false, useCuts=false)

    run_exact=@elapsed ofv_exact, bound_exact, nodes_exact, gap_exact, U_exact, z_exact=solveExact_direct(pitprops, sum(k_split), r, targetSparsity=k_split, timeLimit=7200.0, theGap=1e-4, use_ell1=true, warmStart=zeros(n,r))



        push!(results_run, [n, r, sum(k_split), k_split,
            UB_extended, run_extended, UB_perm, run_perm,
            ofv_exact, nodes_exact, run_exact
        ])


        CSV.write("createTable1.csv", results_run, append=true)
        results_run=similar(results_template, 0)



end


r=3
for k in 6:3:9 
    results_run=similar(results_template, 0)



    run_extended=@elapsed UB_extended, = getSDPUpperBound_gd_multiple_extended(pitprops, r, repeat([k], r), useSOCs=true, usePSDs=true, generateDisjointFeasible=true, useCuts=false, verbose=false, k_tot=k)

    run_exact=@elapsed ofv_exact, bound_exact, nodes_exact, gap_exact, U_exact, z_exact=solveExact_direct(pitprops, k, r, targetSparsity=repeat([k], r), timeLimit=7200.0, theGap=1e-4, use_ell1=true, warmStart=zeros(n,r))


        push!(results_run, [n, r, k, [-1],
            UB_extended, run_extended, -1.0, -1.0,
            ofv_exact, nodes_exact, run_exact
        ])


        CSV.write("createTable1.csv", results_run, append=true)
        results_run=similar(results_template, 0)



end

r=3
for k_split in [[1,1,4], [1,2,3], [2,2,2], [1,1,7], [1,2,6], [1,3,5], [1,4,4], [2,2,5], [2,3,4], [3,3,3]]

    results_run=similar(results_template, 0)



    run_extended=@elapsed UB_extended, = getSDPUpperBound_gd_multiple_extended(pitprops, r, k_split, useSOCs=true, usePSDs=true, generateDisjointFeasible=true, useCuts=false, verbose=false, k_tot=sum(k_split))

    run_perm=@elapsed UB_perm, = getSDPUpperBound_gd_multiple_permutation(pitprops, r, k_split, usePSDs=true, useSOCs=false, useCuts=false)

    run_exact=@elapsed ofv_exact, bound_exact, nodes_exact, gap_exact, U_exact, z_exact=solveExact_direct(pitprops, sum(k_split), r, targetSparsity=k_split, timeLimit=7200.0, theGap=1e-4, use_ell1=true, warmStart=zeros(n,r))


        push!(results_run, [n, r, sum(k_split), k_split,
            UB_extended, run_extended, UB_perm, run_perm,
            ofv_exact, nodes_exact, run_exact
        ])


        CSV.write("createTable1.csv", results_run, append=true)
        results_run=similar(results_template, 0)



end
