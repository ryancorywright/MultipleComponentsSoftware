include("core_julia1p7.jl")
include("sdp_bounds.jl")
include("iterative_deflation3.jl")
include("exact_methods.jl")

include("custombb/branchAndBound.jl")
include("custombb/multiComponents.jl")
include("custombb/utilities.jl")


results_template = DataFrame(n=Int[], r=Int[], k=Int[], upperbound_sdp_perm=Real[], runtime_sdp_perm=Real[], lowerbound_sdp_disj=Real[], violation_gd_disj=Real[],
runtime_sdp_perm_UB=Real[], lowerbound_sdp_disj_UB=Real[], violation_gd_disj_UB=Real[],
ofv_AM=Real[], violation_AM=Real[], runtime_AM=Real[],
ofv_exact=Real[], bound_exact=Real[], violation_exact=Real[], nodes_exact=Real[], time_exact=Real[], gap_exact=Real[],
ofv_exact_UB=Real[], bound_exact_UB=Real[], violation_exact_UB=Real[], nodes_exact_UB=Real[], time_exact_UB=Real[], gap_exact_UB=Real[],
ofv_exact_ws=Real[], bound_exact_ws=Real[], violation_exact_ws=Real[], nodes_exact_ws=Real[], time_exact_ws=Real[], gap_exact_ws=Real[],
ofv_exact_UB_ws=Real[], bound_exact_UB_ws=Real[], violation_exact_UB_ws=Real[], nodes_exact_UB_ws=Real[], time_exact_UB_ws=Real[], gap_exact_UB_ws=Real[]
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

numIters=100

violationPenalty=-1e1


for r=2:1:6
for theK=2:2:10
results_run=similar(results_template, 0)


# Looking into performance without the bound from section 2.5, no symmetry breaking
usePSD=true
useSOC=false
useCuts=false
doDisjoint=true
run_sdp_disj=@elapsed UB_sdp_disj, x_rounded, LB_gd_disj, violation_gd_disj = getSDPUpperBound_gd_multiple_permutation(pitprops, r, repeat([theK], r),  usePSDs=usePSD, useSOCs=useSOC, useCuts=useCuts, generateDisjointFeasible=doDisjoint, verbose=false, UBpenalty_rounding=0.0)
@show UB_sdp_disj, run_sdp_disj, LB_gd_disj


# Looking into performance with the bound from section 2.5, no symmetry breaking
usePSD=true
useSOC=false
useCuts=false
doDisjoint=true
run_sdp_disj_UB=@elapsed UB_sdp_disj_UB, x_rounded_UB, LB_gd_disj_UB, violation_gd_disj_UB = getSDPUpperBound_gd_multiple_permutation(pitprops, r, repeat([theK], r),  usePSDs=usePSD, useSOCs=useSOC, useCuts=useCuts, generateDisjointFeasible=doDisjoint, verbose=false, UBpenalty_rounding=1.0)
@show UB_sdp_disj_UB, run_sdp_disj_UB, LB_gd_disj_UB

# Compute the warm start using the version of Alg 1 with the bound from section 2.5
z_rounded=1.0*(abs.(x_rounded_UB).>1e-4)

# Then solve via the AM heuristic

# Note: need to check that useSOC=false, usePSD=true in the iterative deflation code for this to work
ofv_best, violation_best, runtime_deflation, x_best=findmultPCs_deflation(pitprops, r, repeat([theK], r), numIters=numIters, x_init=zeros(n,r), innerSolveMethod=:BB_warmstart) # Normalizing so use the same penalty as in Table 3
@show ofv_best, violation_best

# Next, try an exact method (without UB)

run_exact=@elapsed ofv_exact, bound_exact, nodes_exact, gap_exact, U_exact, z_exact=solveExact_direct(pitprops, theK*r, r, targetSparsity=repeat([theK], r), timeLimit=600.0, theGap=1e-4, use_ell1=true, warmStart=zeros(n, r), useUB=false, verbose=false)
violation_exact=sum(abs.(U_exact'*U_exact.-Diagonal(ones(r))))
@show ofv_exact, violation_exact


# Next, try an exact method(with UB)
run_exact_UB=@elapsed ofv_exact_UB, bound_exact_UB, nodes_exact_UB, gap_exact_UB, U_exact, z_exact=solveExact_direct(pitprops, theK*r, r, targetSparsity=repeat([theK], r), timeLimit=600.0, theGap=1e-4, use_ell1=true, warmStart=zeros(n, r), useUB=true, verbose=false)
violation_exact_UB=sum(abs.(U_exact'*U_exact.-Diagonal(ones(r))))
@show ofv_exact_UB, violation_exact_UB


# Next, try an exact method (without UB, with WS)
run_exact_ws=@elapsed ofv_exact_ws, bound_exact_ws, nodes_exact_ws, gap_exact_ws, U_exact_ws, z_exact_ws=solveExact_direct(pitprops, theK*r, r, targetSparsity=repeat([theK], r), timeLimit=600.0, theGap=1e-4, use_ell1=true, warmStart=z_rounded, useUB=false, verbose=false)
violation_exact_ws=sum(abs.(U_exact_ws'*U_exact_ws.-Diagonal(ones(r))))
@show ofv_exact_ws, violation_exact_ws


# Next, try an exact method(with UB, with WS)
run_exact_UB_ws=@elapsed ofv_exact_UB_ws, bound_exact_UB_ws, nodes_exact_UB_ws, gap_exact_UB_ws, U_exact_ws, z_exact_ws=solveExact_direct(pitprops, theK*r, r, targetSparsity=repeat([theK], r), timeLimit=600.0, theGap=1e-4, use_ell1=true, warmStart=z_rounded, useUB=true, verbose=false)
violation_exact_UB_ws=sum(abs.(U_exact_ws'*U_exact_ws.-Diagonal(ones(r))))
@show ofv_exact_UB_ws, violation_exact_UB_ws


#     # Remark: only considering scenarios where we have an equal k for each PC, so that the columns line up
push!(results_run, [n, r, theK, UB_sdp_disj, run_sdp_disj, LB_gd_disj, violation_gd_disj, run_sdp_disj_UB, LB_gd_disj_UB, violation_gd_disj_UB, ofv_best, violation_best, runtime_deflation,
    ofv_exact, bound_exact, violation_exact, nodes_exact, run_exact, gap_exact, ofv_exact_UB, bound_exact_UB, violation_exact_UB, nodes_exact_UB, run_exact_UB, gap_exact_UB,
    ofv_exact_ws, bound_exact_ws, violation_exact_ws, nodes_exact_ws, run_exact_ws, gap_exact_ws, ofv_exact_UB_ws, bound_exact_UB_ws, violation_exact_UB_ws, nodes_exact_UB_ws, run_exact_UB_ws, gap_exact_UB_ws])
CSV.write("tableec3_raw.csv", results_run, append=true)

end
end
