include("core_julia1p7.jl")
include("sdp_bounds.jl")
include("iterative_deflation.jl")
include("exact_methods.jl")

results_template = DataFrame(n=Int[], r=Int[], k=Int[], upperbound_sdp_perm=Real[], runtime_sdp_perm=Real[], lowerbound_sdp_disj=Real[], violation_gd_disj=Real[],
ofv_AM=Real[], violation_AM=Real[], runtime_AM=Real[],
ofv_exact=Real[], bound_exact=Real[], violation_exact=Real[], nodes_exact=Real[], time_exact=Real[], gap_exact=Real[],
ofv_deflation=Real[], viol_defl=Real[], t_defl=Real[]
)

# Note for the reader: we also consider running deflation here


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

# Solve full SDP relaxation as p=13, so cheap to do so, and then perform disjoint rounding [we will use a different relaxation in larger-scale settings]
usePSD=true
useSOC=false
useCuts=false
doDisjoint=true
run_sdp_disj=@elapsed UB_sdp_disj, x_rounded, LB_gd_disj, violation_gd_disj = getSDPUpperBound_gd_multiple_permutation(pitprops, r, repeat([theK], r),  usePSD, useSOC, useCuts, doDisjoint)
@show UB_sdp_disj, run_sdp_disj


# Then solve via the AM heuristic

ofv_best, violation_best, runtime_deflation, x_best=findmultPCs_deflation(pitprops, r, repeat([theK], r), numIters, -1.0*n, x_rounded) # Normalizing so use the same penalty as in Table 3
@show ofv_best, violation_best

# Next, try an exact method
run_exact=@elapsed ofv_exact, bound_exact, nodes_exact, gap_exact, U_exact, z_exact=solveExact_direct(pitprops, theK*r, r, repeat([theK], r), 600.0, 1e-2, true)
violation_exact=sum(abs.(U_exact'*U_exact.-Diagonal(ones(r))))
@show ofv_exact, violation_exact
# ofv_exact=-1.0
# bound_exact=-1.0
# violation_exact=-1.0
# nodes_exact=-1
# run_exact=-1.0
# gap_exact=-1.0

# Finally, run deflation (as in Mackey or Berk)
sigma_current=pitprops
x_current=zeros(n,r)
t_deflation=@elapsed for t=1:r
    useSOC=false
    usePSD=true
    ofv_partial, lambda_partial, x_output=getSDPUpperBound_gd_permutation(sigma_current, theK, useSOC, usePSD)
    x_current[:,t].=x_output
    sigma_current=(I-x_output*x_output')*sigma_current*(I-x_output*x_output')
end

ofv_deflation=tr(x_current'*pitprops*x_current)

violation_deflation=sum(abs.(x_current'*x_current.-Diagonal(ones(r))))

    # Remark: only considering scenarios where we have an equal k for each PC, so that the columns line up
push!(results_run, [n, r, theK, UB_sdp_disj, run_sdp_disj, LB_gd_disj, violation_gd_disj, ofv_best, violation_best, runtime_deflation,
    ofv_exact, bound_exact, violation_exact, nodes_exact, run_exact, gap_exact, ofv_deflation, violation_deflation, t_deflation])
CSV.write("tableec1_raw.csv", results_run, append=true)

end
end
