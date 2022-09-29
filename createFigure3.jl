include("core_julia1p7.jl")
include("sdp_bounds.jl")
include("iterative_deflation.jl")
include("exact_methods.jl")

results_template = DataFrame(n=Int[], r=Int[], k_overall=Int[], k_1=Int[], k_2=Int[], k_3=Int[],
                            UB=Real[], ofv=Real[], disj=Real[],runtime_UB=Real[], runtime_LB=Real[]
)


r=3
k_tot=30

# Note: you will need to enumerate over the args between 1 and 74 for k_tot=30 and 1 and 18 for k_tot=15, and run both k_tot=15 and k_tot=30 to fully reproduce our results

for ARG in ["1"]



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

    array_of_ks=[]
    for k_2=2:min(n,14)
        for k_1=k_2:min(n, k_tot-k_2-1)
            if k_2>=k_tot-k_1-k_2
                k_3=k_tot-k_1-k_2
                @show theKs=[k_1, k_2, k_3]
                push!(array_of_ks, [k_1, k_2, k_3])
            end
        end
    end

    k_ind = parse(Int, ARG)
    # for k_1=k_2:min(n, k_tot-k_2-1)
    #     if k_2>=k_tot-k_1-k_2
    #         k_3=k_tot-k_1-k_2
    k_1=array_of_ks[k_ind][1]
    k_2=array_of_ks[k_ind][2]
    k_3=array_of_ks[k_ind][3]
    @show k_1, k_2, k_3

    numIters=200

    violationPenalty=-1e1


        @show theKs=[k_1, k_2, k_3]
            results_run=similar(results_template, 0)

            # Solve a relaxation, and then perform variable fixing
            usePSD=true
            useSOC=false
            useCuts=false
            doDisjoint=true
            run_sdp_disj=@elapsed UB_sdp_disj, x_rounded, LB_gd_disj, violation_gd_disj, z_relax = getSDPUpperBound_gd_multiple_permutation(pitprops, r, theKs,  usePSD, useSOC, useCuts, doDisjoint)
            @show UB_sdp_disj, run_sdp_disj

            indices_reduced=findall(z_relax*ones(r).>=1e-3)
            pitprops_reduced=pitprops[indices_reduced, indices_reduced]

            # Then solve via the AM heuristic

            ofv_deflation, violation_deflation, runtime_deflation, x_best=findmultPCs_deflation(pitprops_reduced, r, theKs, numIters, -1.0*n) # Normalizing so use the same penalty as in Table 3
            @show ofv_deflation, violation_deflation


            # Remark: only considering scenarios where we have an equal k for each PC, so that the columns line up
            push!(results_run, [n, r, k_tot, k_1, k_2, k_3,
            UB_sdp_disj, ofv_deflation, violation_deflation, run_sdp_disj, runtime_deflation])
            CSV.write("figure3raw.csv", results_run, append=true)


#
#
    # Now do ionosphere
    normionosphere=load("data/ionosphere.jld",  "normionosphere")
    @show n=size(normionosphere,1)
    n=size(normionosphere,1)
    B=sqrt(normionosphere);

    numIters=200

    violationPenalty=-1e1
    # Redo array of ks since pitprops has k>p
    array_of_ks=[]
    for k_2=2:min(n,14)
        for k_1=k_2:min(n, k_tot-k_2-1)
            if k_2>=k_tot-k_1-k_2
                k_3=k_tot-k_1-k_2
                #@show theKs=[k_1, k_2, k_3]
                push!(array_of_ks, [k_1, k_2, k_3])
            end
        end
    end

    #@show "Hi"
    #k_ind = parse(Int, ARG)
    # for k_1=k_2:min(n, k_tot-k_2-1)
    #     if k_2>=k_tot-k_1-k_2
    #         k_3=k_tot-k_1-k_2
    k_1=array_of_ks[k_ind][1]
    k_2=array_of_ks[k_ind][2]
    k_3=array_of_ks[k_ind][3]
    @show theKs=[k_1, k_2, k_3]
    results_run=similar(results_template, 0)
    # Solve a relaxation, and then perform variable fixing
    usePSD=true
    useSOC=false
    useCuts=false
    doDisjoint=true
    run_sdp_disj=@elapsed UB_sdp_disj, x_rounded, LB_gd_disj, violation_gd_disj, z_relax = getSDPUpperBound_gd_multiple_permutation(normionosphere, r, theKs,  usePSD, useSOC, useCuts, doDisjoint)
    @show UB_sdp_disj, run_sdp_disj

    indices_reduced=findall(z_relax*ones(r).>=1e-3)
    normionosphere_reduced=normionosphere[indices_reduced, indices_reduced]

    # Then solve via the AM heuristic

    ofv_deflation, violation_deflation, runtime_deflation, x_best=findmultPCs_deflation(normionosphere_reduced, r, theKs, numIters, -1.0*n) # Normalizing so use the same penalty as in Table 3
    @show ofv_deflation, violation_deflation


    # Remark: only considering scenarios where we have an equal k for each PC, so that the columns line up
    push!(results_run, [n, r, k_tot, k_1, k_2, k_3,
    UB_sdp_disj, ofv_deflation, violation_deflation, run_sdp_disj, runtime_deflation])
    CSV.write("figure3raw.csv", results_run, append=true)

# Now do geographical
normGeographic=load("data/geography.jld",  "normGeographic")
@show n=size(normGeographic,1)
n=size(normGeographic,1)
B=sqrt(normGeographic);

results_run=similar(results_template, 0)

numIters=200

violationPenalty=-1e1

usePSD=false
useSOC=true
useCuts=true
doDisjoint=true
run_sdp_disj=@elapsed UB_sdp_disj, x_rounded, LB_gd_disj, violation_gd_disj, z_relax = getSDPUpperBound_gd_multiple_permutation(normGeographic, r, theKs,  usePSD, useSOC, useCuts, doDisjoint)
@show UB_sdp_disj, run_sdp_disj

indices_reduced=findall(z_relax*ones(r).>=1e-3)
normGeographic_reduced=normGeographic[indices_reduced, indices_reduced]

# Then solve via the AM heuristic

ofv_deflation, violation_deflation, runtime_deflation, x_best=findmultPCs_deflation(normGeographic_reduced, r, theKs, numIters, -1.0*n) # Normalizing so use the same penalty as in Table 3
@show ofv_deflation, violation_deflation


# Remark: only considering scenarios where we have an equal k for each PC, so that the columns line up
push!(results_run, [n, r, k_tot, k_1, k_2, k_3,
UB_sdp_disj, ofv_deflation, violation_deflation, run_sdp_disj, runtime_deflation])
CSV.write("figure3raw.csv", results_run, append=true)


normcommunities=load("data/communities.jld", "normCommunities")
@show n=size(normcommunities,1)
n=size(normcommunities,1)
B=sqrt(normcommunities);

numIters=200

violationPenalty=-1e1


results_run=similar(results_template, 0)

# for k_2=2:min(n,14)
#     for k_1=k_2:min(n, k_tot-k_2-1)
#         if k_2>=k_tot-k_1-k_2
#             k_3=k_tot-k_1-k_2
#             @show theKs=[k_1, k_2, k_3]
#             results_run=similar(results_template, 0)

            # Solve a relaxation, and then perform variable fixing
            usePSD=false
            useSOC=true
            useCuts=true
            doDisjoint=true
            run_sdp_disj=@elapsed UB_sdp_disj, x_rounded, LB_gd_disj, violation_gd_disj, z_relax = getSDPUpperBound_gd_multiple_permutation(normcommunities, r, theKs,  usePSD, useSOC, useCuts, doDisjoint)
            @show UB_sdp_disj, run_sdp_disj

            indices_reduced=findall(z_relax*ones(r).>=1e-3)
            normcommunities_reduced=normcommunities[indices_reduced, indices_reduced]

            # Then solve via the AM heuristic

            ofv_deflation, violation_deflation, runtime_deflation, x_best=findmultPCs_deflation(normcommunities_reduced, r, theKs, numIters, -1.0*n) # Normalizing so use the same penalty as in Table 3
            @show ofv_deflation, violation_deflation


            # Remark: only considering scenarios where we have an equal k for each PC, so that the columns line up
            push!(results_run, [n, r, k_tot, k_1, k_2, k_3,
            UB_sdp_disj, ofv_deflation, violation_deflation, run_sdp_disj, runtime_deflation])
            CSV.write("figure3raw.csv", results_run, append=true)
end
