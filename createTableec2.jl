include("core_julia1p7.jl")
include("sdp_bounds.jl")
include("iterative_deflation.jl")
include("exact_methods.jl")

results_template = DataFrame(n=Int[], r=Int[], k_overall=Int[], k_1=Int[], k_2=Int[], k_3=Int[],
                            UB=Real[], ofv=Real[], disj=Real[],runtime=Real[]
)


r=3



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

for k_tot in [15, 30]
    results_run=similar(results_template, 0)

    # Solve a relaxation, and then perform disjoint rounding (no k_t's)
    usePSD=true
    useSOC=false
    useCuts=false
    doDisjoint=true
    k_ts=[k_tot, k_tot, k_tot]
    run_sdp_disj=@elapsed UB_sdp_disj, z_rounded, LB_gd_disj, violation_gd_disj = getSDPUpperBound_gd_multiple_extended(pitprops, k_tot, r,  useSOC, usePSD, doDisjoint, useCuts, k_ts)


    @show UB_sdp_disj, LB_gd_disj, run_sdp_disj
    @show z_rounded
    @show k_1=sum(z_rounded[:,1])
    @show k_2=sum(z_rounded[:,2])
    @show k_3=sum(z_rounded[:,3])

    # Remark: only considering scenarios where we have an equal k for each PC, so that the columns line up
    push!(results_run, [n, r, k_tot, k_1, k_2, k_3,
    UB_sdp_disj, LB_gd_disj, violation_gd_disj, run_sdp_disj])
    CSV.write("tableec2raw.csv", results_run, append=true)
end

    # Now do ionosphere
    normionosphere=load("data/ionosphere.jld",  "normionosphere")
    @show n=size(normionosphere,1)
    n=size(normionosphere,1)
    B=sqrt(normionosphere);
    for k_tot in [15, 30]
        results_run=similar(results_template, 0)

        # Solve a relaxation, and then perform disjoint rounding (no k_t's)
        usePSD=true
        useSOC=false
        useCuts=false
        doDisjoint=true
        k_ts=[k_tot, k_tot, k_tot]
        run_sdp_disj=@elapsed UB_sdp_disj, z_rounded, LB_gd_disj, violation_gd_disj = getSDPUpperBound_gd_multiple_extended(normionosphere, k_tot, r,  useSOC, usePSD, doDisjoint, useCuts, k_ts)


        @show UB_sdp_disj, LB_gd_disj, run_sdp_disj
        @show z_rounded
        @show k_1=sum(z_rounded[:,1])
        @show k_2=sum(z_rounded[:,2])
        @show k_3=sum(z_rounded[:,3])

        # Remark: only considering scenarios where we have an equal k for each PC, so that the columns line up
        push!(results_run, [n, r, k_tot, k_1, k_2, k_3,
        UB_sdp_disj, LB_gd_disj, violation_gd_disj, run_sdp_disj])
        CSV.write("tableec2raw.csv", results_run, append=true)
    end

# Now do geographical
normGeographic=load("data/geography.jld",  "normGeographic")
@show n=size(normGeographic,1)
n=size(normGeographic,1)
B=sqrt(normGeographic);

for k_tot in [15, 30]
    results_run=similar(results_template, 0)

    # Solve a relaxation, and then perform disjoint rounding (no k_t's)
    usePSD=false
    useSOC=true
    useCuts=true
    doDisjoint=true
    k_ts=[k_tot, k_tot, k_tot]
    run_sdp_disj=@elapsed UB_sdp_disj, z_rounded, LB_gd_disj, violation_gd_disj = getSDPUpperBound_gd_multiple_extended(normGeographic, k_tot, r,  useSOC, usePSD, doDisjoint, useCuts, k_ts)


    @show UB_sdp_disj, LB_gd_disj, run_sdp_disj
    @show z_rounded
    @show k_1=sum(z_rounded[:,1])
    @show k_2=sum(z_rounded[:,2])
    @show k_3=sum(z_rounded[:,3])

    # Remark: only considering scenarios where we have an equal k for each PC, so that the columns line up
    push!(results_run, [n, r, k_tot, k_1, k_2, k_3,
    UB_sdp_disj, LB_gd_disj, violation_gd_disj, run_sdp_disj])
    CSV.write("tableec2raw.csv", results_run, append=true)
end


# Now do communities
normcommunities=load("data/communities.jld", "normCommunities")
@show n=size(normcommunities,1)
n=size(normcommunities,1)
B=sqrt(normcommunities);

for k_tot in [15, 30]
    results_run=similar(results_template, 0)

    # Solve a relaxation, and then perform disjoint rounding (no k_t's)
    usePSD=false
    useSOC=true
    useCuts=true
    doDisjoint=true
    k_ts=[k_tot, k_tot, k_tot]
    run_sdp_disj=@elapsed UB_sdp_disj, z_rounded, LB_gd_disj, violation_gd_disj = getSDPUpperBound_gd_multiple_extended(normcommunities, k_tot, r,  useSOC, usePSD, doDisjoint, useCuts, k_ts)


    @show UB_sdp_disj, LB_gd_disj, run_sdp_disj
    @show z_rounded
    @show k_1=sum(z_rounded[:,1])
    @show k_2=sum(z_rounded[:,2])
    @show k_3=sum(z_rounded[:,3])

    # Remark: only considering scenarios where we have an equal k for each PC, so that the columns line up
    push!(results_run, [n, r, k_tot, k_1, k_2, k_3,
    UB_sdp_disj, LB_gd_disj, violation_gd_disj, run_sdp_disj])
    CSV.write("tableec2raw.csv", results_run, append=true)
end
