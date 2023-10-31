include("core_julia1p7.jl")
include("sdp_bounds.jl")
include("iterative_deflation3.jl")
include("exact_methods.jl")

include("covariance_thresholding.jl")
include("custombb/branchAndBound.jl")
include("custombb/multiComponents.jl")
include("custombb/utilities.jl")



results_template = DataFrame(DataSetName=String[], n=Int[], r=Int[], k=Int[], upperbound_sdp_perm=Real[], runtime_sdp_perm=Real[], lowerbound_sdp_disj=Real[], violation_gd_disj=Real[]
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

# Start by running results on pitprops
for r=2:3
    for theK=5:5:10
        results_run=similar(results_template, 0)

        usePSD=true
        useSOC=false
        useCuts=false
        doDisjoint=true
        run_sdp_disj=@elapsed UB_sdp_disj, x_rounded, LB_gd_disj, violation_gd_disj, z_relax = getSDPUpperBound_gd_multiple_permutation(pitprops, r, repeat([theK], r),  usePSDs=usePSD, useSOCs=useSOC, useCuts=useCuts, generateDisjointFeasible=doDisjoint)
        @show UB_sdp_disj, run_sdp_disj

       


        # Remark: only considering scenarios where we have an equal k for each PC, so that the columns line up
        push!(results_run, ["Pitprops", n, r, theK, UB_sdp_disj, run_sdp_disj, LB_gd_disj, violation_gd_disj])
        CSV.write("table2_alg1.csv", results_run, append=true)




    end
end

# # Now do wine UCI datatset
#
normwine=[[1.0       ,  0.0943969 ,  0.211545   , -0.310235  ,  0.270798  ,  0.289101 ,   0.236815 , -0.155929 ,  0.136698    , 0.546364   ,-0.0717472 ,  0.0723432 ,   0.64372   ];
[0.0943969 ,  1.0       ,  0.164045   ,  0.2885    , -0.0545751 , -0.335167 ,  -0.411007 ,  0.292977 , -0.220746    , 0.248985   ,-0.561296  , -0.36871   ,  -0.192011  ];
[0.211545  ,  0.164045  ,  1.0        ,  0.443367  ,  0.286587  ,  0.12898  ,   0.115077 ,  0.18623  ,  0.00965194  , 0.258887   ,-0.0746669 ,  0.00391123,   0.223626  ];
[-0.310235 ,   0.2885   ,   0.443367  ,   1.0      ,  -0.0833331,  -0.321113,   -0.35137 ,   0.361922,  -0.197327   ,  0.018732  , -0.273955 ,  -0.276769 ,   -0.440597 ];
[0.270798  , -0.0545751 ,  0.286587   , -0.0833331 ,  1.0       ,  0.214401 ,   0.195784 , -0.256294 ,  0.236441    , 0.19995    , 0.0553982 ,  0.0660039 ,   0.393351  ];
[0.289101  , -0.335167  ,  0.12898    , -0.321113  ,  0.214401  ,  1.0      ,   0.864564 , -0.449935 ,  0.612413    ,-0.0551364  , 0.433681  ,  0.699949  ,   0.498115  ];
[0.236815  , -0.411007  ,  0.115077   , -0.35137   ,  0.195784  ,  0.864564 ,   1.0      , -0.5379   ,  0.652692    ,-0.172379   , 0.543479  ,  0.787194  ,   0.494193  ];
[-0.155929 ,   0.292977 ,   0.18623   ,   0.361922 ,  -0.256294 ,  -0.449935,   -0.5379  ,   1.0     ,  -0.365845   ,  0.139057  , -0.26264  ,  -0.50327  ,   -0.311385 ];
[0.136698  , -0.220746  ,  0.00965194 , -0.197327  ,  0.236441  ,  0.612413 ,   0.652692 , -0.365845 ,  1.0         ,-0.0252499  , 0.295544  ,  0.519067  ,   0.330417  ];
[0.546364  ,  0.248985  ,  0.258887   ,  0.018732  ,  0.19995   , -0.0551364,  -0.172379 ,  0.139057 , -0.0252499   , 1.0        ,-0.521813  , -0.428815  ,   0.3161    ];
[-0.0717472,  -0.561296 ,  -0.0746669 ,  -0.273955 ,   0.0553982,   0.433681,    0.543479,  -0.26264 ,   0.295544   , -0.521813  ,  1.0      ,   0.565468 ,    0.236183 ];
[0.0723432 , -0.36871   ,  0.00391123 , -0.276769  ,  0.0660039 ,  0.699949 ,   0.787194 , -0.50327  ,  0.519067    ,-0.428815   , 0.565468  ,  1.0       ,   0.312761  ];
[0.64372   , -0.192011  ,  0.223626   , -0.440597  ,  0.393351  ,  0.498115 ,   0.494193 , -0.311385 ,  0.330417    , 0.3161     , 0.236183  ,  0.312761  ,   1.0	];]
normwine=reshape(normwine, (13,13))

#
n=size(normwine,1)
B=sqrt(normwine);
for r=2:3
    for theK=5:5:10
        results_run=similar(results_template, 0)

        # Solve SOC relaxation+cuts for rounding for consistency
        usePSD=true
        useSOC=false
        useCuts=false
        doDisjoint=true
        run_sdp_disj=@elapsed UB_sdp_disj, x_rounded, LB_gd_disj, violation_gd_disj, z_relax = getSDPUpperBound_gd_multiple_permutation(normwine, r, repeat([theK], r),  usePSDs=usePSD, useSOCs=useSOC, useCuts=useCuts, generateDisjointFeasible=doDisjoint)
        @show UB_sdp_disj, run_sdp_disj

      

        # Remark: only considering scenarios where we have an equal k for each PC, so that the columns line up
        push!(results_run, ["Wine", n, r, theK, UB_sdp_disj, run_sdp_disj, LB_gd_disj, violation_gd_disj])
        CSV.write("table2_alg1.csv", results_run, append=true)

    end
end


# Now do ionosphere
normionosphere=load("data/ionosphere.jld",  "normionosphere")
@show n=size(normionosphere,1)
for r=2:3
    for theK in [5, 10, 20]
        results_run=similar(results_template, 0)

        # Solve SOC relaxation+cuts for rounding for consistency
        usePSD=true
        useSOC=false
        useCuts=false
        doDisjoint=true
        run_sdp_disj=@elapsed UB_sdp_disj, x_rounded, LB_gd_disj, violation_gd_disj, z_relax = getSDPUpperBound_gd_multiple_permutation(normionosphere, r, repeat([theK], r),  usePSDs=usePSD, useSOCs=useSOC, useCuts=useCuts, generateDisjointFeasible=doDisjoint)
        @show UB_sdp_disj, run_sdp_disj

        


        # Remark: only considering scenarios where we have an equal k for each PC, so that the columns line up
        push!(results_run, ["Ionosphere", n, r, theK, UB_sdp_disj, run_sdp_disj, LB_gd_disj, violation_gd_disj])
        CSV.write("table2_alg1.csv", results_run, append=true)

    end
end



normGeographic=load("data/geography.jld",  "normGeographic")
@show n=size(normGeographic,1)
for r=2:3
    for theK in [5, 10, 20]
        results_run=similar(results_template, 0)

        usePSD=false
        useSOC=true
        useCuts=true
        doDisjoint=true
        run_sdp_disj=@elapsed UB_sdp_disj, x_rounded, LB_gd_disj, violation_gd_disj, z_relax = getSDPUpperBound_gd_multiple_permutation(normGeographic, r, repeat([theK], r),  usePSDs=usePSD, useSOCs=useSOC, useCuts=useCuts, generateDisjointFeasible=doDisjoint)
        @show UB_sdp_disj, run_sdp_disj


        # Remark: only considering scenarios where we have an equal k for each PC, so that the columns line up
        push!(results_run, ["Geography", n, r, theK, UB_sdp_disj, run_sdp_disj, LB_gd_disj, violation_gd_disj])
        CSV.write("table2_alg1.csv", results_run, append=true)

    end
end


# Now do Communities
normcommunities=load("data/communities.jld", "normCommunities")
@show n=size(normcommunities,1)
for r=2:3
    for theK in [5, 10, 20]
        results_run=similar(results_template, 0)

        # Solve SOC relaxation+cuts for rounding for consistency
        usePSD=false
        useSOC=true
        useCuts=true
        doDisjoint=true
        run_sdp_disj=@elapsed UB_sdp_disj, x_rounded, LB_gd_disj, violation_gd_disj, z_relax = getSDPUpperBound_gd_multiple_permutation(normcommunities, r, repeat([theK], r),  usePSDs=usePSD, useSOCs=useSOC, useCuts=useCuts, generateDisjointFeasible=doDisjoint)
        @show UB_sdp_disj, run_sdp_disj

    


        # Remark: only considering scenarios where we have an equal k for each PC, so that the columns line up
        push!(results_run, ["communities", n, r, theK, UB_sdp_disj, run_sdp_disj, LB_gd_disj, violation_gd_disj])
        CSV.write("table2_alg1.csv", results_run, append=true)

    end
end


#Now do Arrythmia, switching to not adding cutting planes
normarrhythmia=load("data/arrhythmia.jld", "normArrhythmia");
@show size(normarrhythmia,1)
for i=1:size(normarrhythmia,1)
      for j=1:size(normarrhythmia,2)
            if isnan(normarrhythmia[i,j])
                  normarrhythmia[i,j]=0.0
            end
      end
end                                          # N.b. can also remove these columns apriori #we don't gain anything from keeping extra digits after this
B=sqrt(normarrhythmia)
B=real(B); #imaginary part of order 1e-12, due to round-off errors
n=size(normarrhythmia, 1)

for r=2:3
    for theK in [5, 10, 20]
        results_run=similar(results_template, 0)

        # Solve SOC relaxation+cuts for rounding for consistency
        usePSD=false
        useSOC=true
        useCuts=false
        doDisjoint=true
        run_sdp_disj=@elapsed UB_sdp_disj, x_rounded, LB_gd_disj, violation_gd_disj, z_relax = getSDPUpperBound_gd_multiple_permutation(normarrhythmia, r, repeat([theK], r),  usePSDs=usePSD, useSOCs=useSOC, useCuts=useCuts, generateDisjointFeasible=doDisjoint)
        @show UB_sdp_disj, run_sdp_disj



        # Remark: only considering scenarios where we have an equal k for each PC, so that the columns line up
        push!(results_run, ["Arrhythmia", n, r, theK, UB_sdp_disj, run_sdp_disj, LB_gd_disj, violation_gd_disj])
        CSV.write("table2_alg1.csv", results_run, append=true)

    end
end

# # Beyond this will likely require the cluster, so ok to not run on laptop
micromass=load("data/microMass.jld", "normMicroMass");
n=size(micromass, 1)

for r=2:3
    for theK in [5, 10, 20]
        results_run=similar(results_template, 0)

        # Solve SOC relaxation+cuts for rounding for consistency
        usePSD=false
        useSOC=true
        useCuts=false
        doDisjoint=true
        run_sdp_disj=@elapsed UB_sdp_disj, z_rounded, LB_gd_disj, violation_gd_disj, z_relax = getSDPUpperBound_gd_multiple_permutation(micromass, r, repeat([theK], r),  usePSDs=usePSD, useSOCs=useSOC, useCuts=useCuts, generateDisjointFeasible=doDisjoint, highDim=true)
        @show UB_sdp_disj, run_sdp_disj, LB_gd_disj


        # Remark: only considering scenarios where we have an equal k for each PC, so that the columns line up
        push!(results_run, ["micromass", n, r, theK, UB_sdp_disj, run_sdp_disj, LB_gd_disj, violation_gd_disj])
        CSV.write("table2_alg1.csv", results_run, append=true)

    end
end

lung=load("data/lung.jld", "normlung");
n=size(lung, 1)

for r=2:3
    for theK in [5, 10, 20]
        results_run=similar(results_template, 0)

        # Solve SOC relaxation+cuts for rounding for consistency
        usePSD=true
        useSOC=false
        useCuts=false
        doDisjoint=true
        run_sdp_disj=@elapsed UB_sdp_disj, z_rounded, LB_gd_disj, violation_gd_disj, z_relax = getSDPUpperBound_gd_multiple_permutation(lung, r, repeat([theK], r),  usePSDs=usePSD, useSOCs=useSOC, useCuts=useCuts, generateDisjointFeasible=doDisjoint, highDim=false)
        @show UB_sdp_disj, run_sdp_disj, LB_gd_disj



        # Remark: only considering scenarios where we have an equal k for each PC, so that the columns line up
        push!(results_run, ["lung", n, r, theK, UB_sdp_disj, run_sdp_disj, LB_gd_disj, violation_gd_disj])
        CSV.write("table2_alg1.csv", results_run, append=true)

    end
end

gait=load("data/gait.jld", "normgait");
n=size(gait, 1)

for r=2:3
    for theK in [5, 10, 20]
        results_run=similar(results_template, 0)

        # Solve SOC relaxation+cuts for rounding for consistency
        usePSD=false
        useSOC=true
        useCuts=false
        doDisjoint=true
        run_sdp_disj=@elapsed UB_sdp_disj, z_rounded, LB_gd_disj, violation_gd_disj, z_relax = getSDPUpperBound_gd_multiple_permutation(gait, r, repeat([theK], r),  usePSDs=usePSD, useSOCs=useSOC, useCuts=useCuts, generateDisjointFeasible=doDisjoint, highDim=false)
        @show UB_sdp_disj, run_sdp_disj, LB_gd_disj

        # Remark: only considering scenarios where we have an equal k for each PC, so that the columns line up
        push!(results_run, ["gait", n, r, theK, UB_sdp_disj, run_sdp_disj, LB_gd_disj, violation_gd_disj])
        CSV.write("table2_alg1.csv", results_run, append=true)

    end
end

voice=load("data/voice.jld", "normvoice");
n=size(voice, 1)

for r=2:3
    for theK in [5, 10, 20]
        results_run=similar(results_template, 0)

        # Solve SOC relaxation+cuts for rounding for consistency
        usePSD=false
        useSOC=true
        useCuts=false
        doDisjoint=true
        run_sdp_disj=@elapsed UB_sdp_disj, z_rounded, LB_gd_disj, violation_gd_disj, z_relax = getSDPUpperBound_gd_multiple_permutation(voice, r, repeat([theK], r),  usePSDs=usePSD, useSOCs=useSOC, useCuts=useCuts, generateDisjointFeasible=doDisjoint, highDim=false)
        @show UB_sdp_disj, run_sdp_disj, LB_gd_disj


        # Remark: only considering scenarios where we have an equal k for each PC, so that the columns line up
        push!(results_run, ["voice", n, r, theK, UB_sdp_disj, run_sdp_disj, LB_gd_disj, violation_gd_disj])
        CSV.write("table2_alg1.csv", results_run, append=true)

    end
end


gastro=load("data/gastro.jld", "normgastro");
n=size(gastro, 1)

for r=2:3
    for theK in [5, 10, 20]
        results_run=similar(results_template, 0)

        # Solve SOC relaxation+cuts for rounding for consistency
        usePSD=false
        useSOC=true
        useCuts=false
        doDisjoint=true
        run_sdp_disj=@elapsed UB_sdp_disj, z_rounded, LB_gd_disj, violation_gd_disj, z_relax = getSDPUpperBound_gd_multiple_permutation(gastro, r, repeat([theK], r),  usePSDs=usePSD, useSOCs=useSOC, useCuts=useCuts, generateDisjointFeasible=doDisjoint, highDim=false)
        @show UB_sdp_disj, run_sdp_disj, LB_gd_disj


        # Remark: only considering scenarios where we have an equal k for each PC, so that the columns line up
        push!(results_run, ["gastro", n, r, theK, UB_sdp_disj, run_sdp_disj, LB_gd_disj, violation_gd_dis])
        CSV.write("table2_alg1.csv", results_run, append=true)

    end
end