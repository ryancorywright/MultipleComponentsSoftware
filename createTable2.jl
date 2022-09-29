include("core_julia1p7.jl")
include("sdp_bounds.jl")
include("exact_methods.jl")

results_template = DataFrame(Name=String[], n=Int[], r=Int[], k=Int[], k_split=Array{Int}[],
    Bound_PSD=Real[], time_PSD=Real[], Bound_cuts=Real[], time_cuts=Real[], Bound_soc=Real[], time_soc=Real[],
    bound_highdim=Real[], time_highdim=Real[]
)


# Note to the reader: we also compute a "high dimensional" bound which wasn't included in the paper in this set of experiments. Please ignore it.

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


for r=2:3
    for k_each=5:5:10
    results_run=similar(results_template, 0)
    k_split=repeat([k_each], r)

        usePSDs=true
        useCuts=false
        useSOCs=false
        run_psd=@elapsed UB_psd, = getSDPUpperBound_gd_multiple_permutation(pitprops, r, k_split, usePSDs, useSOCs, useCuts)
        @show run_psd, UB_psd

        usePSDs=false
        useCuts=true
        useSOCs=true
        run_cuts=@elapsed UB_cuts, = getSDPUpperBound_gd_multiple_permutation(pitprops, r, k_split, usePSDs, useSOCs, useCuts)
        @show run_cuts, UB_cuts

        usePSDs=false
        useCuts=false
        useSOCs=true
        run_soc=@elapsed UB_soc, = getSDPUpperBound_gd_multiple_permutation(pitprops, r, k_split, usePSDs, useSOCs, useCuts)
        @show run_soc, UB_soc

        usePSDs=false
        useCuts=false
        useSOCs=false
        run_highdim=@elapsed UB_highdim, = getSDPUpperBound_gd_multiple_permutation(pitprops, r, k_split, usePSDs, useSOCs, useCuts)
        @show run_highdim, UB_highdim


        push!(results_run, ["Pitprops", n, r, sum(k_split), k_split,
            UB_psd, run_psd ,UB_cuts, run_cuts, UB_soc, run_soc, UB_highdim, run_highdim
        ])


        CSV.write("exploreUCI_bounds_table2.csv", results_run, append=true)
        results_run=similar(results_template, 0)



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
#
n=size(normwine,1)
B=sqrt(normwine);

for r=2:3
    for k_each=5:5:10
    results_run=similar(results_template, 0)
    k_split=repeat([k_each], r)

        usePSDs=true
        useCuts=false
        useSOCs=false
        run_psd=@elapsed UB_psd, = getSDPUpperBound_gd_multiple_permutation(normwine, r, k_split, usePSDs, useSOCs, useCuts)
        @show run_psd, UB_psd

        usePSDs=false
        useCuts=true
        useSOCs=true
        run_cuts=@elapsed UB_cuts, = getSDPUpperBound_gd_multiple_permutation(normwine, r, k_split, usePSDs, useSOCs, useCuts)
        @show run_cuts, UB_cuts

        usePSDs=false
        useCuts=false
        useSOCs=true
        run_soc=@elapsed UB_soc, = getSDPUpperBound_gd_multiple_permutation(normwine, r, k_split, usePSDs, useSOCs, useCuts)
        @show run_soc, UB_soc

        usePSDs=false
        useCuts=false
        useSOCs=false
        run_highdim=@elapsed UB_highdim, = getSDPUpperBound_gd_multiple_permutation(normwine, r, k_split, usePSDs, useSOCs, useCuts)
        @show run_highdim, UB_highdim


        push!(results_run, ["Wine", n, r, sum(k_split), k_split,
            UB_psd, run_psd ,UB_cuts, run_cuts, UB_soc, run_soc, UB_highdim, run_highdim
        ])


        CSV.write("exploreUCI_bounds_table2.csv", results_run, append=true)
        results_run=similar(results_template, 0)



    end
end

normionosphere=load("data/ionosphere.jld",  "normionosphere")
@show n=size(normionosphere,1)

for r=3:-1:2
    for k_each=20:5:20
    results_run=similar(results_template, 0)
    k_split=repeat([k_each], r)

            usePSDs=true
            useCuts=false
            useSOCs=false
            run_psd=@elapsed UB_psd, = getSDPUpperBound_gd_multiple_permutation(normionosphere, r, k_split, usePSDs, useSOCs, useCuts)
            @show run_psd, UB_psd

            usePSDs=false
            useCuts=true
            useSOCs=true
            run_cuts=@elapsed UB_cuts, = getSDPUpperBound_gd_multiple_permutation(normionosphere, r, k_split, usePSDs, useSOCs, useCuts)
            @show run_cuts, UB_cuts

            usePSDs=false
            useCuts=false
            useSOCs=true
            run_soc=@elapsed UB_soc, = getSDPUpperBound_gd_multiple_permutation(normionosphere, r, k_split, usePSDs, useSOCs, useCuts)
            @show run_soc, UB_soc

            usePSDs=false
            useCuts=false
            useSOCs=false
            run_highdim=@elapsed UB_highdim, = getSDPUpperBound_gd_multiple_permutation(normionosphere, r, k_split, usePSDs, useSOCs, useCuts)
            @show run_highdim, UB_highdim


            push!(results_run, ["ionosphere", n, r, sum(k_split), k_split,
                UB_psd, run_psd ,UB_cuts, run_cuts, UB_soc, run_soc, UB_highdim, run_highdim
            ])


            CSV.write("exploreUCI_bounds_table2.csv", results_run, append=true)
            results_run=similar(results_template, 0)


    end
end


# Now do miniBooNE

normminiboone=load("data/miniBoone.jld",  "normMiniBooNE")
@show n=size(normminiboone,1)

for r=2:3
    for k_each=20:5:20
    results_run=similar(results_template, 0)
    k_split=repeat([k_each], r)

        usePSDs=true
        useCuts=false
        useSOCs=false
        run_psd=@elapsed UB_psd, = getSDPUpperBound_gd_multiple_permutation(normminiboone, r, k_split, usePSDs, useSOCs, useCuts)
        @show run_psd, UB_psd

        usePSDs=false
        useCuts=false # Cuts not needed for this dataset since it's an exact relaxation without them, and ArPack sometimes throws weird errors.
        useSOCs=true
        run_cuts=@elapsed UB_cuts, = getSDPUpperBound_gd_multiple_permutation(normminiboone, r, k_split, usePSDs, useSOCs, useCuts)
        @show run_cuts, UB_cuts

        usePSDs=false
        useCuts=false
        useSOCs=true
        run_soc=@elapsed UB_soc, = getSDPUpperBound_gd_multiple_permutation(normminiboone, r, k_split, usePSDs, useSOCs, useCuts)
        @show run_soc, UB_soc

        usePSDs=false
        useCuts=false
        useSOCs=false
        run_highdim=@elapsed UB_highdim, = getSDPUpperBound_gd_multiple_permutation(normminiboone, r, k_split, usePSDs, useSOCs, useCuts)
        @show run_highdim, UB_highdim


        push!(results_run, ["miniBooNE", n, r, sum(k_split), k_split,
            UB_psd, run_psd ,UB_cuts, run_cuts, UB_soc, run_soc, UB_highdim, run_highdim
        ])


        CSV.write("exploreUCI_bounds_table2.csv", results_run, append=true)
        results_run=similar(results_template, 0)



    end
end


# Now do geographic-location-of-music
normGeographic=load("data/geography.jld",  "normGeographic")
@show n=size(normGeographic,1)

for r=2:3
    for k_each=20:5:20
    results_run=similar(results_template, 0)
    k_split=repeat([k_each], r)

            usePSDs=true
            useCuts=false
            useSOCs=false
            run_psd=@elapsed UB_psd, = getSDPUpperBound_gd_multiple_permutation(normGeographic, r, k_split, usePSDs, useSOCs, useCuts)
            @show run_psd, UB_psd

            usePSDs=false
            useCuts=true
            useSOCs=true
            run_cuts=@elapsed UB_cuts, = getSDPUpperBound_gd_multiple_permutation(normGeographic, r, k_split, usePSDs, useSOCs, useCuts)
            @show run_cuts, UB_cuts

            usePSDs=false
            useCuts=false
            useSOCs=true
            run_soc=@elapsed UB_soc, = getSDPUpperBound_gd_multiple_permutation(normGeographic, r, k_split, usePSDs, useSOCs, useCuts)
            @show run_soc, UB_soc

            usePSDs=false
            useCuts=false
            useSOCs=false
            run_highdim=@elapsed UB_highdim, = getSDPUpperBound_gd_multiple_permutation(normGeographic, r, k_split, usePSDs, useSOCs, useCuts)
            @show run_highdim, UB_highdim


            push!(results_run, ["Geographical", n, r, sum(k_split), k_split,
                UB_psd, run_psd ,UB_cuts, run_cuts, UB_soc, run_soc, UB_highdim, run_highdim
            ])


            CSV.write("exploreUCI_bounds_table2.csv", results_run, append=true)
            results_run=similar(results_template, 0)


    end
end

# Now do communities

normcommunities=load("data/communities.jld", "normCommunities")

n=size(normcommunities,1)
B=sqrt(normcommunities);

for r=2:3
    for k_each=20:5:20
    results_run=similar(results_template, 0)
    k_split=repeat([k_each], r)

        usePSDs=true
        useCuts=false
        useSOCs=false
        run_psd=-1.0#@elapsed UB_psd, = getSDPUpperBound_gd_multiple_permutation(normcommunities, r, k_split, usePSDs, useSOCs, useCuts)
        UB_psd=-1.0 # Doesn't run on laptop due to excessive memory requirements, but you can uncomment this if you have a very powerful computer
        @show run_psd, UB_psd

        usePSDs=false
        useCuts=true
        useSOCs=true
        run_cuts=@elapsed UB_cuts, = getSDPUpperBound_gd_multiple_permutation(normcommunities, r, k_split, usePSDs, useSOCs, useCuts)
        @show run_cuts, UB_cuts

        usePSDs=false
        useCuts=false
        useSOCs=true
        run_soc=@elapsed UB_soc, = getSDPUpperBound_gd_multiple_permutation(normcommunities, r, k_split, usePSDs, useSOCs, useCuts)
        @show run_soc, UB_soc

        usePSDs=false
        useCuts=false
        useSOCs=false
        run_highdim=@elapsed UB_highdim, = getSDPUpperBound_gd_multiple_permutation(normcommunities, r, k_split, usePSDs, useSOCs, useCuts)
        @show run_highdim, UB_highdim


        push!(results_run, ["Communities", n, r, sum(k_split), k_split,
            UB_psd, run_psd ,UB_cuts, run_cuts, UB_soc, run_soc, UB_highdim, run_highdim
        ])


        CSV.write("exploreUCI_bounds_table2.csv", results_run, append=true)
        results_run=similar(results_template, 0)



    end
end

# Now do arrythmia

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
    for k_each=10:5:10
    results_run=similar(results_template, 0)
    k_split=repeat([k_each], r)

        usePSDs=true
        useCuts=false
        useSOCs=false
        run_psd=@elapsed -1.0#UB_psd, = getSDPUpperBound_gd_multiple_permutation(normarrhythmia, r, k_split, usePSDs, useSOCs, useCuts)
        UB_psd=-1.0
        @show run_psd, UB_psd

        usePSDs=false
        useCuts=true
        useSOCs=true
        run_cuts=@elapsed UB_cuts, = getSDPUpperBound_gd_multiple_permutation(normarrhythmia, r, k_split, usePSDs, useSOCs, useCuts)
        @show run_cuts, UB_cuts

        usePSDs=false
        useCuts=false
        useSOCs=true
        run_soc=@elapsed UB_soc, = getSDPUpperBound_gd_multiple_permutation(normarrhythmia, r, k_split, usePSDs, useSOCs, useCuts)
        @show run_soc, UB_soc

        usePSDs=false
        useCuts=false
        useSOCs=false
        run_highdim=@elapsed UB_highdim, = getSDPUpperBound_gd_multiple_permutation(normarrhythmia, r, k_split, usePSDs, useSOCs, useCuts)
        @show run_highdim, UB_highdim


        push!(results_run, ["Arrythmia", n, r, sum(k_split), k_split,
            UB_psd, run_psd ,UB_cuts, run_cuts, UB_soc, run_soc, UB_highdim, run_highdim
        ])


        CSV.write("exploreUCI_bounds_table2.csv", results_run, append=true)
        results_run=similar(results_template, 0)



    end
end


# Finally, do microMass
micromass=load("data/microMass.jld", "normMicroMass");
B=sqrt(micromass)
B=real(B); #imaginary part of order 1e-12, due to round-off errors
n=size(micromass, 1)

for r=2:3
    for k_each=5:5:10
    results_run=similar(results_template, 0)
    k_split=repeat([k_each], r)

        usePSDs=true
        useCuts=false
        useSOCs=false
        run_psd=-1.0#@elapsed UB_psd, = getSDPUpperBound_gd_multiple_permutation(normarrhythmia, r, k_split, usePSDs, useSOCs, useCuts)
        UB_psd=-1.0 # Note: this will not run with a feasible amount of RAM
        @show run_psd, UB_psd

        usePSDs=false
        useCuts=true
        useSOCs=true
        run_cuts=-1.0#@elapsed UB_cuts, = getSDPUpperBound_gd_multiple_permutation(normarrhythmia, r, k_split, usePSDs, useSOCs, useCuts)
        UB_cuts=-1.0 # Note: this will not run within a reasonable amount of time
        @show run_cuts, UB_cuts

        usePSDs=false
        useCuts=false
        useSOCs=true
        run_soc=-1.0#@elapsed UB_soc, = getSDPUpperBound_gd_multiple_permutation(micromass, r, k_split, usePSDs, useSOCs, useCuts)
        UB_soc=-1.0 # You can uncomment this to run this one, but it requires a lot of memory
        @show run_soc, UB_soc

        usePSDs=false
        useCuts=false
        useSOCs=false
        run_highdim=@elapsed UB_highdim, = getSDPUpperBound_gd_multiple_permutation(micromass, r, k_split, usePSDs, useSOCs, useCuts)
        @show run_highdim, UB_highdim


        push!(results_run, ["MicroMass", n, r, sum(k_split), k_split,
            UB_psd, run_psd ,UB_cuts, run_cuts, UB_soc, run_soc, UB_highdim, run_highdim
        ])


        CSV.write("exploreUCI_bounds_table2.csv", results_run, append=true)
        results_run=similar(results_template, 0)



    end
end
