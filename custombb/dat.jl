using DataFrames, DelimitedFiles

mutable struct problem
    data::Array{Float64}
    Sigma::Array{Float64}
end


################ Pitprops
Dir = "../Data/"
Sigma = DelimitedFiles.readcsv(string(Dir,"pitprops.csv"),header=false);
pitprops = problem(LinearAlgebra.chol(Sigma),  Sigma);
