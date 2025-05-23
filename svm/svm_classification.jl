# import Pkg
# Pkg.add("LIBSVM")
# Pkg.add("Glob")
# Pkg.add("Random")
# Pkg.add("Statistics")
# Pkg.add("Printf")
using LIBSVM
using Glob
using Random
using Statistics
using Printf

include("../feature_extraction/file_operations.jl")
include("../feature_extraction/feature_extraction_operations.jl")

function load_features_for_class(folder_path::String, file_pattern::String, label::Int)
    all_feature_vectors = Vector{Float64}[]
    filenames = glob(file_pattern, folder_path)

    if isempty(filenames)
        @warn "No files found matching the pattern: $file_pattern"
        return Matrix{Float64}(undef, 0, 0), Vector{Int}()
    end

    num_features_expected = -1
    for (idx, full_file_path) in enumerate(filenames)
        file_name = basename(full_file_path)
        segment_name = get_segment_name(file_name)
        eeg = load_eeg(full_file_path, segment_name)
        features = extract_statistical_and_frequency_features_from_all_channels(eeg)

        if idx == 1
            num_features_expected = length(features)
        elseif length(features) != num_features_expected
            @warn "Feature vector length mismatch in file: $file_name"
            continue
        end
        push!(all_feature_vectors, features)
    end
    feature_matrix = hcat(all_feature_vectors...)
    labels = fill(label, size(feature_matrix, 2))
    return feature_matrix, labels
end

function main_svm(X_interictal, Y_interictal, X_preictal, Y_preictal)
    @assert size(X_interictal, 1) == size(X_preictal, 1)
    @assert size(X_interictal, 2) > 0 && size(X_preictal, 2) > 0

    X = hcat(X_interictal, X_preictal)
    Y = vcat(Y_interictal, Y_preictal)

    Random.seed!(42)
    indicies = shuffle(1:size(X, 2))
    X = X[:, indicies]
    Y = Y[indicies]

    num_train = floor(Int, 0.8 * size(X, 2))
    X_train = X[:, 1:num_train]
    Y_train = Y[1:num_train]
    X_test = X[:, (num_train + 1):end]
    Y_test = Y[(num_train + 1):end]

    train_means = mean(X_train, dims=2)
    train_stds = std(X_train, dims=2)
    train_stds_stable = train_stds .+ 1e-12
    X_train = (X_train .- train_means) ./ train_stds_stable   
    X_test = (X_test .- train_means) ./ train_stds_stable

    println("Training samples: $(size(X_train, 2)), Test samples: $(size(X_test, 2))")
    println("Training class distribution: $(count(==(0), Y_train)) interictal, $(count(==(1), Y_train)) preictal")
    println("Test class distribution: $(count(==(0), Y_test)) interictal, $(count(==(1), Y_test)) preictal") 

    # ---------------------------------------

    w = Dict{Float64, Float64}()
    w[0.0] = length(Y_train) / (count(==(0), Y_train))
    w[1.0] = length(Y_train) / (count(==(1), Y_train))
    @show w

    #ma, mc, mg = 0.0, 0.0, 0.0
    #for c in [10.0^p for p in -2:1:8]
    #   for g in [10.0^p for p in -10:1:0]
        
        model = svmtrain(X_train, Float64.(Y_train); 
        svmtype=SVC, 
        kernel=Kernel.RadialBasis,
        cost=0.1, #c,
        weights=w,
        gamma=0.001 #g
        )
        
        # ---------------------------------------
        y_pred_test, decision_values = svmpredict(model, X_test)
        accuracy = mean(y_pred_test .== Y_test)
        @printf "Test accuracy: %.2f%%\n" (accuracy * 100)
        
        tp = sum((y_pred_test .== 1) .& (Y_test .== 1))
        tn = sum((y_pred_test .== 0) .& (Y_test .== 0))
        fp = sum((y_pred_test .== 1) .& (Y_test .== 0))
        fn = sum((y_pred_test .== 0) .& (Y_test .== 1))
        
        @printf "True Positives: %d\n" tp
        @printf "True Negatives: %d\n" tn
        @printf "False Positives: %d\n" fp
        @printf "False Negatives: %d\n\n" fn

        #if accuracy > ma
        #    ma = accuracy
        #    mc = c
        #    mg = g
        #end
       #end
    #end
    #@show mg, mc, ma
end

folder_paths = ["/data/kaggle_data/Dog_1", "/data/kaggle_data/Dog_2"]
interictal_pattern = "*interictal_segment_*.mat"
preictal_pattern = "*preictal_segment_*.mat"

X_is = []
Y_is = []
X_ps = []
Y_ps = []

for folder_path in folder_paths
    X_i, Y_i = load_features_for_class(folder_path, interictal_pattern, 0)
    X_p, Y_p = load_features_for_class(folder_path, preictal_pattern, 1)

    push!(X_is, X_i)
    push!(Y_is, Y_i)
    push!(X_ps, X_p)
    push!(Y_ps, Y_p)
end

X_interictal = hcat(X_is...)
Y_interictal = vcat(Y_is...)
X_preictal = hcat(X_ps...)
Y_preictal = vcat(Y_ps...)

main_svm(X_interictal[:, 1:200], Y_interictal[1:200], X_preictal, Y_preictal)

#=
Example output for predictions after training on Dog 1 and Dog 2:

Training samples: 212, Test samples: 54
Training class distribution: 163 interictal, 49 preictal
Test class distribution: 37 interictal, 17 preictal
w = Dict(0.0 => 1.3006134969325154, 1.0 => 4.326530612244898)
Test accuracy: 88.89%
True Positives: 11
True Negatives: 37
False Positives: 0
False Negatives: 6
=#