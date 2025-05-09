include("file_operations.jl")
include("feature_extraction_operations.jl")

function main()
    folder_path = "/media/eszter/G/machine_learning_cs/kaggle_data/Dog_1/Dog_1/"
    interictal_file_name = "Dog_1_interictal_segment_0047.mat"
    interictal_file_path = get_file_path(folder_path, interictal_file_name)
    interictal_eeg= load_eeg(interictal_file_path, get_segment_name(interictal_file_name));
    
    interictal_features = extract_statistical_and_frequency_features_from_all_channels(interictal_eeg)
    println(length(interictal_features))
    println(interictal_features)
end

main()
