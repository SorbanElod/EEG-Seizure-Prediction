include("feature_extraction_operations.jl")


function print_features(signal_identifier, min, max, mean, std, band_powers)
    @printf("%s: min = %s, max = %s, mean = %s, std = %s \n", signal_identifier, min, max, mean, std)
    for (idx, (band_name, band_range)) in enumerate(frequency_bands)
        band_start = band_range[1]
        band_end = band_range[2]
        println("\t\t$band_name ([$band_start - $band_end]): $(band_powers[idx])")
    end
end


function compare_features_of_two_signal(signal_1::TimeSeries, signal_2::TimeSeries)
    min_of_first, max_of_first, mean_of_first, std_of_first = extract_statistical_features(signal_1.x)
    band_powers_1 = extract_frequency_bands_from_three_segments(signal_1.x, signal_1.fs)

    min_of_second, max_of_second, mean_of_second, std_of_second = extract_statistical_features(signal_2.x)
    band_powers_2 = extract_frequency_bands_from_three_segments(signal_2.x, signal_2.fs)

    print_features("Interictal segment", min_of_first, max_of_first, mean_of_first, std_of_first, band_powers_1[1:6])
    print_features("Preictal segment", min_of_second, max_of_second, mean_of_second, std_of_second, band_powers_2[1:6])
end


function main()
    folder_path = "/media/eszter/G/machine_learning_cs/kaggle_data/Dog_1/Dog_1/"
    interictal_file_name = "Dog_1_interictal_segment_0047.mat"
    interictal_file_path = get_file_path(folder_path, interictal_file_name)
    interictal_eeg= load_eeg(interictal_file_path, get_segment_name(interictal_file_name));

    preictal_file_name = "Dog_1_preictal_segment_0014.mat"
    preictal_file_path = get_file_path(folder_path, preictal_file_name)
    preictal_eeg= load_eeg(preictal_file_path, get_segment_name(preictal_file_name));

    channel = "NVC1202_32_002_Ecog_c013"

    # print(interictal_eeg.signals[channel].fs) - sampling frequency
    # print(interictal_eeg.signals[channel].x) - signal

    compare_features_of_two_signal(interictal_eeg.signals[channel], preictal_eeg.signals[channel])
    
end

main()