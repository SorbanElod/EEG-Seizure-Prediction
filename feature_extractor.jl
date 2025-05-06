import Pkg;
Pkg.add("MAT");
Pkg.add("EEGToolkit");
Pkg.add("Plots");
Pkg.add("Statistics");
using MAT
using EEGToolkit
using Plots
using Statistics
using Printf

struct CustomEEG
    signals::Dict{String, TimeSeries}
end


function load_eeg(file::String, segment_name::String)
    file = matopen(file)
    segment = read(file, segment_name)
    close(file)
  
    eeg_data         = segment["data"]
    data_length_sec  = segment["data_length_sec"]
    sampling_freq    = segment["sampling_frequency"]
    channels         = segment["channels"]
    sequence         = segment["sequence"]
  
    signals = Dict{String, TimeSeries}()
  
    for i in 1:size(channels, 2)
        channel_label = channels[1, i]
        signals[channel_label] = TimeSeries(eeg_data[i, :], round(sampling_freq))
    end
  
    return CustomEEG(signals)
end


function get_file_path(folder_path::String, file_name::String)
    return folder_path * "/" * file_name
end


function get_segment_name(file_name::String)
    file_name_without_extension = replace(file_name, ".mat" => "")
    splitted_file_name = split(file_name_without_extension, "_")
    number_without_padding = parse(Int, splitted_file_name[5])
    number_string = string(number_without_padding)
    return splitted_file_name[3] * "_segment_" * number_string
end



function get_features_of_a_signal(signal::Vector{Float64})
    min_value = min(signal...)
    max_value = max(signal...)
    mean = Statistics.mean(signal)
    std = Statistics.std(signal)
    return min_value, max_value, mean, std
end


function print_features(signal_identifier, min, max, mean, std)
    @printf("%s: min = %s, max = %s, mean = %s, std = %s \n", signal_identifier, min, max, mean, std)
end


function compare_two_signal_metadata(signal_1::Vector{Float64}, signal_2::Vector{Float64})
    min_of_first, max_of_first, mean_of_first, std_of_first = get_features_of_a_signal(signal_1)

    min_of_second, max_of_second, mean_of_second, std_of_second = get_features_of_a_signal(signal_2)

    print_features("Interictal segment", min_of_first, max_of_first, mean_of_first, std_of_first)
    print_features("Preictal segment", min_of_second, max_of_second, mean_of_second, std_of_second)
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
    compare_two_signal_metadata(interictal_eeg.signals[channel].x, preictal_eeg.signals[channel].x)

    # print(channels)
    # channels = keys(eeg.signals)
    # print(eeg.signals["NVC1202_32_002_Ecog_c013"])
end

main()