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


function main()
    folder_path = "/media/eszter/G/machine_learning_cs/kaggle_data/Dog_1/Dog_1/"
    file_path = get_file_path(folder_path, "Dog_1_interictal_segment_0052.mat")
    print(file_path)
    eeg = load_eeg(file_path, "interictal_segment_52");
    channels = keys(eeg.signals)
    # print(eeg.signals)
    print(channels)
    print(eeg.signals["NVC1202_32_002_Ecog_c013"])
end

main()