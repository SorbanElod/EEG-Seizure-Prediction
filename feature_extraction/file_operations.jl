import Pkg;
Pkg.add("MAT");
Pkg.add("EEGToolkit");
Pkg.add("Plots");
using MAT
using EEGToolkit
using Plots

struct CustomEEG
    signals::OrderedDict{String, TimeSeries}
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
  
    signals = OrderedDict{String, TimeSeries}()
  
    for i in 1:size(channels, 2)
        channel_label = channels[1, i]
        signals[channel_label] = TimeSeries(eeg_data[i, :], round(sampling_freq))
    end
  
    return CustomEEG(signals)
end


function get_file_path(folder_path::String, file_name::String)
    return joinpath(folder_path, file_name)
end


function get_segment_name(file_name::String)
    file_name_without_extension = replace(file_name, ".mat" => "")
    splitted_file_name = split(file_name_without_extension, "_")
    number_without_padding = parse(Int, splitted_file_name[5])
    number_string = string(number_without_padding)
    return splitted_file_name[3] * "_segment_" * number_string
end
