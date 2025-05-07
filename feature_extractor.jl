import Pkg;
Pkg.add("MAT");
Pkg.add("EEGToolkit");
Pkg.add("Plots");
Pkg.add("Statistics");
Pkg.add("FFTW");
using MAT
using EEGToolkit
using Plots
using Statistics
using Printf
using FFTW


struct CustomEEG
    signals::Dict{String, TimeSeries}
end


frequency_bands = Dict(
    "delta" => (0, 4),
    "theta" => (4, 8),
    "alpha" => (8, 14),
    "beta" => (15, 30),
    "low_gamma" => (30, 100),
    "high_gamma" => (100, 200)
)


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


function print_features(signal_identifier, min, max, mean, std, band_powers)
    @printf("%s: min = %s, max = %s, mean = %s, std = %s \n", signal_identifier, min, max, mean, std)
    for (idx, (band_name, band_range)) in enumerate(frequency_bands)
        band_start = band_range[1]
        band_end = band_range[2]
        println("\t\t$band_name ([$band_start - $band_end]): $(band_powers[idx])")
    end
end


function compare_two_signal_metadata(signal_1::TimeSeries, signal_2::TimeSeries)
    min_of_first, max_of_first, mean_of_first, std_of_first = get_features_of_a_signal(signal_1.x)
    band_powers_1 = get_frequency_bands(signal_1.x, signal_1.fs)

    min_of_second, max_of_second, mean_of_second, std_of_second = get_features_of_a_signal(signal_2.x)
    band_powers_2 = get_frequency_bands(signal_2.x, signal_2.fs)


    print_features("Interictal segment", min_of_first, max_of_first, mean_of_first, std_of_first, band_powers_1)
    print_features("Preictal segment", min_of_second, max_of_second, mean_of_second, std_of_second, band_powers_2)
end


function get_frequency_bands(signal, sampling_frequency)
    fft_signal = fft(signal)
    len_of_signal = length(signal)
    n = 1:len_of_signal
    T = len_of_signal/sampling_frequency
    freqs = n / T
    half_size = div(len_of_signal,2)
    one_sided_freqs = freqs[1:half_size]
    amplitudes = abs.(fft_signal[1:half_size])

    num_bands = length(frequency_bands)
    band_powers = zeros(Float64, num_bands)
    for (idx, (band_name, band_range)) in enumerate(frequency_bands)
        band_start = band_range[1]
        band_end = band_range[2]
        band_indices= findall(f -> band_start <= f < band_end, one_sided_freqs)
        band_power = mean(amplitudes[band_indices])^2
        band_powers[idx] = band_power
    end
    return band_powers
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
    compare_two_signal_metadata(interictal_eeg.signals[channel], preictal_eeg.signals[channel])

    # print(channels)
    # channels = keys(eeg.signals)
    # print(eeg.signals["NVC1202_32_002_Ecog_c013"])
end

main()