import Pkg;
Pkg.add("MAT");
Pkg.add("EEGToolkit");
Pkg.add("Plots");
Pkg.add("Statistics");
Pkg.add("FFTW");
Pkg.add("OrderedCollections");
Pkg.add("DSP") # for Hanning windowing
using MAT
using EEGToolkit
using Plots
using Statistics
using Printf
using FFTW
using OrderedCollections
using DSP


struct CustomEEG
    signals::OrderedDict{String, TimeSeries}
end

const frequency_bands = OrderedDict(
    "delta" => (0, 4),
    "theta" => (4, 8),
    "alpha" => (8, 14),
    "beta" => (15, 30),
    "low_gamma" => (30, 100),
    "high_gamma" => (100, 200)
)

number_of_features_per_channel = 10

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


function extract_statistical_features(signal::Vector{Float64})
    min_value = min(signal...)
    max_value = max(signal...)
    mean = Statistics.mean(signal)
    std = Statistics.std(signal)
    return min_value, max_value, mean, std
end


function extract_frequency_bands_from_signal(signal, sampling_frequency)
    # Applying Hanning window to the signal before the FFT 
    hanning_window = hanning(length(signal))
    windowed_signal = signal .* hanning_window

    fft_signal = fft(windowed_signal)
    len_of_signal = length(windowed_signal)
    n = 0:len_of_signal
    T = len_of_signal/sampling_frequency
    freqs = n ./ T
    half_size = div(len_of_signal,2)
    one_sided_freqs = freqs[1:half_size]
    amplitudes = abs.(fft_signal[1:half_size])

    num_bands = length(frequency_bands)
    band_powers = zeros(Float64, num_bands)
    for (idx, (band_name, band_range)) in enumerate(frequency_bands)
        band_start = band_range[1]
        band_end = band_range[2]
        band_indices= findall(f -> band_start <= f < band_end, one_sided_freqs)
        band_power = mean(amplitudes[band_indices].^2)
        logarithmic_band_power = log10(band_power)
        band_powers[idx] = logarithmic_band_power
    end
    return band_powers
end


function extract_frequency_bands_from_three_segments(signal, sampling_frequency)
    length_of_signal = length(signal)
    mid_point = div(length_of_signal, 2)
    length_of_segment_sec = 0.5
    values_in_segment = map(Int, length_of_segment_sec * sampling_frequency)
    half_seg = div(values_in_segment, 2)
    segment_from_start = signal[1:values_in_segment]
    segment_from_middle = signal[(mid_point - half_seg + 1):(mid_point + half_seg + 1)]
    segment_from_end = signal[(length_of_signal-values_in_segment + 1):length_of_signal]

    frequency_bands = Vector{Float64}()
    append!(frequency_bands, extract_frequency_bands_from_signal(segment_from_start, sampling_frequency))
    append!(frequency_bands, extract_frequency_bands_from_signal(segment_from_middle, sampling_frequency))
    append!(frequency_bands, extract_frequency_bands_from_signal(segment_from_end, sampling_frequency))
    return frequency_bands
end



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

    print_features("Interictal segment", min_of_first, max_of_first, mean_of_first, std_of_first, band_powers_1)
    print_features("Preictal segment", min_of_second, max_of_second, mean_of_second, std_of_second, band_powers_2)
end


function extract_features_from_signal(signal::TimeSeries)
    min, max, mean, std = extract_statistical_features(signal.x)
    band_powers = extract_frequency_bands_from_three_segments(signal.x, signal.fs)    
    features = [min, max, mean, std, band_powers...]
    return features
end


function extract_statistical_and_frequency_features_from_all_channels(eeg::CustomEEG)
    channels = collect(keys(eeg.signals))
    features = Vector{Float64}()
    for channel in channels
        append!(features, extract_features_from_signal(eeg.signals[channel]))        
    end
    return features
end


function extract_statistical_features_from_all_channels(eeg::CustomEEG)
    channels = collect(keys(eeg.signals))
    features = Vector{Float64}()
    for channel in channels
        append!(features, extract_statistical_features(eeg.signals[channel].x))        
    end
    return features
end


function extract_frequency_features_from_all_channels(eeg::CustomEEG)
    channels = collect(keys(eeg.signals))
    features = Vector{Float64}()
    for channel in channels
        append!(features, extract_frequency_bands_from_three_segments(eeg.signals[channel].x, eeg.signals[channel].fs))        
    end
    return features
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
    interictal_features = extract_statistical_and_frequency_features_from_all_channels(interictal_eeg)
    preictal_features = extract_statistical_and_frequency_features_from_all_channels(preictal_eeg)
    println(preictal_features)
    # print(channels)
    # channels = keys(eeg.signals)
    # print(eeg.signals["NVC1202_32_002_Ecog_c013"])
end

main()