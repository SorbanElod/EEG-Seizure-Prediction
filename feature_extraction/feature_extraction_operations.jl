import Pkg;
Pkg.add("Statistics");
Pkg.add("FFTW");
Pkg.add("OrderedCollections");
Pkg.add("DSP") # for Hanning windowing
using Statistics
using Printf
using FFTW
using OrderedCollections
using DSP


frequency_bands = OrderedDict(
    "delta" => (0, 4),
    "theta" => (4, 8),
    "alpha" => (8, 14),
    "beta" => (15, 30),
    "low_gamma" => (30, 100),
    "high_gamma" => (100, 200)
)


function extract_statistical_features(signal::Vector{Float64})
    min_value = min(signal...)
    max_value = max(signal...)
    mean = Statistics.mean(signal)
    std = Statistics.std(signal)
    return min_value, max_value, mean, std
end


function fast_fourier_transform(signal::Vector{Float64}, sampling_frequency::Int)
    fft_signal = fft(signal)
    len_of_signal = length(signal)
    freqs = fftfreq(len_of_signal, sampling_frequency)
    half_size = div(len_of_signal, 2)
    one_sided_freqs = freqs[1:half_size]
    amplitudes = abs.(fft_signal[1:half_size])
    return one_sided_freqs, amplitudes
end



function extract_frequency_bands_from_signal(signal, sampling_frequency)
    # Applying Hanning window to the signal before the FFT 
    hanning_window = hanning(length(signal))
    windowed_signal = signal .* hanning_window

    frequencies, amplitudes = fast_fourier_transform(windowed_signal, sampling_frequency)

    num_bands = length(frequency_bands)
    band_powers = zeros(Float64, num_bands)
    for (idx, (band_name, band_range)) in enumerate(frequency_bands)
        band_start = band_range[1]
        band_end = band_range[2]
        band_indices = findall(f -> band_start <= f < band_end, frequencies)
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


function extract_statistical_and_frequency_features_from_signal(signal::TimeSeries)
    min, max, mean, std = extract_statistical_features(signal.x)
    band_powers = extract_frequency_bands_from_three_segments(signal.x, signal.fs)    
    features = [min, max, mean, std, band_powers...]
    return features
end


function extract_statistical_and_frequency_features_from_all_channels(eeg::CustomEEG)
    channels = collect(keys(eeg.signals))
    features = Vector{Float64}()
    for channel in channels
        append!(features, extract_statistical_and_frequency_features_from_signal(eeg.signals[channel]))        
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
