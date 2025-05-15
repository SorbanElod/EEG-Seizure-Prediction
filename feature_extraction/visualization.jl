function plot_single_eeg(eeg::CustomEEG, s::Integer, e::Integer,
    channel::String, spacing::AbstractFloat = 1.5)

    if channel ∉ keys(eeg.signals)
    throw(ArgumentError("Attempting to plot a non-existent channel."))
    end

    p = plot(ylabel="Amplitude (µV)", xlabel="Time(s)", legend=false, ylimits=(-2, 2))

    signal = eeg.signals[channel]
    segment_start = s * signal.fs + 1
    segment_end = e * signal.fs
    if segment_start > length(signal.x) || segment_end > length(signal.x)
    throw(ArgumentError("Specified epoch range is out of bounds for signal length."))
    end
    segment = signal.x[segment_start:segment_end]
    t = range(start=s, stop=e, length=length(segment))
    normalized_signal = (segment .- mean(segment)) ./ (spacing * 2 * std(segment))

    plot!(p, t,normalized_signal, title=channel, color="blue")
    return p
end


function plot_custom_eeg(eeg::CustomEEG, s::Integer, e::Integer;
    channels::Vector{String} = String[],
    spacing::AbstractFloat = 1.5)

    if !isempty(channels) && any(ch -> ch ∉ keys(eeg.signals), channels)
    throw(ArgumentError("Attempting to plot a non-existent channel."))
    end

    signal_dict = isempty(channels) ? eeg.signals : Dict(k => v for (k, v) in eeg.signals if k in channels)
    L = length(signal_dict)

    p = plot(ylabel="Amplitude (µV)", xlabel="Time(s)", yticks=(1:L, collect(keys(signal_dict))), legend=false)

    for (i, (label, signal)) in enumerate(signal_dict)
    segment_start = s * signal.fs + 1
    segment_end = e * signal.fs
    if segment_start > length(signal.x) || segment_end > length(signal.x)
    throw(ArgumentError("Specified epoch range is out of bounds for signal length."))
    end
    segment = signal.x[segment_start:segment_end]

    t = range(start=s, stop=e, length=length(segment))

    normalized_signal = (segment .- mean(segment)) ./ (spacing * (2 * L) * std(segment))

    plot!(p, t, i .+ normalized_signal, label=label, color="blue")
    end

    return p
end



function plot_and_save_signal(signal::Vector{Float64}, filename::String, fs::Int, spacing::AbstractFloat = 1.5)
    p = plot(ylabel="Amplitude (µV)", xlabel="Time(s)", legend=false, ylimits=(-2, 2))

    t = (1:length(signal)) / fs
    normalized_signal = (signal .- mean(signal)) ./ (spacing * 2 * std(signal))

    plot!(p, t, normalized_signal, color="blue")
    savefig(p, filename)
    return p
end