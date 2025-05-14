using DSP

include("visualization.jl")


function main()
    folder_path = "/media/eszter/G/machine_learning_cs/kaggle_data/Dog_1/Dog_1/"
    interictal_file_name = "Dog_1_interictal_segment_0047.mat"
    interictal_file_path = get_file_path(folder_path, interictal_file_name)
    interictal_eeg= load_eeg(interictal_file_path, get_segment_name(interictal_file_name));
    channel = "NVC1202_32_002_Ecog_c013"
    
    signal = interictal_eeg.signals[channel]
    end_idx = signal.fs 
    segment = signal.x[1:end_idx]
    hanning_window = hanning(length(segment))
    windowed_segment = segment .* hanning_window

    # display(plot_and_save_signal(windowed_segment, "signal_with_hanning.png", round(signal.fs)))
    display(plot_and_save_signal(segment, "signal_without_hanning.png", round(signal.fs)))

end

main()