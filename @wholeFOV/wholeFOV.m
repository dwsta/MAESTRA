classdef wholeFOV

methods (Static = true)


    whole_FOV_analysis(rootdir,alias,cfg_data,output_file_connection,whichPass);

    [raw_signal,frames] = compute_traces_whole_well(cfg_data, rootdir,alias,whichPass);

    [peakID,time,signal] = split_peaks(raw_time,raw_signal,varargin);

    [metrics] = measure_peaks(peakIDs,time_interp,signal_interp);
end





end