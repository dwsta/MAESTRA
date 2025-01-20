function [raw_signal,frames] = compute_traces_whole_well(cfg_data, rootdir,alias,whichpass)

arguments

    cfg_data (1,1) struct
    rootdir (1,1) string
    alias (1,1) string
    whichpass (1,1) double = length(cfg_data.Deformation);

end

if cfg_data.TFM.DoTFM
    pivdir = fullfile(rootdir,'output',alias,'piv');
    tfmdir = fullfile(rootdir,'output',alias,'tfm');
    try
        % Try to find the right pass
        filtfm = sprintf('traction_stresses_pass%d.bin',whichpass);
    catch
        error(sprintf('Couldnot find the traction stress for pass %d',whichpass));
        filtfm = 'traction_stresses.bin';
    end
    ref_frames = load(fullfile(pivdir,'reference_frames.txt'));
    [xvec,yvec,tvec,Tx,Ty] = contPIV.readPIV_bin(fullfile(tfmdir,filtfm));
    % Remove the measurements at reference frames, as it can introduce errors
    % later in the analysis
    vecind = load(fullfile(pivdir,'vecind.txt'));
    [tf,idx] = ismember(ref_frames,vecind); % This finds the instances we need to remove from tvec, and slices in U and V
    idx(idx ==0)=[];
    if ~isempty(idx)
        tvec(idx)=[];
        Tx(:,:,idx,:)=[];
        Ty(:,:,idx,:)=[];
    end
    frames = tvec;
    [X,Y] = meshgrid(xvec,yvec);
    raw_signal = squeeze(mean(sqrt(Tx.^2+Ty.^2),[1 2],'omitnan'));
else
    % Ngrids = length(cfg_data.Deformation);
    pivdir = fullfile(rootdir,'output',alias,'piv');
    filpiv = ['smooth_deformations_pass',num2str(whichpass),'.bin'];
    ref_frames = load(fullfile(pivdir,'reference_frames.txt'));
    [xvec,yvec,tvec,U,V] = contPIV.readPIV_bin(fullfile(pivdir,filpiv));
    % Remove the measurements at reference frames, as it can introduce errors
    % later in the analysis
    vecind = load(fullfile(pivdir,'vecind.txt'));
    [tf,idx] = ismember(ref_frames,vecind); % This finds the instances we need to remove from tvec, and slices in U and V
    idx(idx ==0)=[];
    if ~isempty(idx)
        tvec(idx)=[];
        U(:,:,idx,:)=[];
        V(:,:,idx,:)=[];
    end
    frames = tvec;
    [X,Y] = meshgrid(xvec,yvec);
    div = contPIV.divergence_rik(X,Y,U,V);
    raw_signal = squeeze(mean(sqrt(div.^2),[1 2],'omitnan'));
end
raw_signal = raw_signal';
frames = frames';
end