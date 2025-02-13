clearvars;
restoredefaultpath;

addpath 'C:\Users\MercolaLab\Documents\Git\MAESTRA';

% Read the jobfile
rootdir = "H:\T2_1224_MSP\TMRM_bluebeads_interleaved_AK_OM_20241226190508\contractility_run_20241229_190803";
t = readJobFile2('jobfile.csv', rootdir);

% Search for deformations_pass_2.bin 
parentPath = fullfile(rootdir, 'output');
[files] = appUtils.searchForFilesRecursively(parentPath, 'smooth_deformations_pass1.bin','',1);
% files(contains(files,"_001") | contains(files,"_012")) = []; % No idea what 001 and 012 wells are 
pivFolders = cellfun(@(x) fileparts(x),files,'UniformOutput',false); 
fileNames = cellfun(@(x) getFileName(x),files,'UniformOutput',false);

% Rewrite config
cfg_data =  contPIV.loadJsonConfig(fullfile(rootdir,'config.json'));

whichPass = 2;

% Read stiffness map
msp = plateReader("H:\T2_1224_MSP\MSP_stiffness_T2_1224.csv");
cfg_data.TFM.GelStifness = msp; 

warning('off');

N = length(files);
w = waitbar(0,'TFM...');
w.UserData = N;

for ifile = 1 : length(files)
    ifile
    waitbar(ifile/length(files),w);
    
    alias = fileNames{ifile};
    % Read the images
    
    % TFM_step(rootdir,alias,cfg_data);
    % f(ifile) = parfeval(backgroundPool,@contTFM.TFM_step_nPass,0,rootdir,alias,cfg_data,whichPass);
    
    contTFM.TFM_step_nPass(rootdir,alias,cfg_data,whichPass);
    waitbar(ifile/N,w)
    % Uncomment for PIV
    % disp([alias,' Loading Images']);
    % 
    % indIMG = find(contains(t.Alias,alias));
    % if length(indIMG) > 1
    %     keyboard;
    % end
    % 
    % imgpath = t.Location{indIMG};
    % ext = 'tif';
    % 
    % [IMAGES, IMAGE_NAMES] = imageLoader(imgpath,ext);
    % 
    % % Read ref_frames
    % reference_frames = load(fullfile(pivFolders{ifile},'reference_frames.txt'));
    % vecind = load(fullfile(pivFolders{ifile},'vecind.txt'));
    % 
    % [output] = compute_deformations(IMAGES,cfg_data,pivFolders{ifile},[],[],reference_frames,vecind,[],[],[3],files{ifile});
    % 

end

% afterEach(f,@(~)updateWaitbar(w),0);
% afterAll(f,@(~)delete(w),0);


function updateWaitbar(w)
    % Update a waitbar using the UserData property.

    % Check if the waitbar is a reference to a deleted object
    if isvalid(w)
        % Increment the number of completed iterations 
        w.UserData(1) = w.UserData(1) + 1;

        % Calculate the progress
        progress = w.UserData(1) / w.UserData(2);

        % Update the waitbar
        waitbar(progress,w);

        % disp(sprintf('%0.3f completed',progress));
    end
end

function [aa] = getFileName(x)
[~,aa] = fileparts(fileparts(fileparts(x)));
end

function [output] = compute_deformations(IMAGES,cfg_data,filpiv,alias,iexp,reference_frames,vecind,xdrift,ydrift,runWhichPass,pivPrevRes)

% Multipass PIV with UQ
Ngrids = runWhichPass;
output = struct;

for igrid = Ngrids
    dum_cfg = cfg_data;
    dum_cfg.Deformation = dum_cfg.Deformation(igrid);
    if igrid == 1
        [X0,Y0,T0,U0,V0] = deal([]);
    else
        % [X0,Y0,T0] = meshgrid(xvec,yvec,tvec); U0 = nanmean(U,4); V0 = nanmean(V,4);
        [X0,Y0,T0,U0,V0,xdrift,ydrift] = readPIV_bin(pivPrevRes);
        [X0,Y0,T0] = meshgrid(X0,Y0,T0); 
    end
    %     updateStatus(app,iexp,{sprintf('Image quality pass %s %s',num2str(igrid),datestr(now))});
    %     writeToLog(app.Logfile,[alias,' Image quality pass ',num2str(igrid)]);

    % Compute cross-correlation correction
    %     CIQ = image_quality(IMAGES,dum_cfg,reference_frames,vecind,X0,Y0,T0,U0,V0);

    [xvec,yvec,tvec,repvec,U,V] = PIV_multipass_UQ(IMAGES,dum_cfg,reference_frames,vecind,X0,Y0,T0,U0,V0,[]);
    
    [NY,NX,NF,NR] = size(U);
    fid2 = fopen(fullfile(filpiv,['deformations_pass',num2str(igrid),'.bin']),'w');
    fwrite(fid2,'v3','uchar');
    fwrite(fid2,[NX,NY,NF,NR],'single');
    fwrite(fid2,xvec + dum_cfg.Deformation.wdw_size/2,'single');
    fwrite(fid2,yvec + dum_cfg.Deformation.wdw_size/2,'single');
    fwrite(fid2,tvec,'single');
    fwrite(fid2,repvec,'single');
    % fwrite(fid2,xdrift(vecind),'single');
    % fwrite(fid2,ydrift(vecind),'single');
    fwrite(fid2,U,'single');
    fwrite(fid2,V,'single');
    fclose(fid2);

    % PIV post-processing with smoothn

    if cfg_data.Performance.UseGPU % CPU -> GPU
        U = gpuArray(U);
        V = gpuArray(V);
    end

    U = smoothn_rik(U,0.001,'robust');
    V = smoothn_rik(V,0.001,'robust');

    if cfg_data.Performance.UseGPU % GPU -> CPU
        U = gather(U);
        V = gather(V);
    end

    output(igrid).xvec = xvec + dum_cfg.Deformation.wdw_size/2;
    output(igrid).yvec = yvec + dum_cfg.Deformation.wdw_size/2;
    output(igrid).tvec = tvec;
    output(igrid).repvec = repvec;
    output(igrid).U = U;
    output(igrid).V = V;

    [NY,NX,NF,NR] = size(U);
    fid2 = fopen(fullfile(filpiv,['smooth_deformations_pass',num2str(igrid),'.bin']),'w');
    fwrite(fid2,'v3','uchar');
    fwrite(fid2,[NX,NY,NF,NR],'single');
    fwrite(fid2,xvec + dum_cfg.Deformation.wdw_size/2,'single');
    fwrite(fid2,yvec + dum_cfg.Deformation.wdw_size/2,'single');
    fwrite(fid2,tvec,'single');
    fwrite(fid2,repvec,'single');
    % fwrite(fid2,xdrift(vecind),'single');
    % fwrite(fid2,ydrift(vecind),'single');
    fwrite(fid2,U,'single');
    fwrite(fid2,V,'single');
    fclose(fid2);

end
end