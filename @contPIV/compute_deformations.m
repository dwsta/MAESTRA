
function [output] = compute_deformations(IMAGES,cfg_data,filpiv,alias,iexp,reference_frames,vecind,xdrift,ydrift)

% Multipass PIV with UQ
Ngrids = length(cfg_data.Deformation);
output = struct;

% Send images to GPU
useGPU = cfg_data.Performance.UseGPU;

% Keep only the frames to analyze
IMAGES = IMAGES(:,:,[reference_frames vecind]);

% Push images to GPU
if useGPU
    IMAGES = gpuArray(IMAGES);
end

disp(whos("IMAGES"));
disp(any(gpuDeviceTable().DeviceSelected));

% Iterate through grids
for igrid = 1:Ngrids
    dum_cfg = cfg_data;
    dum_cfg.Deformation = dum_cfg.Deformation(igrid);
    if igrid == 1
        [X0,Y0,T0,U0,V0] = deal([]);
    else
        [X0,Y0,T0] = meshgrid(xvec,yvec,tvec); U0 = nanmean(U,4); V0 = nanmean(V,4);
    end
    
    %     updateStatus(app,iexp,{sprintf('Image quality pass %s %s',num2str(igrid),datestr(now))});
    %     writeToLog(app.Logfile,[alias,' Image quality pass ',num2str(igrid)]);

    % Compute cross-correlation correction
    %     CIQ = image_quality(IMAGES,dum_cfg,reference_frames,vecind,X0,Y0,T0,U0,V0);

    % updateStatus(app,iexp,{sprintf('PIV pass %s %s',num2str(igrid),datestr(now))});
    % writeToLog(app.Logfile,[alias,' PIV pass ',num2str(igrid)]);

    disp(any(gpuDeviceTable().DeviceSelected));
    if useGPU
        [xvec,yvec,tvec,~,U,V] = contPIV.PIV_multipass_UQ(IMAGES,dum_cfg,1:length(reference_frames),...
            [1:length(vecind)],X0,Y0,T0,U0,V0,[]);
    else 
        % N_FRAMES = length(vecind);
        % N_REFS = length(reference_frames);
        % [IM_HEIGHT,IM_WIDTH, ~] = size(IMAGES);
        % WDW_SZ = cfg_data.Deformation.wdw_size;                  % window size;
        % WDW_SPC = cfg_data.Deformation.wdw_spacing;              % window distance
        % numelementsy=floor((IM_HEIGHT-WDW_SZ)/WDW_SPC)+1;
        % numelementsx=floor((IM_WIDTH-WDW_SZ)/WDW_SPC)+1;
        % N_WDWS = numelementsx*numelementsy;
        % sx = 1:WDW_SPC:IM_WIDTH-WDW_SZ+1;
        % sy = 1:WDW_SPC:IM_HEIGHT-WDW_SZ+1;
        % xvec = squeeze(single(sx));
        % yvec = squeeze(single(sy));
        % tvec = squeeze(single(vecind));
        % 
        % % U = zeros(numelementsy,numelementsx,N_FRAMES);
        % % V = zeros(numelementsy,numelementsx,N_FRAMES);

        [xvec,yvec,tvec,~,U,V] = contPIV.PIV_multipass_UQ(IMAGES,dum_cfg,1:length(reference_frames),...
            [1:length(vecind)],X0,Y0,T0,U0,V0,[]);
        % return;
        % 
        % for iframe = 1 : N_FRAMES
        %     vecindt = iframe;
        %     F=parfeval(@contPIV.PIV_multipass_UQ,6,IMAGES(:,:,[reference_frames,vecindt]),dum_cfg,[1:length(reference_frames)],[length(reference_frames)+1],X0,Y0,T0,U0,V0,[]);
        %     % [~,~,~,~,Ut,Vt] = contPIV.PIV_multipass_UQ(IMAGES(:,:,[reference_frames,vecindt]),dum_cfg,[1:length(reference_frames)],[length(reference_frames)+1],X0,Y0,T0,U0,V0,[]);
        %     % [~,~,~,~,Ut,Vt] = contPIV.PIV_multipass_UQ(IMAGES,dum_cfg,reference_frames,vecindt,X0,Y0,T0,U0,V0,[]);
        %     % U(:,:,iframe) = Ut;
        %     % V(:,:,iframe) = Vt;
        % end
        % wait(F);
    end
          
    disp(any(gpuDeviceTable().DeviceSelected));

    [NY,NX,NF,NR] = size(U);
    fid2 = fopen(fullfile(filpiv,['deformations_pass',num2str(igrid),'.bin']),'w');
    fwrite(fid2,'v3','uchar');
    fwrite(fid2,[NX,NY,NF,NR],'single');
    fwrite(fid2,xvec + dum_cfg.Deformation.wdw_size/2,'single');
    fwrite(fid2,yvec + dum_cfg.Deformation.wdw_size/2,'single');
    fwrite(fid2,tvec,'single');
    % fwrite(fid2,repvec,'single');
    % fwrite(fid2,xdrift(vecind),'single');
    % fwrite(fid2,ydrift(vecind),'single');
    fwrite(fid2,U,'single');
    fwrite(fid2,V,'single');
    fclose(fid2);

    % PIV post-processing with smoothn
    % updateStatus(app,iexp,{sprintf('post-processing PIV pass %s %s',num2str(igrid),datestr(now))});
    % writeToLog(app.Logfile,[alias,' post-processing PIV pass ',num2str(igrid)]);
    %         [U,V] = filter_PIV(U,V);

    % if useGPU% CPU -> GPU
    %     U = gpuArray(U);
    %     V = gpuArray(V);
    % end

    U = contPIV.smoothn_rik(U,0.001,'robust');
    V = contPIV.smoothn_rik(V,0.001,'robust');

    if useGPU % GPU -> CPU
        U = gather(U);
        V = gather(V);
    end

    % disp(any(gpuDeviceTable().DeviceSelected));
    
    output(igrid).xvec = xvec + dum_cfg.Deformation.wdw_size/2;
    output(igrid).yvec = yvec + dum_cfg.Deformation.wdw_size/2;
    output(igrid).tvec = tvec;
    % output(igrid).repvec = repvec;
    output(igrid).U = U;
    output(igrid).V = V;

    [NY,NX,NF,NR] = size(U);
    fid2 = fopen(fullfile(filpiv,['smooth_deformations_pass',num2str(igrid),'.bin']),'w');
    fwrite(fid2,'v3','uchar');
    fwrite(fid2,[NX,NY,NF,NR],'single');
    fwrite(fid2,xvec + dum_cfg.Deformation.wdw_size/2,'single');
    fwrite(fid2,yvec + dum_cfg.Deformation.wdw_size/2,'single');
    fwrite(fid2,tvec,'single');
    % fwrite(fid2,repvec,'single');
    % fwrite(fid2,xdrift(vecind),'single');
    % fwrite(fid2,ydrift(vecind),'single');
    fwrite(fid2,U,'single');
    fwrite(fid2,V,'single');
    fclose(fid2);

    end
        
        
        
        
            
        
end


