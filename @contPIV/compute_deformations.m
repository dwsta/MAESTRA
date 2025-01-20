
function [output] = compute_deformations(IMAGES,cfg_data,filpiv,alias,iexp,reference_frames,vecind,xdrift,ydrift)

% Multipass PIV with UQ
Ngrids = length(cfg_data.Deformation);
output = struct;
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
    [xvec,yvec,tvec,repvec,U,V] = contPIV.PIV_multipass_UQ(IMAGES,dum_cfg,reference_frames,vecind,X0,Y0,T0,U0,V0,[]);
    
    [NY,NX,NF,NR] = size(U);
    fid2 = fopen(fullfile(filpiv,['deformations_pass',num2str(igrid),'.bin']),'w');
    fwrite(fid2,'v3','uchar');
    fwrite(fid2,[NX,NY,NF,NR],'single');
    fwrite(fid2,xvec + dum_cfg.Deformation.wdw_size/2,'single');
    fwrite(fid2,yvec + dum_cfg.Deformation.wdw_size/2,'single');
    fwrite(fid2,tvec,'single');
    fwrite(fid2,repvec,'single');
    fwrite(fid2,xdrift(vecind),'single');
    fwrite(fid2,ydrift(vecind),'single');
    fwrite(fid2,U,'single');
    fwrite(fid2,V,'single');
    fclose(fid2);

    % PIV post-processing with smoothn
    % updateStatus(app,iexp,{sprintf('post-processing PIV pass %s %s',num2str(igrid),datestr(now))});
    % writeToLog(app.Logfile,[alias,' post-processing PIV pass ',num2str(igrid)]);
    %         [U,V] = filter_PIV(U,V);

    if cfg_data.Performance.UseGPU % CPU -> GPU
        U = gpuArray(U);
        V = gpuArray(V);
    end

    U = contPIV.smoothn_rik(U,0.001,'robust');
    V = contPIV.smoothn_rik(V,0.001,'robust');

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
    fwrite(fid2,xdrift(vecind),'single');
    fwrite(fid2,ydrift(vecind),'single');
    fwrite(fid2,U,'single');
    fwrite(fid2,V,'single');
    fclose(fid2);

end
end