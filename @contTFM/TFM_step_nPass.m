function  TFM_step_nPass(rootdir,alias,cfg_data,nPass)

% Load PIV
nMax = length(cfg_data.Deformation);
if nPass > nMax
    error(sprintf('Pass number should be lower than %d',nMax));
end

Ngrids = nPass; 
pivdir = fullfile(rootdir,'output',alias,'piv');
PIV_file = fullfile(pivdir,['smooth_deformations_pass',num2str(Ngrids),'.bin']);
E = cfg_data.TFM.GelStiffness*1e3;
h = cfg_data.TFM.GelThickness;
fcalx = cfg_data.TFM.CalibrationFactor;
[xvec,yvec,tvec,U,V,Xdrift,Ydrift] = readPIV_bin(PIV_file);
Um = nanmean(U,4);
Vm = nanmean(V,4);
Ntimepoints = length(tvec);
nx = length(xvec);
ny = length(yvec);
nx0 = nx + mod(nx,2);
ny0 = ny + mod(ny,2);

[Xgrid,Ygrid] = meshgrid(xvec,yvec);
Nyq = 2;
TX = nan(ny0,nx0,Ntimepoints);
TY = TX;
SXX = TX; SYY = SXX; SXY = SXX;

parfor itimepoint = 1 : Ntimepoints
    u = Um(:,:,itimepoint);
    v = Vm(:,:,itimepoint);
    [tx,ty,nxx,nxy,nyy,ut0,vt0,X0,Y0,Lx,Ly]...
        = mytrac_msm_2D_rik_from_dp(E,h,h,nx0,ny0,Xgrid,Ygrid,u,v,fcalx,Nyq);
    TX(:,:,itimepoint) = tx;
    TY(:,:,itimepoint) = ty;
    SXX(:,:,itimepoint) = nxx; 
    SYY(:,:,itimepoint) = nyy;
    SXY(:,:,itimepoint) = nxy;

end

TFM_xvec = linspace(xvec(1),xvec(end),nx0);
TFM_yvec = linspace(yvec(1),yvec(end),ny0);
TFM_tvec = tvec;
TFM_repvec = 1;
[NY,NX,NF,NR] = size(TX);

% Output
TFM_dir = fullfile(rootdir,'output',alias,'tfm');
if ~exist(TFM_dir,'dir'); mkdir(TFM_dir); end

fid2 = fopen(fullfile(TFM_dir,sprintf('traction_stresses_pass%d.bin',nPass)),'w');
fwrite(fid2,'v3','uchar');
fwrite(fid2,single([NX,NY,NF,NR]),'single');
fwrite(fid2,single(TFM_xvec),'single');
fwrite(fid2,single(TFM_yvec),'single');
fwrite(fid2,single(TFM_tvec),'single');
fwrite(fid2,single(TFM_repvec),'single');
fwrite(fid2,single(Xdrift),'single');
fwrite(fid2,single(Ydrift),'single');
fwrite(fid2,single(TX),'single');
fwrite(fid2,single(TY),'single');
fclose(fid2);

fid2 = fopen(fullfile(TFM_dir,sprintf('monolayer_stresses_pass%d.bin',nPass)),'w');
fwrite(fid2,'v3','uchar');
fwrite(fid2,single([NX,NY,NF,NR]),'single');
fwrite(fid2,single(TFM_xvec),'single');
fwrite(fid2,single(TFM_yvec),'single');
fwrite(fid2,single(TFM_tvec),'single');
fwrite(fid2,single(TFM_repvec),'single');
fwrite(fid2,single(Xdrift),'single');
fwrite(fid2,single(Ydrift),'single');
fwrite(fid2,single(SXX),'single');
fwrite(fid2,single(SYY),'single');
fwrite(fid2,single(SXY),'single');
fclose(fid2);
