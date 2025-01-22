classdef contPIV

methods (Static = true)
    
    [output] = compute_deformations(IMAGES,cfg_data,filpiv,alias,iexp,reference_frames,vecind,xdrift,ydrift);
    
    [IMAGES] = imageLoader(imgpath,ext,precision,parRun);

    [reference_frames,vecind,IMAGES] = reference_selection(IMAGES,imgpath,ext,cfg_data,filpiv,filoutVEL,iexp,alias)

    [xvec,yvec,tvec,repvec,U,V,varargout] = PIV_multipass_UQ(IMAGES,cfg_data,reference_frames,vecind,X0,Y0,T0,U0,V0,C0);

    [z,s,exitflag] = smoothn_rik(varargin);

    [DIV] = divergence_rik(X,Y,U,V);
    
    function stat = wrapperBatchPIV(rootdir,cfg_data,filoutVEL,t,filInd)
     
        if any(gpuDeviceTable().DeviceSelected) % Check for GPU
            cfg_data.Performance.UseGPU = true; 
            stat = 'GPU';
        else
            cfg_data.Performance.UseGPU = false;
            stat = 'CPU';
        end

      
        for ifile = filInd
            alias = t.Alias{ifile};
            imgPath = t.Location{ifile}; ext = 'tif';
        
            contPIV.wrapperPIV(rootdir,alias,cfg_data,filoutVEL,ifile,imgPath,ext);
        end

     end

    function wrapperPIV(rootdir,alias,cfg_data,filoutVEL,ifile,imgPath,ext)
        disp(rootdir)
        disp(alias)
        pivdir = fullfile(rootdir,'output',alias,'piv');
        if ~exist(pivdir,'dir'); mkdir(pivdir); end
        
        % Load all images to RAM
        [IMAGES] = contPIV.imageLoader(imgPath,ext);
        
        [reference_frames, vecind,IMAGES] = contPIV.reference_selection(IMAGES,imgPath,ext,cfg_data,pivdir,filoutVEL,ifile,alias);
        vecind = [1:10]; %TESTING
        xdrift = zeros(1,1,size(IMAGES,3));
        ydrift = zeros(1,1,size(IMAGES,3));            
        contPIV.compute_deformations(IMAGES,cfg_data,pivdir,alias,ifile,reference_frames,vecind,xdrift,ydrift);
        clear all;

    end

    function wrapperPIVSPMD(rootdir,alias,cfg_data,filoutVEL,ifile,IMAGES,imgPath)
        
        pivdir = fullfile(rootdir,'output',alias,'piv');
        if ~exist(pivdir,'dir'); mkdir(pivdir); end
        
        ext = 'tif';
        [reference_frames, vecind,IMAGES] = contPIV.reference_selection(IMAGES,imgPath,ext,cfg_data,pivdir,filoutVEL,ifile,alias);
        vecind = [1:10]; %TESTING
        xdrift = zeros(1,1,size(IMAGES,3));
        ydrift = zeros(1,1,size(IMAGES,3));            
        contPIV.compute_deformations(IMAGES,cfg_data,pivdir,alias,ifile,reference_frames,vecind,xdrift,ydrift);
        clear IMAGES xdrift ydrift; 

    end


    
    [X,Y,T,U,V,Xdrift,Ydrift] = readPIV_bin(filename);


end





end