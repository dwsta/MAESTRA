clearvars;
% restoredefaultpath;
% addpath 'D:\Cloud\Git\MAESTRA';

rootdir = 'D:\Data\Contractility_batch_test\contractility_run_20250121_210036\';
% rootdir = 'I:\T3_1224_CAT\contractility_run_20250119_123809\test';
% rootdir = 'E:\Contractility_3\2023-03-25\pretreatment_2023-03-25\test_batch';
% cd(rootdir);

t = readJobFile2('jobfile.csv', rootdir);

% Make sure the path in jobfile is correct
t = appUtils.checkVideoLocation(t,rootdir);

cfg_data = contUtils.loadJsonConfig(fullfile(rootdir,'config.json'));

filoutVEL = fullfile(rootdir,'vel');
if (cfg_data.ReferenceSelection.ReferenceSelectionMode ~= 0)
    if ~exist(filoutVEL,'dir'); mkdir(filoutVEL); end
end

warning('off');
% files = appUtils.searchForFilesRecursively('D:\Data\Contractility_batch_test\scan\Well__B_004\r_0005_c_0005\beads','.tif','',0);
% files = files{1:10};


p=parpool('local',3);

%%

timeStart = tic; 

spmd 

    if spmdIndex == 1
        maxNumCompThreads(8);
        ifile = 1;
        % Read images and balance load

        while 1 

            % fprintf('Received package in worker 1 from %d for file = %d \n',labID,ifile);
            % If worker idling
                imgPath = t.Location{ifile}; ext = 'tif';
                tt=tic;
                [IMAGES] = contPIV.imageLoader(imgPath,'tif');
                fprintf('Image loaded in %f s for file = %d \n',toc(tt),ifile);

                % Wait for request from compute workers
                package = spmdReceive('any'); % Receive {labindex,status}
            
                
                labID = package;

                tt=tic;
                IMAGES = gpuArray(IMAGES);
                spmdSend({IMAGES,ifile,tt},labID);
                fprintf('IMAGES sent back from worker 1 to worker %d \n',labID);
                % IMAGES = [];

                

            ifile = ifile + 1;
            if ifile > 20 
                % Exit signal
                spmdSend(777,[2:spmdSize]);
                break;
            end

            % disp(package{1});
            
        end

    else
        % GPU worker
        % Set as free

        if labindex == 2
            gpuDevice(1);
            cfg_data.Performance.UseGPU = true; 
            maxNumCompThreads(8);
        else
            gpuDevice([]);
            cfg_data.Performance.UseGPU = false; 
            maxNumCompThreads(8);
        end

        

        while 1 
            
            spmdSend(labindex,1);
            package = spmdReceive(1);

            % Exit signal
            if length(package) == 1
                if package == 777
                    break;
                end
            end
            IMAGES = package{1}; ifile = package{2}; tt=toc(package{3});
            fprintf('Images received from worker 1 on worker %d in time %f \n',labindex,tt);
            alias = t.Alias{ifile};
            imgPath = t.Location{ifile};
            tt=tic;
            contPIV.wrapperPIVSPMD(rootdir,alias,cfg_data,filoutVEL,ifile,IMAGES,imgPath);
            fprintf('PIV completed in %f s \n',toc(tt));
            IMAGES = [];
        end

 

    % elseif labindex == 3
    %     % CPU worker
    % 
    % 
    %     for ii = 1 : 10
    %         spmdSend({labindex,files{ii}},1);
    % 
    %         A = spmdReceive(1);
    % 
    %         disp(mean(    abs(A(:))));
    %     end
    end
end

timeStop = toc(timeStart)


%% Serial test

serialtimestart = tic;
for ifile = 21: 40
    alias = t.Alias{ifile};
    imgPath = t.Location{ifile};
    tt=tic;
    [IMAGES] = contPIV.imageLoader(imgPath,'tif');
    IMAGES = gpuArray(IMAGES);
    fprintf('Image loaded in %f s \n',toc(tt));

        
    tt=tic;
    contPIV.wrapperPIVSPMD(rootdir,alias,cfg_data,filoutVEL,ifile,IMAGES,imgPath);
    fprintf('PIV completed in %f s \n',toc(tt));
end
serialtimeduration = toc(serialtimestart)


return;
%% New code version


PIV = PIVUQ('config.json');

newtimestart = tic;

for ifile = 1 : 2
    alias = t.Alias{ifile};
    imgPath = t.Location{ifile};
    tt=tic;
    [IMAGES] = contPIV.imageLoader(imgPath,'tif');
    
    fprintf('Image loaded in %f s \n',toc(tt));


    pivdir = fullfile(rootdir,'output',alias,'piv');
        if ~exist(pivdir,'dir'); mkdir(pivdir); end
        
        ext = 'tif';

    [reference_frames,vecind,IMAGES] = contPIV.reference_selection(IMAGES,imgPath,ext,cfg_data,pivdir,filoutVEL,ifile,alias);
    
    PIV.IMAGES = IMAGES;
    tt=tic;
    PIV.refFrames = reference_frames;
    PIV.sessFrames = vecind;
    PIV = PIV.runPIVMPUQ;
    % contPIV.wrapperPIVSPMD(rootdir,alias,cfg_data,filoutVEL,ifile,IMAGES,imgPath);
    fprintf('PIV completed in %f s \n',toc(tt));
end
newtimeduration = toc(newtimestart)

