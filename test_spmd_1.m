clearvars;
% restoredefaultpath;
% addpath 'D:\Cloud\Git\MAESTRA';

rootdir = 'D:\Data\Contractility_batch_test\contractility_run_20250121_210036\';
% rootdir = 'I:\T3_1224_CAT\contractility_run_20250119_123809\test';
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


p=parpool('local',2);

%%

clearAllMemoizedCaches

timeStart = tic; 


spmd 

    % if labindex == 1
    % 
    %     ifile = 1;
    %     % Read images and balance load
    % 
    %     while 1 
    %         package = spmdReceive('any'); % Receive {labindex,status}
    % 
    %         if package == 777
    %             break;
    %         end
    % 
    %         % fprintf('%s',fpath);
    %         labID = package;
    %         fprintf('Received package in worker 1 from %d for file = %d \n',labID,ifile);
    %         % If worker idling
    %             imgPath = t.Location{ifile}; ext = 'tif';
    %             tt=tic;
    %             [IMAGES] = contPIV.imageLoader(imgPath,'tif');
    %             fprintf('Image loaded in %f s \n',toc(tt));
    %             tt=tic;
    %             spmdSend({IMAGES,ifile,tt},labID);
    %             fprintf('IMAGES sent back from worker 1 to worker %d \n',labID);
    %             IMAGES = [];
    % 
    %         ifile = ifile + 1;
    %         if ifile > 10 
    %             break;
    %         end
    % 
    %         % disp(package{1});
    % 
    %     end
    % 
    % else
    %     % GPU worker
    %     % Set as free

        if labindex == 1
            gpuDevice(1);
            cfg_data.Performance.UseGPU = true; 
            
            for ifile = 1 : 10
                imgPath = t.Location{ifile}; ext = 'tif';
                tt=tic;
                [IMAGES] = contPIV.imageLoader_RS(imgPath,'tif');
                fprintf('Image loaded in %f s \n',toc(tt));
                tt=tic;
                contPIV.wrapperPIVSPMD(rootdir,alias,cfg_data,filoutVEL,ifile,IMAGES,imgPath);
                fprintf('PIV completed in %f s \n',toc(tt));
                IMAGES = [];
                
            end



        else
            gpuDevice([]);
            cfg_data.Performance.UseGPU = false; 

            for ifile = 11:20
                imgPath = t.Location{ifile}; ext = 'tif';
                tt=tic;
                [IMAGES] = contPIV.imageLoader_RS(imgPath,'tif');
                fprintf('Image loaded in %f s \n',toc(tt));
                tt=tic;
                contPIV.wrapperPIVSPMD(rootdir,alias,cfg_data,filoutVEL,ifile,IMAGES,imgPath);
                fprintf('PIV completed in %f s \n',toc(tt));
                IMAGES = [];
            end

        end
end

timeStop = toc(timeStart)


%% Serial test
clearAllMemoizedCaches;

serialtimestart = tic;
for ifile = 1 : 20
    alias = t.Alias{ifile};
    imgPath = t.Location{ifile};
    tt=tic;
    [IMAGES] = contPIV.imageLoader_RS(imgPath,'tif');
    fprintf('Image loaded in %f s \n',toc(tt));

    
    tt=tic;
    contPIV.wrapperPIVSPMD(rootdir,alias,cfg_data,filoutVEL,ifile,IMAGES,imgPath);
    fprintf('PIV completed in %f s \n',toc(tt));
end
serialtimeduration = toc(serialtimestart) 