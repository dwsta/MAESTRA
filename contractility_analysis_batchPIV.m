clearvars;
restoredefaultpath;
addpath 'D:\Code\Contractility_GitHub - Adithan - Final';

% Read the jobfile
% rootdir = 'E:\Contractility_3\2023-03-25\treatment_30min_2023-03-25\contractility_run_20241124_130815_WL64\';
% rootdir = 'E:\Contractility_3\2023-03-25\treatment_30min_2023-03-25\contractility_run_20241110_194855_beads_driftCorrected';
rootdir = 'E:\Contractility_3\2023-03-25\pretreatment_2023-03-25\test_batch';
cd(rootdir);

t = readJobFile2('jobfile.csv', rootdir);

% Make sure the path in jobfile is correct
t = appUtils.checkVideoLocation(t,rootdir);

cfg_data = contUtils.loadJsonConfig(fullfile(rootdir,'config.json'));

filoutVEL = fullfile(rootdir,'vel');
if (cfg_data.ReferenceSelection.ReferenceSelectionMode ~= 0)
    if ~exist(filoutVEL,'dir'); mkdir(filoutVEL); end
end

warning('off');
% 
% N = length(files);
% w = waitbar(0,'PIV...');
% w.UserData = N;

c = parcluster('local');
% job = createJob(c);

% try

%%

% Separate CPU and GPU queue 
    for ifile = 1 : 4 % length(files)
        % tt=tic; 
        
        % w=waitbar(ifile/length(files));
        
        alias = t.Alias{ifile};
        imgPath = t.Location{ifile}; ext = 'tif';
        
        % contPIV.wrapperPIV(rootdir,alias,cfg_data,filoutVEL,ifile,imgPath,ext);

        job(ifile) = batch(c,@contPIV.wrapperPIV, 0, {rootdir,alias,cfg_data,filoutVEL,ifile,imgPath,ext},'Pool',3);
        

        % wholeFOV.whole_FOV_analysis(rootdir,alias,cfg_data,csv_whole_FOV);
        % toc(tt)
    end
% end

% submit(job);
    
% fclose(csv_whole_FOV);





function [aa] = getFileName(x)
[~,aa] = fileparts(fileparts(fileparts(x)));
end


