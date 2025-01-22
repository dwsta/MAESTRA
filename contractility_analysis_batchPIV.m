clearvars;
restoredefaultpath;
addpath 'D:\Cloud\Git\MAESTRA';

% Read the jobfile
% rootdir = 'E:\Contractility_3\2023-03-25\treatment_30min_2023-03-25\contractility_run_20241124_130815_WL64\';
% rootdir = 'E:\Contractility_3\2023-03-25\treatment_30min_2023-03-25\contractility_run_20241110_194855_beads_driftCorrected';
% rootdir =
% 'E:\Contractility_3\2023-03-25\pretreatment_2023-03-25\test_batch';
rootdir = 'D:\Data\Contractility_batch_test\contractility_run_20250121_210036\';
cd(rootdir);

t = readJobFile2('jobfile.csv', rootdir);

% Make sure the path in jobfile is correct
t = appUtils.checkVideoLocation(t,rootdir);

cfg_data = contUtils.loadJsonConfig(fullfile(rootdir,'config.json'));
cfg_GPU = cfg_data; 
cfg_CPU = cfg_data;
cfg_CPU.Performance.UseGPU = false; 


filoutVEL = fullfile(rootdir,'vel');
if (cfg_data.ReferenceSelection.ReferenceSelectionMode ~= 0)
    if ~exist(filoutVEL,'dir'); mkdir(filoutVEL); end
end

warning('off');
% 
% N = length(files);
% w = waitbar(0,'PIV...');
% w.UserData = N;

% cCPU = parcluster('local');


% cGPU = parcluster('local');


% job = createJob(c);

p = parpool('local',2);
%%
spmd
    if labindex == 1
        gpuDevice(1);
    elseif labindex == 2
        gpuDevice([]);
    end
end

% try

%%


for ii = 1 : 2
    
    f(ii) = parfeval(@contPIV.wrapperBatchPIV,1, rootdir,cfg_GPU,filoutVEL,t,ii);

end












%% Separate CPU and GPU queue 
% 
% tt=tic; 
% cpuFilInd = [1:10];
% jobCPU = createJob(cCPU);
% createTask(jobCPU,@contPIV.wrapperBatchPIV, 0, {rootdir,cfg_CPU,filoutVEL,t,cpuFilInd});
% % jobCPU = batch(cCPU,@contPIV.wrapperBatchPIV, 0, {rootdir,cfg_CPU,filoutVEL,t,cpuFilInd},'Pool',2);
% 
% 
% gpuFilInd = [11:20];
% jobGPU = createJob(cGPU);
% createTask(jobGPU,@contPIV.wrapperBatchPIV, 0, {rootdir,cfg_GPU,filoutVEL,t,gpuFilInd});
% 
% % batch(cGPU,@contPIV.wrapperBatchPIV, 0, {rootdir,cfg_GPU,filoutVEL,t,gpuFilInd},'Pool',2);
% 
% % end
% 
% submit(jobCPU);
% submit(jobGPU);
% 
% while 1
%    if wait(jobCPU) & wait(jobGPU)
%        break
%    end
%    pause(1);
% end
% toc(tt)
% % fclose(csv_whole_FOV);





function [aa] = getFileName(x)
[~,aa] = fileparts(fileparts(fileparts(x)));
end


