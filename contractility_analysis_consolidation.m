clearvars;
restoredefaultpath;
addpath 'D:\Code\Contractility_GitHub - Adithan - Final';

% Read the jobfile
% rootdir = 'E:\Contractility_3\2023-03-25\treatment_30min_2023-03-25\contractility_run_20241124_130815_WL64\';
% rootdir = 'E:\Contractility_3\2023-03-25\treatment_30min_2023-03-25\contractility_run_20241110_194855_beads_driftCorrected';
rootdir = 'E:\Contractility_3\2023-03-25\treatment_30min_2023-03-25\dynamic_contractility_run_20241015_164036_Mito';

t = readJobFile2('jobfile.csv', rootdir);

whichPass = 2;

% Search for deformations_pass_2.bin 
parentPath = fullfile(rootdir, 'output');
[files] = appUtils.searchForFilesRecursively(parentPath, sprintf('smooth_deformations_pass%d.bin',whichPass),'',1);
files(contains(files,"_001") | contains(files,"_012")) = []; % No idea what 001 and 012 wells are 

% files(~contains(files,'Well__C_009_r_0005_c_0005_beads')) = [];

pivFolders = cellfun(@(x) fileparts(x),files,'UniformOutput',false); 
fileNames = cellfun(@(x) getFileName(x),files,'UniformOutput',false);

% Rewrite config
cfg_data = contUtils.loadJsonConfig(fullfile(rootdir,'config.json'));

% Create output file
output_file_connections = create_csv(cfg_data,rootdir,whichPass);
csv_whole_FOV = output_file_connections(1);

warning('off');

N = length(files);
w = waitbar(0,'Consildating...');
w.UserData = N;
%%
% try
    for ifile = 1 : length(files)
        
        w=waitbar(ifile/length(files));
        
        alias = fileNames{ifile};
    
        wholeFOV.whole_FOV_analysis(rootdir,alias,cfg_data,csv_whole_FOV,whichPass);
    
    

    end
% end
    
fclose(csv_whole_FOV);
close(w)


function [aa] = getFileName(x)
[~,aa] = fileparts(fileparts(fileparts(x)));
end


function output_file_connections = create_csv(cfg_data,rootdir,whichPass)
% Helper function to create the csv output files and remove clutter from the main loop
% For whole FOV
if cfg_data.TFM.DoTFM
    outcsv = fullfile(rootdir,sprintf(['contractility_motion_pass%d_',datestr(now,30),'.csv'],whichPass));
    headers = {'alias','signal','time','peak_ID','mean_amplitude','mean_peak_duration','mean_rise_time', 'mean_fall_time', 'mean_pw90', 'mean_pw50', 'mean_pw30', 'mean_peak_value', 'mean_baseline', 'mean_valley',...
    'mean_E','mean_nu','mean_G', 'mean_K', ...
    'sd_amplitude', 'sd_peak_duration' ,'sd_rise_time', 'sd_fall_time', 'sd_pw90', 'sd_pw50', 'sd_pw30', 'sd_peak_value', 'sd_baseline', 'sd_valley',...
    'sd_E','sd_nu','sd_G', 'sd_K'};
else
    outcsv = fullfile(rootdir,['contractility_motion_',datestr(now,30),'.csv']);
    headers = {'alias','signal','time','peak_ID','mean_amplitude','mean_peak_duration','mean_rise_time', 'mean_fall_time', 'mean_pw90', 'mean_pw50', 'mean_pw30', 'mean_peak_value', 'mean_baseline', 'mean_valley',...
    'sd_amplitude', 'sd_peak_duration' ,'sd_rise_time', 'sd_fall_time', 'sd_pw90', 'sd_pw50', 'sd_pw30', 'sd_peak_value', 'sd_baseline', 'sd_valley'};

end
output_file_connections(1) = fopen(outcsv,'w');
fprintf(output_file_connections(1),'%s',strjoin(headers,','));
% For cell by cell

end