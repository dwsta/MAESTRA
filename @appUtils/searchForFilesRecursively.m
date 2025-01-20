% Search recursively for files by pattern matching down the diretory tree %
% 
% INPUTS :
% folderPath - Char array - Parent directory path
% key - char array - Pattern to search for in file name
% avoidKey - char array - (key) & ~(avoidKey)
% 
% OUTPUTS :
% finalList - Cell array of FULL path to files
% nameList - Cell array of file names with extensions
% 
%
% Adithan Kandasamy
% 2021, UW

function [finalList,nameList] = searchForFilesRecursively(folderPath,key,avoidKey,doRecurse)
    
    if nargin  < 4
        doRecurse = 1;
    end
    finalList = {};
    nameList = {};
    a = dir(folderPath);
    aDir = a(horzcat(a.isdir));
    % Clean up directoies to only on that level
    tmp = struct2cell(aDir);
    % ind = find(~contains(tmp(1,:),["." "@"])); % Remove level-up directories and @ directories
    ind = ~contains(tmp(1,:),["@"]); % Remove level-up directories and @ directories
    ind = ind & ~strcmpi(tmp(1,:),".");
    ind = ind & ~strcmpi(tmp(1,:),"..");
    ind = find(ind);
    aDir = aDir(ind); % Directories
    
    aNotDir = a(~horzcat(a.isdir)); % Files
    
    for i = 1 : length(aNotDir)
        if (regexp(aNotDir(i).name, regexptranslate('wildcard',key))) % If there is a match
            if (nargin == 3 ) % If avoidkey is specified
                if isempty(regexp(aNotDir(i).name, regexptranslate('wildcard',avoidKey)))  
                    finalList = [finalList {fullfile(aNotDir(i).folder,aNotDir(i).name)}];
                    nameList = [nameList {aNotDir(i).name}];
                end
            else
                finalList = [finalList {fullfile(aNotDir(i).folder,aNotDir(i).name)}];
                nameList = [nameList {aNotDir(i).name}];
            end
        end
    end
    
    if ~doRecurse
        return;
    end
    
    if isempty(aDir) 
        return
    else
        for i = 1 : length(aDir)
            if (nargin == 3 ) % If avoidkey is specified
                [tmpFinal, tmpName] = appUtils.searchForFilesRecursively(fullfile(aDir(i).folder,aDir(i).name),key,avoidKey);
                finalList = [finalList tmpFinal];
                nameList = [nameList tmpName];
            else
                [tmpFinal, tmpName] = appUtils.searchForFilesRecursively(fullfile(aDir(i).folder,aDir(i).name),key);
                finalList = [finalList tmpFinal];
                nameList = [nameList tmpName];
            end
        end
    end
    
end

        
        
    
    
    
    