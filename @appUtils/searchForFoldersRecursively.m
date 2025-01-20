% Key is case insensitive
% avoidKey is case sensitive

function [finalList,nameList] = searchForFoldersRecursively(folderPath,key,avoidKey)
    
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
    
    for i = 1 : length(aDir)
        
        if regexpi(aDir(i).name, regexptranslate('wildcard',key)) % If there is a match
            
                if (nargin == 3 ) % If avoidkey is specified
                    if isempty(regexp(aDir(i).name, regexptranslate('wildcard',avoidKey)))  
                        finalList = [finalList {fullfile(aDir(i).folder,aDir(i).name)}];
                        nameList = [nameList {aDir(i).name}];
                    else
                        [tmpFinal, tmpName] = searchForFoldersRecursively(fullfile(aDir(i).folder,aDir(i).name),key,avoidKey);
                        finalList = [finalList tmpFinal];
                        nameList = [nameList tmpName];
                    end
                else
                    finalList = [finalList {fullfile(aDir(i).folder,aDir(i).name)}];
                    nameList = [nameList {aDir(i).name}];
                end               
            
        else % Go deeper
            if (nargin == 3 ) % If avoidkey is specified
                    [tmpFinal, tmpName] = appUtils.searchForFoldersRecursively(fullfile(aDir(i).folder,aDir(i).name),key,avoidKey);
                    finalList = [finalList tmpFinal];
                    nameList = [nameList tmpName];
                else
                    [tmpFinal, tmpName] = appUtils.searchForFoldersRecursively(fullfile(aDir(i).folder,aDir(i).name),key);
                    finalList = [finalList tmpFinal];
                    nameList = [nameList tmpName];
            end    
        end
    end
    
    if isempty(aNotDir)
        return
    end
       
end

    