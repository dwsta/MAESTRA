classdef appUtils

methods (Static = true)


    % Recursively search folders and files 
   
    [finalList,nameList] = searchForFilesRecursively(folderPath,key,avoidKey,doRecurse);
    [finalList,nameList] = searchForFoldersRecursively(folderPath,key,avoidKey);

    % Change path in jobfile to match the right location in the current drive setup
    function jobfile = checkVideoLocation(jobfile,newPath);
        for ii = 1 : height(jobfile)
            [newStr] = appUtils.correctForAbsolutePath(jobfile.Location{ii},newPath);
            jobfile.Location{ii} = newStr;
        end
    end

    [newStr] = correctForAbsolutePath(oldPath,newPath);


end





end