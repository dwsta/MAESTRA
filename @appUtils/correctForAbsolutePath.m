function [newStr] = correctForAbsolutePath(oldPath,newPath)
% Change the path to videos by matching pattern

if contains(oldPath,'\')
    oldFileSep = '\';
else
    oldFileSep = '/';
end

if contains(newPath,'\')
    newFileSep = '\';
else
    newFileSep = '/';
end

oldSplit = split(oldPath, oldFileSep);
newSplit = split(newPath, newFileSep);
mergeWords = intersect( oldSplit , newSplit );
mergeWords(cellfun(@isempty,mergeWords)) = [];

if ~isempty(mergeWords)
    oldInd = find(strcmp(oldSplit,mergeWords{1}));
    newInd = find(strcmp(newSplit,mergeWords{1}));

    newStr = strjoin([newSplit(1:newInd(1)-1); oldSplit(oldInd(1):end)],filesep);
else
    warning('No overlapping path found');
    return;
end



