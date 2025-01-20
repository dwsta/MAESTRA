function [IMAGES] = imageLoader(imgPath,ext,precision,parRun)

arguments

    imgPath (1,1) string {mustBeNonempty}
    ext (1,1) string = "tif";
    precision {mustBeMember(precision,{'single','double'})} = 'single' ;
    parRun (1,1) logical = true; % true will use parfor for loading images, abour 25-40% faster

end

% For now we only accept videos as image sequences
[files] = appUtils.searchForFilesRecursively(imgPath,sprintf('*.%s',ext),'',0);
if isempty(files)
    warning(sprintf('No images found in %s',imgPath));
    IMAGES = []; 
    return;
end

% Find the size of the images
inf = imfinfo(files{1});
% Pre-allocation as single or double, faster than allocating as uint and
% converting later
IMAGES = zeros(inf.Width, inf.Height, length(files), precision);

if parRun
    parfor ii = 1 : length(files)
        IMAGES(:,:,ii) = single(imread(files{ii}));
    end
else
    for ii = 1 : length(files)
        IMAGES(:,:,ii) = single(imread(files{ii}));
    end
end

end