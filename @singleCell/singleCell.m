classdef singleCell

methods (Static = true)


    function folderPath = getCellposeModelPath
        folderPath = mfilename("fullpath");
        folderPath = fullfile( fileparts(folderPath) ,... 
                            'cellpose_model');
    end

end




end