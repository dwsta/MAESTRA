classdef plateReader

    % Plate (csv) of format : 
    %
    % Drug 1 2 3 ... (Title row)
    % A    # # # ... (wells)
    % ...    ...    ... 
    % 
    % Concs(uM) ... (Concs title)
    % A    # # # ... (Wells) 
    % ...    ...    ...



    properties 

        comps = [];
        concs = [];
        wells = [];
        wellNos = [];
        platePath = [];

    end

    methods (Access = public)
        function obj = plateReader(fPath)

            [comps,concs,wells,wellNos] = plateReader.readPlate(fPath);
            obj.comps = comps;obj.concs = concs;
            obj.wells = wells; obj.wellNos = wellNos';

            obj.platePath = fPath;

        end
        

        function tt = getConc(obj,well,wellNo)
            % tt = getConc(obj,well,wellNo)
            
            assert(~isempty(obj.concs),'No concentration platemap found.');
            assert(~isempty(wellNo),'Need to supply well number as the 3rd argument.');

            if nargin == 2 
                [well,wellNo] = plateReader.parseWellText(well);
            end
            
            tt = arrayfun(@(x,y) obj.concs(x,:).(y),well,string(wellNo));
            
        end

        function tt = getCompound(obj, well, wellNo)
            
            if nargin == 2 
                [well,wellNo] = plateReader.parseWellText(well);
            end

            if isempty(wellNo)
               wellNo = [obj.wellNos]; 
               well = repmat(well,length(wellNo),1);
               findUniqueName = true;
            else
               findUniqueName = false;
            end
            
            assert(length(well) == length(wellNo));
            
            tt = arrayfun(@(x,y) obj.comps(x,:).(y),well,string(wellNo));

            if findUniqueName
                tt(ismissing(tt)) = [];
                tt = unique(tt);
            end


        end

        function [wells,wellNos] = getWell(obj, compound, conc)
            
            maskComp = cellfun(@(x) strcmp(x,compound), table2cell(obj.comps));

            if nargin > 2 & ~isempty(obj.concs)
                maskConc = cellfun(@(x) strcmp(x,string(conc)), table2cell(obj.concs));
            else
                assert(isempty(conc),'No concentration platemap found. Use [] to return all rows.');
                maskConc = logical(ones(size(maskComp)));
            end
            
            mask = maskConc & maskComp; 
            
            [r,c] = ind2sub(size(mask), find(mask));

            wells = obj.wells(r); wellNos = obj.wellNos(c);

        end

    end

    methods (Static = true)
        
        function [comps,concs,wells,wellNos] = readPlate(fPath)


            raw = readcell(fPath); 
            
            % Replace empty rows (<missing>) with nans
            % mask = cellfun(@(x) any(isa(x,'missing')), raw);
            raw = string(raw);


            [plateFormat,rows,cols] = findPlateFormat(raw);
            
            wells = char(64+[1:rows]');
            wellNos = [1:cols];

            comps = array2table(  raw(2:1+rows,2:1+cols) ,'VariableNames',string(wellNos));
            comps.Properties.RowNames = string(wells);
            
            if height(raw) > 2 + rows 
                concs = array2table( raw(rows+3:2*(1+rows),2:1+cols) ,'VariableNames',string(wellNos));
                concs.Properties.RowNames = string(wells);
            else
                concs = [];
            end

                function [plateFormat,rows,cols] = findPlateFormat(raw)
        
                    % 96 well plate is 8 x 12
                    % 384 well plate is 16 x 24
                
                    plateFormat = nan; 
                    switch size(raw,2)-1 
                        case 12
                            
                            dupl = determineDupl(size(raw,1),8);
                            
                            if ~isnan(dupl)
                                plateFormat = 96;
                                rows = 8; 
                                cols = 12;
                            end
                
                        case 16
                            dupl = determineDupl(size(raw,1),8);
                            
                            if ~isnan(dupl)
                                plateFormat = 384;
                                rows = 16;
                                cols = 24;
                            end
                        otherwise 
                            error('Wrong number of rows and columns detected');
                    
                    end
                
                
                end
                
                
                function dupl = determineDupl(siz,rows)
                    % dupl - Number of platemaps (i.e. compound name and concs - dupl = 2)
                    
                    % siz/x - 1 == rows
                
                    dupl = siz/(rows+1);
                
                    if mod(dupl,1) ~= 0 % not a whole number
                        dupl = nan;
                        error('Row formatting is incorrect. Cannot determine the number of rows. See determineDupl')
                    end
                
                end

        end

        function [well,wellNo] = parseWellText(well)
            % If format 'Well__B_002' is used (no wellNo supplied)
            pattern = '(\w+)?_(\w+)_(\d+)';
            tt=regexp(well,pattern,'tokens','once');
            well = tt{2};
            wellNo = str2double(tt{3});
        end
    end



methods (Access = private)

    

end


end