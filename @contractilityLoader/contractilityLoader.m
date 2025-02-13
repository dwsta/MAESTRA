classdef contractilityLoader < handle

    properties (Access = public)
        plate = []; % plateReader obj
        dataPath = [];
        data = []; % Data table (Results)
        aliasKey = [];%'(\d+)_(\d+)_(\d+)_(\w+)_(\w+)_Well__(\w)_(\d+)'; % Default pattern
                % E.g alias : 2021_12_13_Pretreatment_BottomLeft_Well__B_002_Fitc_2
                % (YYYY)_(MM)_(DD)_(condId)_(posId)_(Well)__(WellId)_(WellNo)_##' 
                % String inside parantheses are processed as
                % literals 
                %
                % To use a different format, create a children class and write a
                % new method 'parseAlias' to accomodate.

    end


    methods(Access = public)

        % Constructor
        function obj = contractilityLoader(dataPath, plate)

            obj.data = readtable(dataPath,'Delimiter',',');
            obj.data = obj.getTimeSeriesFromCSV(obj.data);
            obj.dataPath = dataPath;

            if nargin == 2
                obj.plate = plate;
                obj = obj.addCompoundsToTable(plate);
            end

        end

        

        function sub = filter(obj, alias, compound, concs,  condId, posId)
            % sub = filter(obj, alias, compound, concs,  condId, posId)

            indsAll = cellfun(@(x) ~isempty(x),obj.data.alias); % All true
            indsNone = zeros(height(obj.data),1); % All false

            inds = indsAll;

            if ~isempty(alias) & sum(ismember(fieldnames(obj.data),'alias'))
                   inds = inds & contains(obj.data.alias,alias);
            end
            
            indsT = indsNone;
            if ~isempty(condId) & sum(ismember(fieldnames(obj.data),'condId'))
                    for ii = 1 : length(condId); indsT = indsT | strcmpi(obj.data.condId,condId{ii}); end;
                    inds = inds & indsT;
            end
            

            indsT = indsNone;
            if ~isempty(posId) & sum(ismember(fieldnames(obj.data),'posId'))
                for ii = 1 : length(posId); indsT = indsT | strcmpi(obj.data.posId,posId{ii}); end;
                inds = inds & indsT; 
            end
            

            if ~isempty(obj.plate) % If platemap is present
                indsT = indsNone;
                if ~isempty(compound)
                    for ii = 1 : length(compound); indsT = indsT | strcmpi(obj.data.compound,compound{ii}); end;
                    inds = inds &indsT; 
                end
                
                indsT = indsNone;
                if ~isempty(concs)
                    for ii = 1 : length(concs); indsT = indsT | strcmpi(obj.data.concs,concs{ii}); end;
                    inds = inds & indsT; 
                end
                

            else % If absent, sub with well and well No
                warning('No platemap present. Trying with well and well # as first two arguments.');

                indsT = indsNone; 
                if ~isempty(compound) 
                    for ii = 1 : length(compound); indsT = indsT | strcmpi(obj.data.wellId,wellId{ii}); end;
                    inds = inds & indsT;
                end

                indsT = indsNone;
                if ~isempty(concs)
                    for ii = 1 : length(concs); indsT = indsT | strcmpi(obj.data.concs,num2str(wellNo{ii})); end;
                    inds= inds & indsT;
                end
                
            end
            
            sub = obj.data(inds,:);
        end

    end


    methods (Access = private)

        function obj = addCompoundsToTable(obj,plate)
            
            obj.data = obj.parseAlias(obj.data, obj.aliasKey);    

            comp = arrayfun(@(x,y) plate.getCompound(x,y), obj.data.wellId,obj.data.wellNo);
            concs = arrayfun(@(x,y) plate.getConc(x,y), obj.data.wellId,obj.data.wellNo);
            % concs = cellfun(@(x,y) getWellConcs2023(x,y,1),string(wellId),num2cell(wellNo));
        
            obj.data.compound = comp;
            obj.data.concs = concs;
        end

       
    end

    methods (Static = true)

            function data = parseAlias(data,pattern,orderNo)
                % pattern - '(\d+)_(\d+)_(\d+)_(\w+)_(\w+)_Well__(\w)_(\d+)'
                % Order : wellId, wellNo, condId, posId
                if ~isempty(pattern)
                    names = regexp( data.alias, pattern, 'tokens','once');
                    if ~isnan(orderNo(1)); data.wellId = cellfun(@(x) x{orderNo(1)},names); end;
                    if ~isnan(orderNo(2)); data.wellNo = cellfun( @(x) str2num( string(x{orderNo(2)})),names); end;
                    if ~isnan(orderNo(3)); data.condId = cellfun( @(x)  string(x{orderNo(3)}),names); end; 
                    if ~isnan(orderNo(4)); data.posId = cellfun( @(x) string(x{orderNo(4)}),names); end;
                else
                    pattern = 'Well__(\w)_(\d+)';
                    names = regexp( data.alias, pattern, 'tokens','once');
                    data.wellId = cellfun(@(x) x{1},names);
                    data.wellNo = cellfun( @(x) str2num( string(x{2})),names);
                end

            end

    
            function res = getTimeSeriesFromCSV(res)
    
                for ii = 1: height(res)
                    try; res.time{ii} = convertCharToCell(res.time{ii}); end;
                    try; res.signal{ii} = convertCharToCell(res.signal{ii}); end;
                    try; res.peak_ID{ii} = convertCharToCell(res.peak_ID{ii}); end;
                end
                    
                    function tt=convertCharToCell(tt)
                        tt=eval(tt);
                        tt=[tt{:}];
                    end
    
            end

            function subplotSignal(sub,m,n)
                if height(sub)>(m*n)
                    warning('More rows than m x n. Plotting first m x n. ');
                    % sub=sub(randperm(height(sub), m*n),:); 
                    sub = sub(1:m*n,:);
                end
                maxY = cellfun(@(x) prctile(x,99),sub.signal);
                maxY = max(maxY) + 0.1*range(maxY);
                % minY = min(cellfun(@(x) prctile(x,10),sub.signal));
                minY = 0;
                clf; hold on; 
                nn = min(height(sub), m*n);
                for ii=1:nn; subplot(m,n,ii); plot(sub.signal{ii}); ylim([minY maxY]);
                % title(sub.alias{ii}); 
                title(sprintf('Well-%s%0.3d - %s',sub.wellId(ii),sub.wellNo(ii), sub.posId{ii}),'Interpreter','none');
                end
            end

    end

    methods (Access = public)
        function xx=addBoxToAxes(obj,metric,alias,compound,conc,condId,posId)
                % addBoxToAxes(obj,metric,alias,compound,conc,condId,posId)
                arguments (Input)
                    obj (1,1) contractilityLoader
                    metric (1,1) string
                    alias (1,:) cell = {}
                    compound (1,:) cell = {}
                    conc = {}
                    condId  (1,:) = {} 
                    posId (1,:) = {}
                end
                sub = filter(obj, alias, compound, conc, condId, posId);
                xx = [sub.(metric)];
                % Xpos = length(get(gca,'Children'))*max(abs(xx(:)))+max(abs(xx(:)))*50;
                Xpos = length(get(gca,'Children'))+1;
                gg = xx.*0 + Xpos;
                boxplot(xx,gg,'Positions',Xpos); hold on; 
                scatter(gg+(rand(size(gg))-0.5)/4,xx,'filled','MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',0.4);
                axis square;
        end
    end

end









