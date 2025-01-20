classdef contUtils 

methods (Static = true)


    function [avMem] = getAvailableRAM
        % In GB
        if ~isunix % Windows
            t = memory;
            avMem = t.MaxPossibleArrayBytes * 1e-9;
            avMem = avMem * 0.9; % Heuristic
        else
            [r,w] = unix('free | grep Mem');
            stats = str2double(regexp(w, '[0-9]*', 'match'));
            avMem = stats(end)/1e6; % Report is in kBs
            avMem = avMem * 0.9;
        end
    end
    
    function value = isNaturalNumber(number)
        value = number>0 & mod(number,1)==0;
    end







end

% I/O methods

methods (Static = true)

    % Logging 
    function writeToLog(logfile,msg)
        fid = fopen(logfile,'a');
        fprintf(fid,[msg,' ',datestr(now),'\n']);
        fclose(fid);
    end

    function [cfg_data] = loadJsonConfig(configfile)
        fid = fopen(configfile);
        raw = fread(fid,inf);
        str = char(raw');
        fclose(fid);
        cfg_data = jsondecode(str);
    end

    % Default contractility config file
    function [cfg_data] = loadDefaultConfig()
        configfile = 'default_config.json';
        fid = fopen(configfile);
        raw = fread(fid,inf);
        str = char(raw');
        fclose(fid);
        cfg_data = jsondecode(str);
    end


    % Read results
        [X,Y,T,U,V,Xdrift,Ydrift] = readPIV_bin(filename);
    

    
    
        [X,Y,T,SXX,SYY,SXY,Xdrift,Ydrift] = readMSM_bin(filename);
    




end





end