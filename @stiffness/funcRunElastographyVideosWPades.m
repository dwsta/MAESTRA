function [EArr,nuArr,GArr,KArr,LinE,NH] = funcRunElastographyVideosWPades(X,Y,UC,VC,UB,VB,inds,runLinE,runNH)

    % Wrapper for elastography of cardiomyocyte videos
    
    
    % Assume the input displacement vector field are already filtered
    dx = X(1,2) - X(1);
    dy = Y(2,1) - Y(1);
    
    UBF = nan([ size(UB,[1 2]) length(inds)]); VBF = UBF;
    UCF = UBF; VCF = VBF;
    
    for ii = 1 : length(inds)
        UBF(:,:,ii) = UB(:,:,inds(ii));
        VBF(:,:,ii) = VB(:,:,inds(ii));
        UCF(:,:,ii) = UC(:,:,inds(ii));
        VCF(:,:,ii) = VC(:,:,inds(ii));
    end

    if nargin <= 5
        runNH = 1;
        runLinE = 0;
    end
    
    [LinE, NH] = runOptimNHLinEWPades(UCF,VCF,UBF,VBF,X,Y,runLinE,runNH);
    
    
    
    EArr=[];nuArr=[];
    GArr=[];KArr=[];
    
    if runLinE
        for i = 1: length(inds)
            EArr(i)=LinE(i).sol.E;
            nuArr(i)=LinE(i).sol.nu;
        end
    end
    
    if runNH
        for i = 1: length(inds)
            GArr(i)=NH(i).sol.G;
            KArr(i)=NH(i).sol.K;
        end 
    end

end
