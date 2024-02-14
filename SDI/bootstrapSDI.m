function [se_si, se_di,data_si,data_di] = bootstrapSDI(FittedIEM, toi,zeropoint, varargin)

    p = inputParser;
    addRequired(p, 'FittedIEM');
    addOptional(p, 'toi', []);
    addParameter(p, 'nIterations', 10000);
    parse(p, FittedIEM, varargin{:});

    nIterations = p.Results.nIterations;
    
    testTime = linspace(toi(1),toi(2), size(FittedIEM,3));
        
    nIEMs =size(FittedIEM,1);
    nTimes = sum(toi(1) <= testTime & testTime <= toi(2));
    bootstrappedSI = zeros(nIterations,nTimes);
    bootstrappedDI = zeros(nIterations,nTimes);
    parfor iIterations = 1:nIterations
        [bootstrappedSI(iIterations,:), bootstrappedDI(iIterations,:)] ...
            = calculateSDI(FittedIEM(randi(nIEMs, nIEMs, 1),:,:), toi,zeropoint);
    end
    se_si = squeeze(std(bootstrappedSI));
    se_di = squeeze(std(bootstrappedDI));
    data_si=bootstrappedSI;
    data_di=bootstrappedDI;
    
end

