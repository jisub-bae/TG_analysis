function [si, di, time] = calculateSDI(FittedIEM, toi, zeropoint,PermutedIEM)
    
    sw = 100;
    corr = 90;
    
    time = linspace(toi(1),toi(2), size(FittedIEM,3));
    FittedIEM-0.5;
    
    if nargin < 4
        PermutedIEM = [];
    end
    
    if isempty(PermutedIEM)
        [H1, H2, H3] = getH123FromT(FittedIEM,zeropoint);
    else
        [H1, H2, H3] = getH123FromPermuted(FittedIEM, PermutedIEM);
    end
    
    index = toi(1) <= time & time <= toi(2);
    time = time(index);
    H1 = H1(index,index);
    H2 = H2(index,index);
    H3 = H3(index,index);
    
    Hd = H1 .* H2;
    Hs = (1-H1) .* (1-H2) .* H3;
    
    Hd = squeeze(Hd);
    Hs = squeeze(Hs);
    
    nSamples = min([size(Hd) size(Hs)]);
    
    nanConv = ones(nSamples); 
    for i1 = 1:nSamples
        for i2 = 1:nSamples
            if abs(i1-i2) <= corr
                nanConv(i1,i2) = nan;
            end
        end
    end
    Hd = Hd .* nanConv;
    Hs = Hs .* nanConv;
    
    Hd = padarray(Hd, [sw sw], nan);
    Hs = padarray(Hs, [sw sw], nan);
    
    di = nan(1,nSamples); si = nan(1,nSamples);
    for i1 = (1:nSamples)+sw
        si(i1-sw) = nanmean(nanmean(Hs(i1-sw:i1+sw, i1-sw:i1+sw),1),2);
        di(i1-sw) = nanmean(nanmean(Hd(i1-sw:i1+sw, i1-sw:i1+sw),1),2);
    end
    
end