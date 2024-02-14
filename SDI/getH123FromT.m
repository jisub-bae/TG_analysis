function [H1, H2, H3] = getH123FromT(TargetIEM,zeropoint)
     
    nTimes = size(TargetIEM,3);

    slopes = TargetIEM;
    diag = repmat(slopes(:,1:nTimes+1:end), [1 1 nTimes]);
    
    alpha = 0.05;
    H1 = squeeze(ttest(diag, slopes, 'tail', 'right', 'alpha', alpha)) == 1;
    H2 = squeeze(ttest(permute(diag, [1 3 2]), slopes, 'tail', 'right', 'alpha', alpha)) == 1;
    H3 = squeeze(ttest(slopes, zeropoint, 'tail', 'right', 'alpha', alpha)) == 1;
    
end