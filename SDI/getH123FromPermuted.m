function [H1, H2, H3] = getH123FromPermuted(ObservedIEM, PermutedIEM)
    
    
    nSubjects = length(ObservedIEM);
    nTimes = length(ObservedIEM(1).testTime);
    
    slopes = nan(nSubjects, nTimes, nTimes, 'single');
    for s1 = 1:nSubjects
        slopes(s1,:,:) = ObservedIEM(s1).slope;
    end
    
    diag = nan(nSubjects, nTimes, 1, 'single');
    for i1 = 1:nTimes
        diag(:,i1,1) = slopes(:,i1,i1);
    end
    diag = repmat(diag, [1 1 nTimes]);
    
    sqrtn = sqrt(nSubjects);
    tObserved = mean(slopes) ./ ( std(slopes)/sqrtn );
    tDiagObserved = mean(diag-slopes) ./ ( std(diag-slopes)/sqrtn );
    diag_t = permute(diag, [1 3 2]);
    tDiagObserved_t = mean(diag_t-slopes) ./ ( std(diag_t-slopes)/sqrtn );
    
    P1 = zeros(size(tObserved), 'single'); P2 = P1; P3 = P1;
    nPermutations = size(PermutedIEM(1).slope,1);
    for iPermutations = 1:nPermutations
        permutedSlopes = nan(size(slopes), 'single');
        for s1 = 1:nSubjects
            permutedSlopes(s1,:,:) = PermutedIEM(s1).slope(iPermutations,:,:);
        end
        tPermuted = mean(permutedSlopes) ./ ( std(permutedSlopes)/sqrtn );
        permutedDiag = nan(nSubjects, nTimes, 1, 'single');
        for i1 = 1:nTimes
            permutedDiag(:,i1,1) = permutedSlopes(:,i1,i1);
        end
        permutedDiag = repmat(permutedDiag, [1 1 nTimes]);
        tDiagPermuted = mean(permutedDiag-permutedSlopes) ./ ( std(permutedDiag-permutedSlopes)/sqrtn );
        permutedDiag_t = permute(permutedDiag, [1 3 2]);
        tDiagPermuted_2 = mean(permutedDiag_t-permutedSlopes) ./ ( std(permutedDiag_t-permutedSlopes)/sqrtn );
        
        P1 = P1 + (tDiagPermuted > tDiagObserved);
        P2 = P2 + (tDiagPermuted_2 > tDiagObserved_t);
        P3 = P3 + (tPermuted > tObserved);
    end
    P1 = P1 / nPermutations;
    P2 = P2 / nPermutations;
    P3 = P3 / nPermutations;
    
    H1 = squeeze(logical(P1 < 0.05));
    H2 = squeeze(logical(P2 < 0.05));
    H3 = squeeze(logical(P3 < 0.05));
end