function [ICPRfactors, harmonicsNumbers, harmonicsIndexes] = ICPRcombinations(harmoNums, maxCombinationsLength)
    ICPRfactors = []; harmonicsNumbers = []; harmonicsIndexes = [];
    if nargin<2, maxCombinationsLength = Inf; end
    for ai = 2:min([length(harmoNums) maxCombinationsLength])
        combinations = nchoosek(1:numel(harmoNums), ai); %Combinate indexes.
        combinations = arrayfun(@(x) combinations(x, :), 1:size(combinations, 1), 'UniformOutput', false)';
        %..Leave only unique combinations with sum of harmonic numbers 0..
        numCombs = cellfun(@(x) harmoNums(x), combinations, 'UniformOutput', false);
        ICPRrels = cellfun(@(x) ones(size(x)), combinations, 'UniformOutput', false);
        validIdxs = false(size(numCombs));
        for bi = 1:length(numCombs)
            %Sort combinations.
            [numCombs{bi}, idxs] = sort(numCombs{bi});
            combinations{bi} = combinations{bi}(idxs);
            if ai == 2, ICPRrels{bi}(2) = numCombs{bi}(1)/numCombs{bi}(2);
            else, ICPRrels{bi}([1 3]) = [0.5 0.5];
            end
            ICPRrels{bi}(2) = -ICPRrels{bi}(2);
            %Check the current combination.
            validIdxs(bi) = sum(numCombs{bi}.*ICPRrels{bi}) == 0;
        end
        ICPRfactors = [ICPRfactors; ICPRrels(validIdxs)]; %#ok
        harmonicsIndexes = [harmonicsIndexes; combinations(validIdxs)]; %#ok
        harmonicsNumbers = [harmonicsNumbers; numCombs(validIdxs)]; %#ok
    end
end