classdef phaseProcessing<polyHarmonic

	methods(Access = public)
		
		function myPhaseProcessing = phaseProcessing(varargin)
			myPhaseProcessing = myPhaseProcessing@polyHarmonic(varargin{:});
        end
	
		function [ICPRmeans, ICPRstds, phasesSTDs, ICPRfactors, harmonicsIndexes, ICPRs, unwrappedPhs, plottingPhases] = ICPRprocessing(myPhaseProcessing)
            relationsError = 5; %percents - allowed error in central freqz relations.
            if nargin<5, stdDisp = true; end
            lineStyles = {'k', 'k-.', 'k--'};
            
            %Obtain full phases; cut off possible runouts at the ends of realisation.
            unwrappedPhs = compFullPhase(myPhaseProcessing, 'unwrap, setnan5%');  %, setnan5%
            unwrappedPhs = cellfun(@(x) reshape(x, [], 1), unwrappedPhs, 'UniformOutput', false);
			%Get phase distortions.
            plottingPhases = compFullPhase(myPhaseProcessing, 'unwrap, defreq'); %For plotting
            plottingPhases = cellfun(@(x, y) x/y, plottingPhases, num2cell(1:numel(plottingPhases)), 'UniformOutput', false);
            
			centrFreqz = arrayfun(@(x) x.getCentralFreq, myPhaseProcessing.components);
            centrFreqz = centrFreqz/centrFreqz(1); %Get relations
            harmoNums = round(centrFreqz);
            %Delete harmonics with fraction relations.
            delIdxs = centrFreqz - harmoNums > relationsError/100;
            harmoNums(delIdxs) = [];
			[ICPRfactors, ~, harmonicsIndexes] = ICPRcombinations(harmoNums, 3);
            
            %Get ICPRs.
            ICPRs = cell(0, 0);
            for ai = 1:size(ICPRfactors, 1)
                relatedPhases = cellfun(@(x, y) x*y, unwrappedPhs(harmonicsIndexes{ai}), num2cell(ICPRfactors{ai}), 'UniformOutput', false);
                ICPRs(end+1) = {reshape(sum(horzcat(relatedPhases{:}), 2), size(myPhaseProcessing.t))}; %#ok
            end
            %Translate ICPRs to unambuguous range [-pi; pi)
            ICPRs = cellfun(@(x) x - repmat(floor(median(x, 'omitnan')/pi)*pi, size(x)), ICPRs, 'UniformOutput', false);
            ICPRmeans = cellfun(@(x) mean(x, 'omitnan'), ICPRs);
            ICPRstds = cellfun(@(x) std(x, 'omitnan'), ICPRs);
            phasesSTDs = cellfun(@(x) std(x), plottingPhases);
			
		end
		
    end

end