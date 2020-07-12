function myFile = processFileStruct(myFile)
    %Put here desirable domains.
    spectraFields = {'acceleration', 'velocity', 'displacement'};
    bandWidth = myFile.bandWidth;
    %Specify informative frequencies or leave empty array [] 2 process all frequencies.
    baseFrequencies = {'shaftFreq'}; %'FTF', 'BSF'.

    basicFreqs = myFile.baseFrequencies;
    for ai = 1:length(spectraFields)
        mySignalObj = signalObj(myFile.Fs, myFile.(spectraFields{1, ai}).signal);
            for bi = 1:length(basicFreqs)
                if ~isempty(baseFrequencies) && ~nnz(ismember(basicFreqs{bi}, baseFrequencies))
                    continue;
                end
                currentBase = myFile.(spectraFields{1, ai}).(basicFreqs{bi});
                freqVect = currentBase.frequencies;
                if isempty(freqVect)
                    warning(['Frequency ' basicFreqs{bi} ' is not assigned!']);
                    continue;
                end
                myPhaseProcessing = phaseProcessing(myFile.Fs, mySignalObj, 'time', mySignalObj.getTimeVector);
                myBand = arrayfun(@(x) [x-bandWidth/2 x+bandWidth/2], freqVect, 'UniformOutput', false);
                myPhaseProcessing = setBand(myPhaseProcessing, myBand, 0, {'fft'});
                myPhaseProcessing = myPhaseProcessing.filtSign; %Get components at the current base frequency harmonics.
                [currentBase.ICPRmeans, currentBase.ICPRstds, currentBase.phasesSTDs, currentBase.ICPRfactors, ...
                    currentBase.harmonicsIndexes, currentBase.ICPRs, currentBase.fullPhases, currentBase.phaseDistortions] = ICPRprocessing(myPhaseProcessing);
                myFile.(spectraFields{1, ai}).(basicFreqs{bi}) = currentBase;
            end
    end
    myFile.baseFrequencies = baseFrequencies;
end