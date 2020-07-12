function plotICPRsAndSpec(files, domain, baseFreq, icprLabels, ICPRindexes)
stdDisp = 1;
    dataset = files(1).dataset;
    myFolder = fullfile(pwd, 'Out', dataset);
    mkdir(myFolder);
    for ai = 1:numel(files)
        load(fullfile(files(ai).folder, files(ai).name), 'myFile');
        ICPRs = vertcat(myFile.(domain).(baseFreq).ICPRs{ICPRindexes});
        dt = 1/myFile.Fs;
        t = 0:dt:length(ICPRs)*dt-dt;
        myFig = figure('Units', 'points', 'Position', [0, 0, 800, 600], 'Visible', 'off', 'Color', 'w');
        hold on; plot(t, ICPRs);
        label = [domain ' - ' baseFreq strjoin(icprLabels, ',')];
        title(label)
        numbers = arrayfun(@(x) num2str(myFile.(domain).(baseFreq).harmonicsIndexes{x}), ICPRindexes, 'UniformOutput', false);
        numbers = cellfun(@(x) regexprep(x, ' +', ','), numbers, 'UniformOutput', false);
        legendos = cellfun(@(x) ['\Phi_{' x '}'], numbers, 'UniformOutput', false);
        if stdDisp
            legendos = arrayfun(@(x, y) sprintf('%s - STD %1.3f', x{1}, y), legendos, ...
                myFile.(domain).(baseFreq).ICPRstds(ICPRindexes), 'UniformOutput', false);
        end
        legend(legendos, 'Location', 'best');
        saveas(myFig, fullfile(myFolder, [files(ai).name ' - ' label '.jpg']), 'jpg');
        saveas(myFig, fullfile(myFolder, [files(ai).name ' - ' label '.fig']), 'fig');
        
        myFig = figure('Units', 'points', 'Position', [0, 0, 800, 600], 'Visible', 'on', 'Color', 'w');
        hold on;
        title(label)
        plot(myFile.frequencyVector, myFile.(domain).spectrum);
        stem(myFile.(domain).(baseFreq).frequencies, myFile.(domain).(baseFreq).linAmplitudes)
        legend('Spectrum', 'Informative frequencies', 'Location', 'best');
        xlim([min(myFile.(domain).(baseFreq).frequencies)*0.85 max(myFile.(domain).(baseFreq).frequencies)*1.15]);
        saveas(myFig, fullfile(myFolder, [files(ai).name ' - spectrum - ' label '.jpg']), 'jpg');
        saveas(myFig, fullfile(myFolder, [files(ai).name ' - spectrum - ' label '.fig']), 'fig');
        
        close all
    end
    
end