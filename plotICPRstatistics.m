close all
mkdir('Out');

% dataset = {'EDP_7_norm' 'EDP_7_8_9_horz_parall'};
% baseFreq = 'shaftFreq';
% domain = 'velocity';
% icprLabel = 'PQI12';
% ICPRindex = 1;

if iscell(dataset)
    files = [];
    for ai = 1:numel(dataset)
        files = [files; addDataset(fullfile(pwd, 'In', dataset{ai}))];
    end
else
    files = addDataset(fullfile(pwd, 'In', dataset));
end

%Plot whole dataset.
allFiles = plotState(files, domain, baseFreq, icprLabel, ICPRindex, 'all states');
conditions = loadConditions(files);

%Vertical misalignment dataset.
if nnz(strcmp('Misalignment_vert_parall', dataset))
    [indexesMin, conditionLabelMin] = getValidIndexes(conditions, {'misalignmentDegree'}, {'0.34'}); %0.34 mm state.
    [indexesMid, conditionLabelMid] = getValidIndexes(conditions, {'misalignmentDegree'}, {'0.8'}); %0.8 mm state.
    [indexesMax, conditionLabelMax] = getValidIndexes(conditions, {'misalignmentDegree'}, {'1.09'}); %1.09 mm state.
    data = plotStates(files, domain, baseFreq, icprLabel, ICPRindex, {indexesMin, indexesMid, indexesMax}, {conditionLabelMin, conditionLabelMid, conditionLabelMax});
    save(fullfile('Out', ['Misalignment_vert_parall - ' domain ' - ' baseFreq ' - ' icprLabel '.mat']), 'data');
end

%Horizontal misalignment dataset.
if nnz(strcmp('Misalignment_horz_parall', dataset))
    [indexesN, conditionLabelN] = getValidIndexes(conditions, {'misalignmentDegree'}, {'0'}); %Normal state.
    [indexesMin, conditionLabelMin] = getValidIndexes(conditions, {'misalignmentDegree'}, {'0.5'}); %0.5 mm state.
    [indexesMax, conditionLabelMax] = getValidIndexes(conditions, {'misalignmentDegree'}, {'1.15'}); %1.15 mm state.
    data = plotStates(files, domain, baseFreq, icprLabel, ICPRindex, {indexesN, indexesMin, indexesMax}, {conditionLabelN, conditionLabelMin, conditionLabelMax});
    save(fullfile('Out', ['Misalignment_horz_parall - ' domain ' - ' baseFreq ' - ' icprLabel '.mat']), 'data');
end

%MPZ bearings dataset.
if nnz(strcmp('MPZ', dataset))
%     [indexesB, conditionLabelB] = getValidIndexes(conditions, {'state'}, {'burnish'});
%     [indexesN, conditionLabelN] = getValidIndexes(conditions, {'state'}, {'normal'});
%     [indexesD, conditionLabelD] = getValidIndexes(conditions, {'state'}, {'defect'});
%     data = plotStates(files, domain, baseFreq, icprLabel, ICPRindex, {indexesB, indexesN, indexesD}, {conditionLabelB, conditionLabelN, conditionLabelD});
%     save(fullfile('Out', ['Bearings - ' domain ' - ' baseFreq ' - ' icprLabel '.mat']), 'data');
    %MPZ bearings dataset - 2nd bearing.
    [indexesB, conditionLabelB] = getValidIndexes(conditions, {'state', 'bearing'}, {'burnish', {'2', '3'} });
    [indexesN, conditionLabelN] = getValidIndexes(conditions, {'state', 'bearing'}, {'normal', {'2', '3'} });
    [indexesD, conditionLabelD] = getValidIndexes(conditions, {'state', 'bearing'}, {'defect', {'2', '3'} });
    data = plotStates(files, domain, baseFreq, icprLabel, ICPRindex, {indexesB, indexesN, indexesD}, {conditionLabelB, conditionLabelN, conditionLabelD});
    save(fullfile('Out', ['Bearing 2 - ' domain ' - ' baseFreq ' - ' icprLabel '.mat']), 'data');
end


function data = plotStates(files, domain, baseFreq, icprLabel, ICPRindex, indexes, conditionLabels)
    plotMarkers = {'go', 'ks', 'r*', 'c+'};
    myFig = figure('Units', 'points', 'Position', [0, 0, 800, 600], 'Visible', 'on', 'Color', 'w');
    myAx = axes ('parent', myFig);
    myFigHist = figure('Units', 'points', 'Position', [0, 0, 800, 600], 'Visible', 'on', 'Color', 'w');
    myAxHist = axes ('parent', myFigHist);
    data = cell(size(conditionLabels));
    succ = false(size(conditionLabels)); %Successful plotting.
    for ai = 1:numel(conditionLabels)
        [data{ai}, succ(ai)] = plotState(files(indexes{ai}), domain, baseFreq, icprLabel, ICPRindex, conditionLabels{ai}, myAx, plotMarkers{ai});
    end
    legend(myAx, conditionLabels(succ));
    label = [domain ' - ' baseFreq ' - ' icprLabel ' - Histogram'];
    plotHist([data{succ}], label, myAxHist)
    legend(myAxHist, conditionLabels(succ));
end

function [data, succ] = plotState(files, domain, baseFreq, icprLabel, ICPRindex, conditionLabel, myAx, plotMarkers)
    if nargin < 8
        plotMarkers = [];
    end
    if nargin < 7
        myAx = [];
    end
%     load(fullfile(files(1).folder, files(1).name), 'myFile');
%     harmonicsNumbers = myFile.(domain).(baseFreq).harmonicsIndexes{ICPRindex};
    harmonicsNumbers = {[1 2], [1 3], [], [], [], [], [1 2 3]};
    harmonicsNumbers = harmonicsNumbers{ICPRindex};
    disp(harmonicsNumbers);
    peaks = getFromFile(files, domain, baseFreq, 'logSNRs');
    SNRs = peaks(:, harmonicsNumbers);
    SNRs = min(SNRs, [], 2);
    ICPRstds = getFromFile(files, domain, baseFreq, 'ICPRstds');
    label = [domain ' - ' baseFreq ' - ' icprLabel]; % conditionLabel ' - '
    succ = plotICPRdistr(SNRs, ICPRstds(:, ICPRindex), label, myAx, plotMarkers);
    data.SNRs = SNRs;
    data.ICPRstds = ICPRstds(:, ICPRindex);
    data.names = {files.name};
    data.conditionLabel = conditionLabel;
end

function result = getFromFile(files, domain, baseFreq, field)
    %load(fullfile(files(1).folder, files(1).name), 'myFile');
    vectNum = 8; %numel(myFile.(domain).(baseFreq).(field));
    result = NaN(numel(files), vectNum);
    for ai = 1:numel(files)
        load(fullfile(files(ai).folder, files(ai).name), 'myFile');
        bi = 1:numel(myFile.(domain).(baseFreq).(field));
        result(ai, bi) = myFile.(domain).(baseFreq).(field)(bi);
    end
end

function conditions = loadConditions(files)
    conditions = cell(size(files));
    for ai = 1:numel(files)
        load(fullfile(files(ai).folder, files(ai).name), 'myFile');
        conditions{ai} = myFile.conditions;
    end
    conditions = [conditions{:}];
end

function [indexes, conditionLabel] = getValidIndexes(conditions, conditionNames, conditionValues)
    conditionLabel = '';
    indexes = true(size(conditions));
    for ai = 1:numel(conditionNames)
        currCond = arrayfun(@(x) x.(conditionNames{ai}), conditions, 'UniformOutput', false);
        desiredValues = conditionValues{ai};
        if ~iscell(desiredValues)
            desiredValues = {desiredValues}; end
        currentIndexes = false(size(indexes));
        for bi = 1:numel(desiredValues)
            currentIndexes = currentIndexes | cellfun(@(x) strcmp(x, desiredValues{bi}), currCond);
        end
        indexes = indexes & currentIndexes;
        conditionLabel = [conditionLabel ' - ' conditionNames{ai} ' ' strjoin(desiredValues, ',')];
    end
end



function succ = plotICPRdistr(SNRs, ICPRstds, label, myAx, plotMarkers, magnif)
succ = false;
    if ~nnz(~isnan(ICPRstds))
        warning(['No valid samples at ' label]);
        return;
    end
    if nargin < 5 || isempty(myAx)
        myFigure = figure('Units', 'points', 'Position', [0, 0, 800, 600], 'Visible', 'off', 'Color', 'w');
        myAx = axes ('parent', myFigure);
    else
        myFigure = myAx.get('Parent');
    end
    if nargin < 5 || isempty(plotMarkers)
        plotMarkers = 'ro';
    end
    if nargin < 6
        magnif = 1;
    end
    hold on;
    plot(myAx, SNRs', ICPRstds', plotMarkers);
    X = linspace(min(SNRs), max(SNRs));
    %plot(X, repmat(0.1, size(X)));
    xlabel(myAx, 'SNR');
    ylabel(myAx, 'ICPR STD value, rad');
    title(myAx, label);
    succ = true;
    saveas(myFigure, fullfile(pwd, 'Out', [label '.jpg']), 'jpg');
    saveas(myFigure, fullfile(pwd, 'Out', [label '.fig']), 'fig');
    %Consider outliers.
    radianHold = 10;
    if max(ICPRstds) > radianHold && magnif
        currPlot = ICPRstds;
        currPlot = currPlot(currPlot < radianHold);
        if nnz(max(currPlot))
            ylim(myAx, [0, max(currPlot)*1.05]);
            saveas(myFigure, fullfile(pwd, 'Out', [label ' __Magnified__ .jpg']), 'jpg');
            saveas(myFigure, fullfile(pwd, 'Out', [label ' __Magnified__ .fig']), 'fig');
        end
    end
end

function plotHist(states, label, myAx, legendoz)
    if nargin < 3 || isempty(myAx)
        myFigure = figure('Units', 'points', 'Position', [0, 0, 800, 600], 'Visible', 'on', 'Color', 'w');
        myAx = axes ('parent', myFigure);
    else
        myFigure = myAx.get('Parent');
    end
    hold on
    maxBin = 0;
    for ai = 1:numel(states)
        myState = sort(states(ai).ICPRstds);
        F = cumsum(myState);
        F = F/max(F);
        %Find 90% levels of probabilitu distribution of the current state.
        index = find(F>0.9, 1, 'first');
        maxBin = max(maxBin, myState(index));
    end
    binRange = linspace(0, maxBin, 10);
    hcs = cell(size(states));
    for ai = 1:numel(states)
        hcx = histcounts(states(ai).ICPRstds, [binRange Inf]);
        hcs{ai}=hcx./sum(hcx);
    end
    bar(myAx, binRange, vertcat(hcs{:})');
    hold on; grid;
    xlabel('ICPR STD, rad');
    ylabel('Probability density function');
    title(label);
    if nargin >=4
        legend(legendoz);
    end
    saveas(myFigure, fullfile(pwd, 'Out', [label '.jpg']), 'jpg');
    saveas(myFigure, fullfile(pwd, 'Out', [label '.fig']), 'fig');
end

function files = addDataset(Root)
    folders = genpath(Root);
    folders = strsplit(folders, ';');
    folders = folders(~cellfun(@(x) isempty(x), folders));
    files = cell(size(folders));
    for ai = 1:numel(files)
        files{ai} = dir(fullfile(folders{ai}, '*.mat'));
    end
    files = vertcat(files{:});
    currentDataSet = strsplit(Root, filesep);
    files = arrayfun(@(x) setfield(x, 'dataset', currentDataSet{end}), files);
end
