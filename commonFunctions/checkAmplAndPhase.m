function Result = checkAmplAndPhase(mySignalParams, mySignal, mode, fileNm)
    %Out initial amplitude and phase signals and processed signal's params.
    %Compare amplitudes and phases.
    if ~exist('mode', 'var'), mode = ''; end; if ~exist('fileNm', 'var'), fileNm = ''; end
    if isempty(mode), mode = 'unwrap, defreq, normze'; end
    strs = {'-', ';', '-.', '--'}; strs = [strs strs strs]; %Line types.
    myFs = mySignal.getFs;
    t = mySignal.getTimeVector;
    if isstruct(mySignalParams)
        phaseIni = (mySignalParams.phi);
        amplIni = mySignalParams.A;
    else
        phaseIni = compFullPhase(mySignalParams, mode);
        phaseIni = signalObj(mySignalParams.getFs, phaseIni);
        if nnz(strfind(mode, 'normze'))
            mySignal = normSgnl(mySignal);
            mySignalParams = normSgnl(mySignalParams);
        end
        Result.sgnIni = double(mySignalParams);
        Result.sgnFltd = double(mySignal);
        Result.t = t;
        amplIni = getAmplitude(mySignalParams);
        amplIni = signalObj(mySignalParams.getFs, amplIni);
        figure('units', 'points', 'Position', [0 ,0 ,400,300], 'Visible', 'on');
        if strfind(mode, 'plotogether'), subplot(2, 1, 1); end
        plot(t, double(mySignalParams)); hold on; plot(t, double(mySignal)); %, 'k-.'
        xlabel('Время, сек'); ylabel('Амплитуда, В');
        sc = max(abs([Result.sgnIni Result.sgnFltd])); ylim([-sc sc]);
        xlim([0 5]);
        legend('Исходный', 'Выделенный'); title('Реализации исходного и выделенного сигналов');
        Result.signDiff = double(mySignalParams)-double(mySignal); Result.signSE = Result.signDiff.^2;
        if ~strfind(mode, 'plotogether')
            if ~isempty(fileNm), saveas(gcf, [fileNm '_signal.jpg'], 'jpg'); end
            figure('units', 'points', 'Position', [0 ,0 ,400,300], 'Visible', 'on');
            plot(t, Result.signDiff); hold on; plot(t, Result.signSE);
        else
            subplot(2, 1, 2); hold on;
            plot(t, Result.signDiff); plot(t, Result.signSE); %, 'k-.'
        end
        xlim([0 5]);
        xlabel('Время, сек'); ylabel('Амплитуда, В');
        sc = max(abs([Result.signDiff Result.signSE])); ylim([-sc sc]);
        legend({'Ошибка', 'Квадратичная ошибка'});
        title('Ошибка и квадратичная ошибка');
        Result.signMSE = sqrt(mean(Result.signSE)); Result.signME = mean(abs(Result.signDiff));
        sgn = [Result.signDiff Result.signSE]; sgn = sgn( floor(numel(sgn)*0.1):floor(numel(sgn)*0.9) );
        text(t( floor(numel(t)*0.1) ), 0.9*max(sgn), ['MSE: ' num2str(Result.signMSE) ...
            '  ME: ' num2str(Result.signME)], 'Color', 'red'); ylim([min(sgn)*1.05 max(sgn)*1.1]);
        strNm = '_signDiff.jpg'; if strfind(mode, 'plotogether'), strNm = '_signAndErr'; end
        if ~isempty(fileNm), saveas(gcf, [fileNm strNm], 'jpg'); end
    end
    figure('units', 'points', 'Position', [0 ,0 ,400,300], 'Visible', 'on');
    if strfind(mode, 'plotogether')
        subplot(2, 1, 1)
    end
    plot(t, double(phaseIni));
    if nnz(strfind(mode, 'sepplot'))
        title('Initial phase');
        if ~isempty(fileNm), saveas(gcf, [fileNm '_phsIni.jpg'], 'jpg'); end
        figure('units', 'points', 'Position', [0 ,0 ,800,600], 'Visible', 'on');
        hold on; title('Detrended full phase of real signal.');
        myFullPhase = compFullPhase(mySignal, mode); %plot(t, myFullPhase);
        if ~isempty(fileNm), saveas(gcf, [fileNm '_phsFltd.jpg'], 'jpg'); end
    else
        hold on; title('Clean and filtrated signals phases comparison');
        myFullPhase = compFullPhase(mySignal, [mode 'gcf']);
        legend({'Initial phase', 'Fltrd phase'});
        legend('Исходный', 'Выделенный'); title('Девиации фазы исходного и выделенного сигналов');
        xlabel('Время, сек'); ylabel('Фаза, рад');
        sgph = [myFullPhase( floor(0.1*numel(myFullPhase)):floor(0.9*numel(myFullPhase)) ) double(phaseIni)];
        ylim([min(sgph), max(sgph)]); %From minimum value to maximum value within assigned range.
%         if ~isempty(fileNm), saveas(gcf, [fileNm '_phases.jpg'], 'jpg'); end
    end
    Result.phaseDiff = double(phaseIni)-myFullPhase; Result.phaseSE = Result.phaseDiff.^2;
    if strfind(mode, 'plotogether')
        subplot(2, 1, 2)
    else
        figure('units', 'points', 'Position', [0 ,0 ,800,600], 'Visible', 'on');
    end
    plot(t, Result.phaseDiff); hold on; plot(t, Result.phaseSE);
    legend({'Difference', 'Squired error'});
    title('Difference between initial phase and phase of real signal.');
        legend({'Ошибка', 'Квадратичная ошибка'});
        title('Ошибка и квадратичная ошибка фазы');
        xlabel('Время, сек'); ylabel('Фаза, рад');
    Result.phaseMSE = sqrt(mean(Result.phaseSE)); Result.phaseME = mean(abs(Result.phaseDiff));
    text(t( floor(numel(t)*0.1) ), 1.1*max([Result.phaseSE Result.phaseDiff]), ['MSE: ' num2str(Result.phaseMSE) ...
            '  ME: ' num2str(Result.phaseME)], 'Color', 'red');
    if ~isempty(fileNm), saveas(gcf, [fileNm '_phsDiff.jpg'], 'jpg'); end
    
    disp('Mean and max phase difference:');
    disp(mean(phaseIni.getSignalVector-myFullPhase))
    disp(max(phaseIni.getSignalVector-myFullPhase))
    %%%===
    myFullAmpl = getAmplitude(mySignal);
%     if nnz(strfind(mode, 'normze'))
%         myFullAmpl = myFullAmpl/max(myFullAmpl);
%         amplIni = normSgnl(amplIni);
%     end
    figure('units', 'points', 'Position', [0 ,0 ,800,600], 'Visible', 'on');
    plot(t, double(amplIni));
    if nnz(strfind(mode, 'sepplot'))
        title('Initial ampl');
        if ~isempty(fileNm), saveas(gcf, [fileNm '_amplIni.jpg'], 'jpg'); end
        figure('units', 'points', 'Position', [0 ,0 ,800,600], 'Visible', 'on');
        hold on; title('Amplitude of real signal.');
    else
        hold on; title('Clean and filtrated signals amplitude comparison');
    end
    plot(t, myFullAmpl);
    if ~nnz(strfind(mode, 'sepplot'))
        legend({'Initial ampl', 'Fltrd ampl'}); 
        if ~isempty(fileNm), saveas(gcf, [fileNm '_ampls.jpg'], 'jpg'); end
    else
        if ~isempty(fileNm), saveas(gcf, [fileNm '_amplFltd.jpg'], 'jpg'); end
    end
    Result.amplDiff = double(amplIni)-myFullAmpl; Result.amplSE = Result.amplDiff.^2;
    figure('units', 'points', 'Position', [0 ,0 ,800,600], 'Visible', 'on');
    plot(t, Result.amplDiff); hold on; plot(t, Result.amplSE);
    legend({'Difference', 'Squired error'});
    title('Difference between initial amplitude and amplitude of real signal.');
    Result.amplMSE = sqrt(mean(Result.amplSE)); Result.amplME = mean(abs(Result.amplDiff));
    text(t(end)/2, 0.9*max([Result.amplSE Result.amplDiff]), ['MSE: ' num2str(Result.amplMSE) ...
            '  ME: ' num2str(Result.amplME)], 'Color', 'red');
    if ~isempty(fileNm), saveas(gcf, [fileNm '_amplDiff.jpg'], 'jpg'); end
    disp('Mean and max amplitude difference:');
    disp(mean(amplIni.getSignalVector-myFullAmpl))
    disp(max(amplIni.getSignalVector-myFullAmpl))
end