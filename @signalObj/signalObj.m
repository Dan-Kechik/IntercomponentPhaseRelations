classdef signalObj
    
    properties(SetAccess = protected, GetAccess = public)
        signal
        signalIni
        t
        Fs
		band
		bandFilter
        Label
    end
    
    methods(Access = public)   
        
        function mySignalObj = signalObj(myFs, mySignal, myTime, mode)
            if ~exist('mySignal', 'var'), mySignal = []; end
            if ~exist('myTime', 'var'), myTime = []; end
            if ~exist('mode', 'var'), mode = []; end
            mySignalObj.Fs = myFs;
            if isempty(mySignal), mySignal = zeros(1, floor(myFs/2)); end
            mySignalObj.signal = reshape(double(mySignal), 1, []);
            %Compute a time vector according to setted mode or a signal sampling and Fs.
            mySignalObj = compTimeVector(mySignalObj, myTime, mode);
            mySignalObj.band = {};
            mySignalObj.bandFilter = {};
        end

        function myTime = getTimeVector(mySignalObj)
            myTime = mySignalObj.t;
        end   

        function myFs = getFs(mySignalObj)
            myFs = mySignalObj.Fs;
        end   

        function mySignalObj = setTimeVector(mySignalObj, myTime)
            mySignalObj.t = myTime;
        end

        function mySignal = getSignalVector(mySignalObj)
            mySignal = mySignalObj.signal;
        end
    
        function mySignalObj = setSignalVector(mySignalObj, mySignal)
            mySignalObj.signal = double(mySignal);
        end

        function mySignalObj = setBand(mySignalObj, myBand, reset, params)
		%Band is a cell array, assign passbands 4 filters.
		%Set/add a band, sort them, build a filters, compute a signal components.
			if ~exist('params', 'var'), params = {''}; end
			if ~exist('reset', 'var'), reset = 0; end
			if ~reset
				myBand = [myBand mySignalObj.band];
				startPos = numel(mySignalObj.band) + 1;  %Process only added elements.
			else
				startPos = 1;  %Process setted elements.
            end
            if numel(params) == 1, params = repmat(params, size( startPos:numel(myBand) )); end
            mySignalObj.band = myBand;
			for i = startPos:numel(myBand)
                if isa(params{i}, 'digitalFilter'), mySignalObj.bandFilter(i) = params(i); continue; end
                if iscell(params{i}) %Heap of arguments is defined - search for filter type label.
                    strings = params{i}( cellfun(@(x) ischar(x), params{i}) );
                    k = strfind(strings, 'filter', 'ForceCellOutput', true);
                else %Just only char array denoting filter type.
                    k = strfind(params, 'filter', 'ForceCellOutput', true);
                end
                if ~nnz(cellfun(@(x) nnz(x), k)) %Set filtering mode if it's unnecessary ...
                    mySignalObj.bandFilter(i) = params(i); continue; %to design filter.
                end
                mySignalObj.bandFilter{i} = eval(['designfilt(' params ');']);
			end
			mySignalObj = sortBand(mySignalObj);
% 			mySignalObj = mySignalObj.filtSign;
        end

        function myBand = getBand(mySignalObj)
            myBand = mySignalObj.band;
        end
        
        function mySignalObj = setLabel(mySignalObj, myLabel)
            mySignalObj.Label = myLabel;
        end
	
        function myFullPhase = compFullPhase(mySignalObj, mode)
			if ~exist('mode', 'var'), mode = 0; end %On/off unwrapping or set up mode of phase computation.
            need2unwrap=0; if isnumeric(mode), need2unwrap=mode; mode=''; end
            if strfind(mode, 'aver')
                myFullPhase = compFullPhase(mySignalObj, 'unwrap, defreq, whole');
                myFullPhase = mean(myFullPhase, 'omitnan');
                if strfind(mode, 'degr'), myFullPhase = myFullPhase*180/pi; end
                if strfind(mode, 'pi_num'), myFullPhase = myFullPhase/pi; end
                return;
            end
            if strfind(mode, 'unwrap'), need2unwrap = 1; end
            mySignal = mySignalObj.getSignalVector;
            if ~nnz(strfind(mode, 'eyeDiag'))
                analytSign = hilbert(mySignal);
                rge = regexp(mode, 'setnan.+%', 'match');
                if numel(rge)
                    num = strrep(strrep(rge{1}, 'setnan', ''), '%', '');
                    num = str2double(num)/100;
                    [val, sam1] = min(abs( mySignalObj.t-mySignalObj.t(end)*num ));
                    [val, sam2] = min(abs( mySignalObj.t-mySignalObj.t(end)*(1-num) ));
                    analytSign([1:sam1, sam2:end]) = NaN(size( analytSign([1:sam1, sam2:end]) ));
                end
                myFullPhase = angle(analytSign);
                if need2unwrap
                    myFullPhase = unwrap(myFullPhase);
                end
            else
                sgnS = sign(mySignal);
                phase = NaN(size(mySignal));
                %Fill in half-period values by crossing zero.
                phaseIterator = 0;
                deriv = [0 diff(sgnS)];
                sampling = find(deriv);
                for ai = sampling
                    phase(ai) = phaseIterator;
                    phaseIterator = phaseIterator + pi;
                end
                myFullPhase = inpaint_nans(phase);
            end
            if strfind(mode, 'detrend'), myFullPhase = detrend(myFullPhase, 'linear', [1 numel(myFullPhase)]); end
            if strfind(mode, 'defreq')
                centralFreq = mySignalObj.getCentralFreq;
                centralFreqComp = 2*pi*centralFreq*mySignalObj.t;
                myFullPhase = myFullPhase - centralFreqComp;
            end
            if nnz(strfind(mode, 'linear'))
                cutRng = ceil(0.05*numel(myFullPhase));
                leftB = myFullPhase(1+cutRng); rightB = myFullPhase(end-cutRng);
                step = (rightB-leftB)/numel(myFullPhase);
                submitTrnd = leftB:step:rightB-step; submitTrnd = submitTrnd(1:numel(myFullPhase))-leftB;
                myFullPhase = myFullPhase-reshape(submitTrnd, size(myFullPhase));
            end
            if nnz(strfind(mode, 'median'))
                myFullPhase = reshape(myFullPhase, [], 1);
                myFullPhase = myFullPhase - repmat(floor( median(myFullPhase, 1, 'omitnan')/pi )*pi, size(myFullPhase, 1), 1);
            end
            if strfind(mode, 'deriv')
                dt = mySignalObj.t(2) - mySignalObj.t(1);
                myFullPhase = ([0 diff(myFullPhase)]/dt)/(2*pi);
            end
            if strfind(mode, 'plot')
                if ~nnz(strfind(mode, 'gcf')), figure('Units', 'points', 'Position', [0 0 800 600], 'Visible', 'on', 'Color', 'w'); end; hold on
                if strfind(mode, 'degr'), myFullPhase = myFullPhase*180/pi; end
                plot(mySignalObj.t, myFullPhase);
                if strfind(mode, 'ylim'), ylim([min([myFullPhase -pi]), max([myFullPhase pi])]); end
                %title(sprintf( 'Full phase. Average value is %1.3f.', compFullPhase(mySignalObj, 'aver') ));
                if ~nnz(strfind(mode, 'gcf'))
                    title(sprintf( 'Full phase. Average value is %1.3f.', compFullPhase(mySignalObj, 'aver') ));
                else
                    l = legend('show'); lStr = l.String;
                    if strcmp(lStr{end}, 'data1'), lStr = []; end %There was no any legend, creqte a new one.
                    newLeg = sprintf('Phs #%d, aver %1.3f, std %1.3f', numel(lStr)+1, compFullPhase(mySignalObj, 'aver, degr'), std(myFullPhase*180/pi, 'omitnan'));
                    legend([lStr, newLeg], 'Location', 'best');
                end
                set(gcf, 'Name', 'Full phase. Average values are displayed');
            end
        end
	
        function myCentralFreq = getCentralFreq(myHarmonic, myFullPhase)
            %Averaging of full phases derivative.
            if nargin<2
                myFullPhase = compFullPhase(myHarmonic, 'unwrap, whole'); end
            dt = myHarmonic.t(2) - myHarmonic.t(1);
            deriv = diff(myFullPhase)/dt;
            myCentralFreq = mean(deriv, 'omitnan')/(2*pi);
        end
	
        function myEnvel = getAmplitude(mySignalObj)
            mySignal = mySignalObj.getSignalVector;
            myEnvel = envelope(mySignal);
        end

        function mySignalObj = compTimeVector(mySignalObj, myTime, mode)
            %Set time vector from Fs and pointed length (in seconds (def) or in samples), if time is a scalar, or pointed time vector.
            if numel(myTime) > 1
                mySignalObj.t = myTime;
                return;
            end
            if ~exist('mode', 'var')
                mode = [];
            end
            if ~exist('myTime', 'var')
                myTime = [];
            end
            if isempty(mode)
                mode = 'seconds';
            end
            if isempty(myTime)  %Set time vector equal to signal length.
                mode = 'signalSampling';
                myTime = numel(mySignalObj.signal);
            end
            dt = 1/mySignalObj.Fs;
            switch mode
                case 'signalSampling'
                    mySignalObj = compTimeVector(mySignalObj, myTime, 'samples');
                    return;
                case 'samples'
                    %Translate to seconds.
                    Tmax = myTime*dt;
                case 'seconds'
                    Tmax = myTime;
                otherwise
                    error('Unknown time setting mode %s.', mode);
            end
            mySignalObj.t = 0:dt:Tmax-dt;  %
        end
		
		function mySignalObj = compSignalVector(mySignalObj)
            mySignal = mySignalObj.signal;
            len = length(mySignalObj.t);
            repeationNumber = floor(len/length(mySignal));
            mySignalResampled = repmat(mySignal, 1, repeationNumber);
            %The rest part of period.
            restSamples = len - repeationNumber*length(mySignal);
            mySignalResampled = [mySignalResampled mySignal(1:restSamples)];
            mySignalObj.signal = mySignalResampled;
        end
        
        function plotSignal(mySignalObj, varargin)
            myTime = mySignalObj.t;
            mySignal = mySignalObj.signal;
            if numel(myTime) ~= numel(mySignal)
                dbstop error
                error('Samples number mismatch.');
            end
            if numel(varargin)
                if ~nnz(strfind([varargin{:}], 'gcf')), figure('Units', 'points', 'Position', [0 0 800 600], 'Visible', 'on', 'Color', 'w'); end
                [mySpectrum, f] = getSignalSpectrum(mySignalObj, varargin{:});
                if numel(strrep([varargin{:}], 'gcf', ''))
                    subplot(2,1,1), plot(mySignalObj.t, mySignalObj.signal);
                    subplot(2,1,2), plot(f, mySpectrum);
                else, hold on; plot(mySignalObj.t, mySignalObj.signal);
                end
            else
                figure; plot(myTime, mySignal);
            end
        end

        function [mySignalObj, anSg] = specWin(mySignalObj, varargin)
            %Filtrate signal through ifft of band limited fft coeffs.
            %Arguments: filtration band; window arguments (optional).
            myBand = varargin{1};
            [Spec, f] = getSignalSpectrum(mySignalObj, 'complex full');
            [~, lowS] = min(abs( f - myBand(1) )); f1 = f(lowS);
            [~, highS] = min(abs( f - myBand(2) )); f2 = f(highS);
            w = window('gausswin', highS-lowS+1)';
            SpN = zeros(size(Spec));
            SpN(lowS:highS) = Spec(lowS:highS).*w;
            SpN = fliplr(SpN);
            anSg = ifft(SpN)*length(Spec);
            mySign = real(anSg);
            mySignalObj.signal = mySign;
        end

        function [mySignalObj, mySign] = spectrumFiltrate(mySignalObj, varargin)
            %Filtrate signal through ifft of band limited fft coeffs.
            %Arguments: filtration band; window arguments (optional).
            myBand = varargin{1};
            [Spec, f] = getSignalSpectrum(mySignalObj, 'complex full');
            [~, fi1] = min(abs( f - myBand(1) )); f1 = f(fi1);
            [~, fi2] = min(abs( f - myBand(2) )); f2 = f(fi2); filtRangeSamp = [fi1; fi2];
            fSp = ones(size(Spec)); fSp2 = fSp;
            fSp(1+(fi1:fi2)) = zeros(size( fSp(1+(fi1:fi2)) ));
            fSp2(end-(fi1:fi2)+1) = zeros(size( fSp2(end-(fi1:fi2))+1 ));
            fSp = ~fSp; fSp2 = ~fSp2;
            spF = zeros(size(Spec));
            %Apply spectral window with assigned parameters.
            params = [varargin repmat({[]}, 1, 3-numel(varargin))];
            if fi1 == 1 %Zero freq - use one-sided window.
                params{3} = 2*(fi2-fi1); end
            if ~numel(params)
                spF(fSp) = Spec(fSp); spF(fSp2) = Spec(fSp2);
            elseif numel(params) < 4
                    [~, spF(fSp)] = windowing(mySignalObj, params{2:3}, Spec(fSp));
                    [~, spF(fSp2)] = windowing(mySignalObj, params{2:3}, Spec(fSp2));
            elseif numel(params) >= 4
                    [~, spF(fSp)] = windowing(mySignalObj, params{2:3}, Spec(fSp), params{4:end});
                    [~, spF(fSp2)] = windowing(mySignalObj, params{2:3}, Spec(fSp2), params{4:end});
            end
            
            mySign = (ifft(spF))*length(fSp);
            mySignalObj.signal = mySign;
        end
        
        function [mySignalObj, myFiltedSign] = getFFTanalyticSignal(mySignalObj, filtRange)
            df = mySignalObj.Fs/length(mySignalObj.signal);
            fVect = 0:df:mySignalObj.Fs-df;
            [mv, filtRange(1)] = min(abs( fVect - filtRange(1) ));
            [mv, filtRange(end)] = min(abs( fVect - filtRange(end) ));
            filtRange = filtRange(1):filtRange(end);
            
            initialSpec = fft(mySignalObj.signal)/length(mySignalObj.signal);
            filteredSpec = zeros(size(initialSpec));
            
            filteredSpec(filtRange) = initialSpec(filtRange);
            myFiltedSign = (ifft(filteredSpec))*length(filteredSpec);
            myFiltedSign = 2*real(myFiltedSign);
            mySignalObj.signal = myFiltedSign;
        end

        function [mySignalObj, myComponents] = filtSign(mySignalObj, myFilt)
			if ~exist('myFilt', 'var')
				myFilt = [];
			end
			if isempty(myFilt)
				myFilt = mySignalObj.bandFilter;
			end
			if isempty(myFilt)
				warning('Your filter is not assigned!');
                myComponents = {};
				return;
			end
			myComponents = cell(size(myFilt));
            mySignal = mySignalObj.signal;
            mySignalProcessed = zeros(size(mySignal));
			for i = 1:numel(myFilt)
                params = [];
                if iscell(myFilt{i}) %Filter type and additional parameters.
                    params = myFilt{i}(2:end); myFilt(i) = myFilt{i}(1); end
                if ~ischar(myFilt{i}), myComponents{i} = filtfilt(myFilt{i}, mySignal); continue; end
                if strcmp(myFilt{i}, 'fft')
                    if ~isempty(params) %Unpack varargin.
                        [~, myComponents{i}] = spectrumFiltrate(mySignalObj, mySignalObj.band{i}, params{:});
                    else
                        [~, myComponents{i}] = spectrumFiltrate(mySignalObj, mySignalObj.band{i});
                    end
                elseif strcmp(myFilt{i}, 'specWin')
                    [~, myComponents{i}] = specWin(mySignalObj, mySignalObj.band{i});
                elseif strcmp(myFilt{i}, 'analyt')
                    [~, myComponents{i}] = getFFTanalyticSignal(mySignalObj, mySignalObj.band{i});
                elseif strcmp(myFilt{i}, 'mod')
                    [~, myComponents{i}] = modFiltrate(mySignalObj, mySignalObj.band{i});
                elseif nnz(strfind(myFilt{i}, 'wavtFFT'))
                    waveletName = strrep(myFilt{i}, 'wavtFFT', '');
                    File = struct('signal', reshape(mySignal, [], 1), 'Fs', mySignalObj.Fs);
                    frequency = mean(mySignalObj.band{i});
                    myComponents{i} = reshape(waveletFilteringAndInstFreq(File, frequency, waveletName), 1, []);
                elseif nnz(strfind(myFilt{i}, 'wavtTF'))
                    waveletName = strrep(myFilt{i}, 'wavtTF', '');
                    myComponents{i} = waveletFiltering(mySignalObj, mySignalObj.band{i}, waveletName);
                else
                    myComponents{i} = filtfilt(myFilt{i}, mySignal);
                end
				mySignalProcessed = mySignalProcessed + myComponents{i};
            end
            mySignalObj.signal = real(mySignalProcessed);
        end

        function [filtedComp, mySign] = modFiltrate(mySignalObj, myBand)
            [~, SSBsignal2] = SSBmod(mySignalObj, myBand(:, 2), mySignalObj.Fs, 'low imag');
            if ~nnz(myBand(:, 1)) %If the first band is zero, just translate signal back.
                [~, filtedComp] = SSBmod(SSBsignal2, myBand(:, 2), SSBsignal2.getFs, 'low imag'); %filtedComp = SSBsignal2;
                mySign = filtedComp.getSignalVector; return; 
            end %Low pass filtration.
            [~, SSBsignal] = SSBmod(mySignalObj, myBand(:, 1), mySignalObj.Fs, 'low imag');
            %Transmit freqs back.
            [~, SSBsignalHi] = SSBmod(SSBsignal, myBand(:, 1), mySignalObj.Fs, 'low imag');
            [~, SSBsignalHi2] = SSBmod(SSBsignal2, myBand(:, 2), mySignalObj.Fs, 'low imag');
            %Submit signal with deleted lower freqs and the signal band from 2 rest only signal band.
            filtedComp = SSBsignalHi2 - SSBsignalHi; %filtedComp = ldivide(SSBsignalHi2 - SSBsignalHi, 2);
            mySign = filtedComp.getSignalVector;
        end
        
        function [SSBsignal, SSBobj] = SSBmod(mySignalObj, carrier, Fs, mode, carrAngle)
            %Getting AM signals. Default - suppressed carrier. Possible simple AM,
            %low/high SSB, suppress image channel or image only. If it' necessary - restore carrier.
            %Carrier is frequency (generate sinus with def params) or vector.
                if ~exist('mode', 'var'), mode = ''; end
                if ~exist('carrAngle', 'var'), carrAngle = 0; end
                if isnumeric(carrier) && numel(carrier) == 1 %Carrier generation.
                    carrier = struct('A', 1, 'f', carrier, 'phi', carrAngle);
                end
                if isstruct(carrier) %Carrier generation.
                    carrier = quasiHarmonic(mySignalObj.Fs, carrier, 'time', mySignalObj.t);
                end
                carrier = setSignalVector(mySignalObj, carrier); %if isnumeric(carrier), carrier = setSignalVector(mySignalObj, carrier); end
                %Image channel checking: 'imag' - suppress imag channel; info signal is the higher heterodyne freq.
                %Image channel division based on different phase shift gotten by imag and
                %useful channels during phases subtraction with carrier.
                if strfind(mode, 'imag')
                    md1 = strrep(mode, 'imag', ''); %Modulate with necessary params.
                    %Carrier's shift will be added with different signs.
                    phi1 = 0;
                    phi2 = pi/2-phi1; %Sum of channel shifts should be pi/2.
                    [~, carrierSh] = shiftSignal(carrier, phi1);
                    [~, S1] = SSBmod(mySignalObj, carrierSh, Fs, md1);
                    [~, carrierSh] = shiftSignal(carrier, -phi2);
                    [~, S2] = SSBmod(mySignalObj, carrierSh, Fs, md1);
                    [~, S1] = shiftSignal(S1, -phi1);
                    [~, S2] = shiftSignal(S2, phi2);
                    %It's possible 2 get upper and under bands. Under is default.
                    mul = 1;
                    if strfind(mode, 'imagGetUpper'), mul = -1; end
                    SSBobj = ldivide((S1 + mul.*S2), 2); SSBsignal = SSBobj.getSignalVector; %ldivide((S1 + mul.*S2), 4)
                    return;
                end
                mult = 0; %Def - both sides.
                if strfind(mode, 'low'), mult = 1; end
                if strfind(mode, 'high')
                   mult = -1; %Subtract quaternary signal 2 get high side band.
                end
                [~, signalShifted] = shiftSignal(mySignalObj);
                modSignal = mySignalObj.*carrier; %modSignal = mySignalObj.signal.*carrier.getSignalVector;
                [~, carrierSh] = shiftSignal(carrier);
                modSignalShift = signalShifted.*carrierSh;
                SSBobj = modSignal + mult.*modSignalShift;
                if strfind(mode, 'carrier')
                    SSBobj = SSBobj + carrier; %Restore carrier if it's necessary.
                end
                SSBsignal = SSBobj.getSignalVector;
        end
        
        function [signalShifted, shiftedObj] = shiftSignal(mySignalObj, myAngle)
        %Shift all harmonics on assigned (def. pi/2) angle.
            if ~exist('myAngle', 'var'), myAngle = pi/2;  end
            if ~nnz(myAngle), shiftedObj = mySignalObj; signalShifted = double(mySignalObj); return; end
            mySignal = double(mySignalObj);
            if isreal(mySignal)
                analytSign = hilbert(mySignal);
            else
                analytSign = mySignal;
            end
            sinphaseSign = real(analytSign);
            quadrSign = imag(analytSign);
            ratio = abs(myAngle/(pi/2)); %Ratio sin/cos, or real/imag, sinphase/quadr. components.
            signalShifted = (1-ratio)*sinphaseSign - sign(myAngle)*ratio*quadrSign;
            shiftedObj = mySignalObj.setSignalVector(signalShifted);
        end
        
        function [mySpectrum, f] = getSignalSpectrum(mySignalObj, varargin)
            spectrumTypes = {'amplitude', 'phase', 'power', 'real', 'image', 'complex'};
            spec = fft(mySignalObj.signal)/length(mySignalObj.signal);
            lenMax = ceil(length(spec)/2);
            cl = strfind(varargin, 'full');
            if cellfun(@(x) ~isempty(x), cl, 'UniformOutput', true)
                lenMax = numel(spec);
            end
            spec = spec(1:lenMax);
            for ai = 1:numel(spectrumTypes)
                typeIdx = strfind(varargin, spectrumTypes{ai});
                typeIdx = cellfun(@(x) ~isempty(x), typeIdx, 'UniformOutput', true);
                typeIdx = find(typeIdx);
                if ~isempty(typeIdx), typeIdx = ai; break; end %Only ane assigned spectrum type.
            end
            if isempty(typeIdx)
                typeIdx = 1;
            end
            switch typeIdx
                case 1
                    mySpectrum = abs(spec);
                case 2
                    mySpectrum = angle(spec);
                case 3
                    mySpectrum = abs(spec).^2;
                case 4
                    mySpectrum = real(spec);
                case 5
                    mySpectrum = imag(spec);
                case 6
                    mySpectrum = spec;
                otherwise
                    mySpectrum = [];
            end
            [f, ~] = getFreqVector(mySignalObj);
            f = f(1:lenMax);
            spectrumTypes = {'log'};
            logArg = strfind(spectrumTypes, varargin);
            if ~isempty(logArg{:})
                multiplier = 10;  %Amplitude/power multiplier for log spectrum.
                if typeIdx == 2
                    multiplier = 20;
                end
                mySpectrum = multiplier*log(mySpectrum);
            end
        end
        
        function [f, df] = getFreqVector(mySignalObj)
            df = mySignalObj.Fs/length(mySignalObj.signal);
            f = 0:df:mySignalObj.Fs-df;
        end
		
		function mySignalObj = compSlowlyChangedComp(mySignalObj, varargin)
			testSign = sin(2*pi*1000*mySignalObj.t);
			noiseSignal = awgn(testSign, 1, 'measured', 'linear');
			if strcmp(varargin{1}, 'maxFreq')
				maxFreq = varargin{2};
            elseif strcmp(varargin{1}, 'linTrnd')
                stp = (varargin{3}-varargin{2})/(numel(mySignalObj.t)-1);
                mySignalObj.signal = varargin{2}:stp:varargin{3};
                return;
            end
            reset = true;
            if numel(varargin) > 2
                reset = false;
            end
            %Create an object of noise signal and get it's complex spectrum.
            NoiseObj = signalObj(mySignalObj.Fs, noiseSignal, mySignalObj.t);
            [mySpectrum, f] = getSignalSpectrum(NoiseObj, 'complex');
            %Find the closest element to assigned max freq and rest only lower samples; count filtered signal.
            [~, idx] = min(abs(f - maxFreq));
            mySpectrum(idx:end) = zeros(size(mySpectrum(idx:end)));  %Delete samples and save Fs.
            mySpectrum = [mySpectrum mySpectrum(end:-1:1)];  %Symmetrical spectrum.   end-1
            sn = real(ifft(mySpectrum));
            if reset
                mySignalObj.signal = sn/max(sn);
            else
                %Normalise slowly changed component to the signal with normalization coeff.
                sn = (sn/max(abs(sn)))*max(mySignalObj.signal)*varargin{3};
                sn = sn - mean(sn);
                if std(sn) > 10
                   disp(std(sn)); 
                end
                mySignalObj.signal = mySignalObj.signal + sn;
            end
		end
	
        function [mySignalObj, noisyComps] = setNoise(mySignalObj, SNR, myBand, filtMode)
            if ~exist('myBand', 'var'), myBand = []; end
            if ~exist('filtMode', 'var'), filtMode = ''; end
            if isempty(filtMode), filtMode = 'fft'; end
            if ~isempty(myBand) %Recompute SNR considering band limitation and filtrate noise.
                if ~iscell(myBand), myBand = {myBand}; end %Cells contain 2 numbers - band limits.
                bandWidth = sum(cellfun(@(x) diff(x), myBand));
                coeff = (mySignalObj.Fs/2)/bandWidth;
                SNR = SNR/coeff;
            end
            noisyComps = awgn(mySignalObj.signal, SNR, 'measured', 'linear');
            noisyComps = noisyComps-mySignalObj.signal;
            mySignalObjN = signalObj(mySignalObj.Fs, noisyComps); %mySignalObj.setSignalVector(noisyComps);
            if ~isempty(myBand) 
                mySignalObjN = setBand(mySignalObjN, myBand, 0, {filtMode});
                [~, noisyComps] = filtSign(mySignalObjN);
                mySignalObj.signal = mySignalObj.signal + sum(vertcat(noisyComps{:}), 1);
            else
                mySignalObj.signal = mySignalObj.signal + noisyComps; noisyComps = {noisyComps};
            end
        end
		
		function mySignalObj = sortBand(mySignalObj)
			lowFreqs = reshape(cellfun(@(x) x(1, 1), mySignalObj.band), [], 1);
			[~, idxsSorted] = sort(lowFreqs, 'ascend');
			mySignalObj.band = mySignalObj.band(idxsSorted);
			mySignalObj.bandFilter = mySignalObj.bandFilter(idxsSorted);
        end
        
        %====Windowing functions====
        
        function [window, fVect] = getFreqWindow(mySignalObj, freqRange)
            l = length(mySignalObj.signal);
            df = mySignalObj.Fs/l;
            f = 0:df:mySignalObj.Fs-df; window = 1:numel(f);
            fVect = zeros(size(mySignalObj.signal));
            if ischar(freqRange)
                if strfind(freqRange, 'Half')
                    midSamp = (l-1)/2+1; %The first smpl is DC, find middle of the rest spectrum.
                    if strfind(freqRange, 'lowHalf')
                        window = 1:floor(midSamp); fVect = zeros(size(mySignalObj.signal));
                    elseif strfind(freqRange, 'highHalf')
                        window = [1, ceil(midSamp):numel(window)]; fVect = zeros(size(mySignalObj.signal));
                    end
                end
            end
            fVect(window) = f(window);
        end
        
        %====Math functions====
        
        function mSO = plus(mySignalObj, mySignalObj2, mode)
            if exist('mode', 'var'), [mSO, s1, s2] = getOperands(mySignalObj, mySignalObj2, mode);
            else, [mSO, s1, s2] = getOperands(mySignalObj, mySignalObj2); end
            mSO = setSignalVector(mSO, s1+s2);
        end
        
        function mSO = uminus(mySignalObj, mode)
            if exist('mode', 'var')
                mSO = polyHarmonic(mySignalObj.getFs, {mySignalObj}, mySignalObj.getTimeVector);
                mSO = mSO.uminus(mode); return;
            end
            mSO = mySignalObj.setSignalVector(-mySignalObj.getSignalVector);
        end
        
        function mSO = minus(mySignalObj, mySignalObj2, mode)
            md = ''; if exist('mode', 'var'), md = ', mode'; end %Insert mode in call string.
            command = sprintf('mySignalObj2 = uminus(mySignalObj2%s);\n', md);
            command = [command sprintf('mSO = plus(mySignalObj, mySignalObj2%s);', md)];
            eval(command);
        end
        
        function mSO = ldivide(mySignalObj, mySignalObj2, mode)
            if exist('mode', 'var'), [mSO, s1, s2] = getOperands(mySignalObj, mySignalObj2, mode);
            else, [mSO, s1, s2] = getOperands(mySignalObj, mySignalObj2); end
            mSO = mSO.setSignalVector(s1./s2);
        end
        
        function mSO = times(mySignalObj, mySignalObj2, mode)
            if exist('mode', 'var'), [mSO, s1, s2] = getOperands(mySignalObj, mySignalObj2, mode);
            else, [mSO, s1, s2] = getOperands(mySignalObj, mySignalObj2); end
            mSO = setSignalVector(mSO, s1.*s2);
        end
        
        function mSO = power(mySignalObj, mySignalObj2, mode)
            if exist('mode', 'var'), [mSO, s1, s2] = getOperands(mySignalObj, mySignalObj2, mode);
            else, [mSO, s1, s2] = getOperands(mySignalObj, mySignalObj2); end
            mSO = setSignalVector(mSO, s1.^s2);
        end
        
        function myDouble = double(mySignalObj)
            myDouble = mySignalObj.getSignalVector;
        end
        
        function [mSO, s1, s2] = getOperands(mySignalObj, mySignalObj2, mode)
           %Function choosing signal object of probable SO or double variables
           %creates polyHarmonic if it's necessary, returns signal object and doubles. 
            s1 = double(mySignalObj); s2 = double(mySignalObj2);
            %Set result 2 one of objects 2 save data.
            if isa(mySignalObj, 'signalObj'), mSO = mySignalObj; end
            if isa(mySignalObj2, 'signalObj'), mSO = mySignalObj2; end
            if exist('mode', 'var')
                mSO = polyHarmonic(mSO.getFs, {mySignalObj, mySignalObj2}, mSO.getTimeVector); return;
            end
        end
        
        %Normalize signal.
        function mySignalObj = normSgnl(mySignalObj)
            mySignalObj.signal = mySignalObj.signal/max(mySignalObj.signal);
        end
        
        function l = length(mySignalObj)
            l = length(mySignalObj.signal);
        end
        
        function [anaSgn, reaCmp, imgCmp] = hilbert(mySignalObj)
            anaSgn = hilbert(double(mySignalObj));
            mySignalObj.signal = real(anaSgn); reaCmp = mySignalObj;
            mySignalObj.signal = imag(anaSgn); imgCmp = mySignalObj;
        end
        
        function sm = sum(mySignalObj)
            sm = sum(double(mySignalObj));
        end
        
        function mySignalObj = abs(mySignalObj)
            mySignalObj.signal = abs(mySignalObj.signal);
        end
        
        function pwr = signPower(mySignalObj)
            pwr = sum(double(mySignalObj).^2)/length(mySignalObj);
        end
        
        function [signObj, compt] = waveletFiltering(mySignalObj, myBand, waveletName)
            [wavCoeffs, frequencies] = cwt(mySignalObj.signal, waveletName, mySignalObj.Fs, 'VoicesPerOctave', 48);
            coefIxs = (frequencies < myBand(2))&(frequencies > myBand(1));
            cwtCoefs = wavCoeffs(coefIxs, :); %Rest only indexes according to assigned band.
            compt = icwt(cwtCoefs, waveletName);
            signObj = mySignalObj.setSignalVector(compt);
        end
        
        function [idxsMatrix, signMtx, timeMtx, signObjs] = overlpdWndws(mySignalObj, windRange, windOverlapPercent)
            if ~exist('windOverlapPercent', 'var'), windOverlapPercent = []; end
            [idxsMatrix, signMtx] = vectOverlpng(mySignalObj, 's', windRange, windOverlapPercent);
            [~, timeMtx] = vectOverlpng(mySignalObj, 't', windRange, windOverlapPercent);
            signObjs = arrayfun(@(x) signalObj(mySignalObj.Fs, signMtx(x, :), timeMtx(x, :)), 1:size(timeMtx, 1), 'UniformOutput', false);
            signObjs = cellfun(@(x) polyHarmonic(x.getFs, x, 'time', x.getTimeVector), signObjs, 'UniformOutput', false);
            signObjs = [signObjs{:}];
        end
        
        function [idxsMatrix, vectMtx] = vectOverlpng(mySignalObj, vector, windRange, windOverlapPercent)
            % Create matrix of overlapping ranges from vector assinged by windRange in percents of vector length. 
            if ischar(vector)
                switch vector
                    case 't'
                        vector = mySignalObj.t;
                    case 's'
                        vector = mySignalObj.signal;
                end
            end
            if ~exist('windOverlapPercent', 'var'), windOverlapPercent = []; end
            if isempty(windOverlapPercent), windOverlapPercent = 0; end
            Range = floor(numel(vector)*windRange/100);
            overlapSampls = floor(windOverlapPercent*Range/100);
            framesNumber = floor(numel(vector)/(Range-overlapSampls));
            idxsMatrix = zeros(framesNumber, Range); vectMtx = idxsMatrix;
            increment = Range - overlapSampls;
            startPosition = 1;

            idVect = 1:numel(vector);
            for i=1:1:framesNumber
                if startPosition+Range > idVect(end), break; end
               idxsMatrix(i,:) = idVect(startPosition:startPosition+Range-1);
               vectMtx(i,:) = vector(idxsMatrix(i,:));
               startPosition = startPosition + increment;
            end
        end
        
        %==Window functions==
        
        function [frms, sampleFrames] = timeFrames(mySignalObj, lngth, overlap, mode)
            if ~exist('mode', 'var'), mode = ''; end
            %Translate from percents to seconds.
            if nnz(strfind(mode, 'lenPercent'))
                lngth = (mySignalObj.t(end)-mySignalObj.t(1))*lngth/100; end
            if nnz(strfind(mode, 'lapPercent'))
                overlap = lngth*overlap/100; end
            %Translate from seconds to samples.
            dt = mySignalObj.t(2)-mySignalObj.t(1);
            lngth = ceil(lngth/dt);
            overlap = ceil(overlap/dt);
            shft = lngth-overlap;
            frms = zeros(ceil(length(mySignalObj)/shft)-1, lngth);
            sampleFrames = zeros(ceil(length(mySignalObj)/shft)-1, 2);
            for ai = 1:size(frms, 1)
                frms(ai, :) = mySignalObj.signal(shft*(ai-1)+1:shft*(ai-1)+lngth); 
                sampleFrames(ai, :) = [shft*(ai-1)+1 shft*(ai-1)+lngth];
            end
        end
        
        function [frms, sampleFrames] = winAverage(mySignalObj, lngth, overlap, mode)
            if ~exist('mode', 'var'), mode = ''; end
            %Process by quasiconstant period.
            phD = compFullPhase(mySignalObj, 'deriv');
            if iscell(phD), phD = phD{1}; end
            mySignal = mySignalObj.signal;
            sampling = find(phD<0);
            cutO = sampling(1):sampling(end);
            if nnz(strfind(mode, 'cutOff'))
                mySignal = mySignal(cutO);
                phD = phD(cutO);
                sampling = find(phD<0);
            end
            if nnz(strfind(mode, 'maxPeriod'))
                %Distances between derivative peaks; the first sample and the first peak.
                sampling(1) = 0;
                dt = diff(sampling);
                dPrev = 0; %The last previous window sample.
                frms = NaN(numel(dt), max(dt));
                sampleFrames = NaN(numel(dt), 2);
                for ai = 1:numel(dt)
                    sampleFrames(ai, 1) = sampling(ai)+1;
                    sampleFrames(ai, 2) = sampling(ai+1);
                    dPrev = sampleFrames(ai, 2);
                    frms(ai, 1:dPrev-sampleFrames(ai, 1)+1) = mySignal(sampleFrames(ai, 1):dPrev);
                end
            end
        end
        
        function [frms] = freqFrames(mySignalObj, lngth, overlap, mode)
            if ~exist('mode', 'var'), mode = ''; end
            [f, df] = getFreqVector(mySignalObj);
            %Translate from percents to Hz.
            if nnz(strfind(mode, 'lenPercent'))
                lngth = mySignalObj.Fs*lngth/100; end
            if nnz(strfind(mode, 'lapPercent'))
                overlap = lngth*overlap/100; end
            %Translate from seconds to samples.
            lngth = ceil(lngth/df);
            overlap = ceil(overlap/df);
            shft = lngth-overlap;
            frms = zeros(ceil(length(mySignalObj)/shft)-1, lngth);
            for ai = 1:size(frms, 1)
                frms(ai, :) = f(shft*(ai-1)+1:shft*(ai-1)+lngth); end
        end
        
        function [frqWind, sampRng, wind] = getFreqWinds(mySignalObj, range, wind)
            if ~exist('wind', 'var'), wind = ''; end
            [mySpectrum, f] = getSignalSpectrum(mySignalObj);
            [~, sampRng] = min(abs(range(1)-f));
            [~, sampRng(2)] = min(abs(range(end)-f));
            frqWind = mySpectrum(sampRng(1):sampRng(2));
            if ~isempty(wind)
                [wind, frqWind] = windowing(mySignalObj, wind, [], frqWind); end
        end
        
        function [wind, signal, sObjOut] = windowing(mySignalObj, wind, winLen, signal, varargin)
            if ~exist('signal', 'var'), signal = []; end
            if isempty(signal), signal = mySignalObj.signal; end
            if ~exist('winLen', 'var'), winLen = []; end
            if isempty(winLen), winLen = length(signal); end
            if ~exist('wind', 'var'), wind = ''; end
            if isempty(wind), wind = 'rectwin'; end
            if ~isnumeric(wind)
                wind = window(wind, winLen, varargin{:}); end
            wind = wind(end-max(size(signal))+1:end); %GAG: use right window half default if window is longer.
            wind = reshape(wind, size(signal));
            if nargout > 1
                signal = signal.*wind; end
            if nargout > 2
                sObjOut = mySignalObj.setSignalVector(signal); end
        end
        
        function [SO, signal] = corr(mSO1, mSO2)
            if nargin<2, mSO2 = mSO1; end
            signal = xcorr(double(mSO1), double(mSO2), 'biased');
            start = (numel(signal)-1)/2;
            signal = signal(start+1:end);
            SO = signalObj(mSO1.Fs, signal, mSO1.t);
        end
        
        function [spec, f, crossCorr, crosSpec, spec1, spec2, epsn] = coherenceFun(mSO1, mSO2)
            crossCorr = corr(mSO1, mSO2);
            [spec, f] = crossCorr.getSignalSpectrum('complex');
            crosSpec = abs(spec).^2;
            spec1 = getSignalSpectrum(mSO1.corr, 'complex');
            spec2 = getSignalSpectrum(mSO2.corr, 'complex');
            spec1 = abs(spec1); spec2 = abs(spec2);
            [~, noiseLevelVector, ~, ~] = logSpectrum(struct('signal', double(mSO1), 'spectrum', spec1', 'frequencies', f', 'Fs', mSO1.Fs), struct('rmsFactor', '2'), 'correl');
            idxs1=spec1<noiseLevelVector'; spec1(idxs1) = zeros(size( spec1(idxs1) ));
            [~, noiseLevelVector, ~, ~] = logSpectrum(struct('signal', double(mSO2), 'spectrum', spec2', 'frequencies', f', 'Fs', mSO2.Fs), struct('rmsFactor', '2'), 'correl');
            idxs2=spec2<noiseLevelVector'; spec2(idxs2) = zeros(size( spec2(idxs2) ));
            [~, noiseLevelVector, ~, ~] = logSpectrum(struct('signal', double(mSO1), 'spectrum', crosSpec', 'frequencies', f', 'Fs', mSO1.Fs), struct('rmsFactor', '2'), 'correl');
            idxs3=crosSpec<noiseLevelVector'; crosSpec(idxs3) = zeros(size( crosSpec(idxs3) ));
            idxs = idxs1 | idxs2 | idxs3;
            denomina = spec1.*spec2; epsn = 0;
            spec = crosSpec./(denomina+eps); %+rms(denomina)*2
            spec(idxs) = zeros(size( spec(idxs) ));
        end
        
        function signal = getStep(mSO, situation, ampl, mode)
            [~, sampleFrames] = timeFrames(mSO, situation, 0, mode);
            signal = zeros(size(mSO.t));
            signal(sampleFrames(2, 1):end) = ones(size( signal(sampleFrames(2, 1):end) ))*ampl;
        end

    end
    
end