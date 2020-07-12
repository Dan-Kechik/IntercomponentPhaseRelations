classdef polyHarmonic < quasiHarmonic

	properties(SetAccess = protected, GetAccess = public)
		%Set of discrete quasiharmonical components.
		components
        componentsOrig
		%Some components may be modulated by a several other - ampl., angle mod.
		modulations

		nonLinearOperator
		%Pairs paramName-value: mod, filt, non-lin. Specify operations and their order.
		operations
	end


	methods(Access = public)
        
		function myHarmonic = polyHarmonic(myFs, mySignalComponents, varargin)
            %Set Fs, compute a time vector, add a signal components,
            %resample all signans, compute a signal.
            myHarmonic = myHarmonic@quasiHarmonic(myFs, [], varargin{:}); %Don't set a parameters, but set a time if it's parameters assigned.
            if ~exist('mySignalComponents', 'var')
                mySignalComponents = [];  %Create a signal with standart filling, that will be created in quasiHarmonic constructor.
            end
            myHarmonic = addSignalComponents(myHarmonic, mySignalComponents);
        end
        
		function myHarmonic = addSignalComponents(myHarmonic, mySignalComponents, reset)
            %Add a signal components, resample added signans, compute a signal.
            %Components may be a structs, or a signal objects, or vectors, lying in
            %arrays or arrays of cells contain a different data types.
            %Original data lying in original field, internal format is object.
			if ~exist('reset', 'var')
				reset = false;
			end
			if reset
				myHarmonic.signal = [];
                myHarmonic.componentsOrig = [];
                myHarmonic.components = [];
            else
                myHarmonic.componentsOrig = myHarmonic; nO = numel(myHarmonic.components); nC = numel(mySignalComponents);
                myHarmonic.operations = [myHarmonic.operations sprintf('New components %d-%d were added.', nO+1, nC)];
			end
			for i = 1:numel(mySignalComponents)
				if iscell(mySignalComponents(i)) %Extract from cell.
					myCurrSignal = mySignalComponents{i};
                else
                    myCurrSignal = mySignalComponents(i);
                end
                %Add a original data.
				%myHarmonic.componentsOrig = [myHarmonic.componentsOrig, {myCurrSignal}];
				if isnumeric(myCurrSignal) %If it's a signal vector.
					myCurrSignal = signalObj(myHarmonic.Fs, myCurrSignal, myHarmonic.t);
				end
				if isstruct(myCurrSignal) %If it's struct with component parameters.
					myCurrSignal = quasiHarmonic(myHarmonic.Fs, myCurrSignal, 'time', myHarmonic.t);
                end
                if isa(myCurrSignal, 'polyHarmonic')  %sum( strcmp({'signalObj', 'quasiHarmonic'}, class(myCurrSignal)) )
                    myHarmonic.t = myCurrSignal.t;
                    myHarmonic.Fs = myCurrSignal.Fs;
                    myCurrSignal = reshape(myCurrSignal.components, 1, []);
                end
                %Add a object of signal component.
                if ~strcmp(class(myCurrSignal), 'signalObj')
                    myCurrSignal = signalObj(myHarmonic.Fs, double(myCurrSignal), myHarmonic.t);
                end
				myHarmonic.components = [myHarmonic.components, myCurrSignal];
            end
			myHarmonic = sortComponents(myHarmonic);
            myHarmonic = compSignalVector(myHarmonic);
        end
        
        function myHarmonic = setBaseFreqLabel(myHarmonic, myLabel)
            for ai = 1:numel(myHarmonic.components)
                myHarmonic.components(ai) = myHarmonic.components(ai).setLabel(sprintf('%s - %d', myLabel, ai));
            end
        end
        
		function myComponents = getSignalComponents(myHarmonic, mode)
            if ~exist('mode', 'var')
                mode = '';
            end
            if strcmp(mode, 'orig')
                myComponents = myHarmonic.componentsOrig;
            else
                myComponents = myHarmonic.components;
            end
        end
        
        function myHarmonic = compTimeVector(myHarmonic, myTime, mode)
            %Set time vector from Fs and pointed length (in seconds (def), in periods, in samples), if time is a scalar, or pointed time vector.
            if ~exist('mode', 'var')
                mode = [];
            end
            if strcmp(mode,  'periods')
    %             %Translate to seconds.
    %             myHarmonic = compTimeVector(myHarmonic, myHarmonic.Fs, 'samples'); %For period measurement.
    %             F0 =  getCentralFreq(myHarmonic);
    %             T = 1/F0;  %Seconds in the signal period.
    %             myTime = myTime*T; %Time in pointed periods number.
                mode = 'seconds';
            end
            myHarmonic = compTimeVector@signalObj(myHarmonic, myTime, mode);
            myHarmonic = resampleToTime(myHarmonic);
        end
        
        function myHarmonic = resampleToTime(myHarmonic)
            myComponents = myHarmonic.components;
            for i = 1:numel(myComponents)
                myComponents(i) = compTimeVector(myComponents(i), myHarmonic.t);
                myComponents(i) = resampleToTime(myComponents(i));
            end
            myHarmonic = compSignalVector(myHarmonic);
        end

		function myHarmonic = compSignalVector(myHarmonic)
			mySignal = zeros(size(myHarmonic.t));
			for i = 1:numel(myHarmonic.components)
                if numel(myHarmonic.components(i).getTimeVector) ~= numel(mySignal)
                    myHarmonic.components(i) = compTimeVector(myHarmonic.components(i), myHarmonic.t);
                end
				myHarmonic.components(i) = myHarmonic.components(i).compSignalVector;
				%mySignal = mySignal + myHarmonic.components(i).getSignalVector;
                mySignal = bsxfun(@plus, mySignal, myHarmonic.components(i).getSignalVector);
            end
            myHarmonic.signal = mySignal;
        end
        
        function myHarmonic = setSignalVector(myHarmonic, mySignal)
            if ~iscell(mySignal), mySignal = {mySignal}; end
            myHarmonic = addSignalComponents(myHarmonic, mySignal, 1);
        end
		
        function myHarmonic = processSignal(myHarmonic)
			for i = 1:size(myHarmonic.operations, 2)
				myHarmonic = clearOperations(myHarmonic);
				currOperator = myHarmonic.operations{1, i};
				myHarmonic.(currOperator) = myHarmonic.operations{2, i};
				myHarmonic = compNonLinearOperator(myHarmonic);
				myHarmonic = compModulations(myHarmonic);
				myHarmonic = filtSign(myHarmonic);
			end
		end

        function myHarmonic = filtSign(myHarmonic, myFilt)
			if ~exist('myFilt', 'var')
				myFilt = [];
			end
			%Filter signal and get a separate components accord to filter set.
			[myHarmonic, myComponents] = filtSign@signalObj(myHarmonic, myFilt);
		    for i = 1:numel(myComponents)
				myComponentSignals(i) = signalObj(myHarmonic.Fs, myComponents{i}, myHarmonic.t);
            end
            if ~isempty(myComponents)
                myHarmonic.componentsOrig = myHarmonic;
                myHarmonic.operations = [myHarmonic.operations sprintf('Signal was filted.')];
                myHarmonic.components = myComponentSignals;
            end
        end
        
        function [myHarmonic, noisyComps] = setNoise(myHarmonic, SNR, myBand, filtMode)
            if ~exist('myBand', 'var'), myBand = []; end
            if ~exist('filtMode', 'var'), filtMode = ''; end
            if isempty(filtMode), filtMode = 'fft'; end
            %-Check if a separate component SNR was assigned-
            cmpNm = 0; if numel(SNR) == 2, cmpNm = SNR(2); SNR = SNR(1); end
            if nnz(cmpNm) %Use as reference assigned component.
                [~, noisyComps] = setNoise(myHarmonic.getComponent(cmpNm), SNR, myBand, filtMode);
            else
                [~, noisyComps] = setNoise@signalObj(myHarmonic, SNR, myBand, filtMode);
            end
            myHarmonic = addSignalComponents(myHarmonic, noisyComps, 0);
            myHarmonic.operations = [myHarmonic.operations sprintf('Noise was added.')];
        end
		
		function myHarmonic = compNonLinearOperator(myHarmonic)
			mySignal = myHarmonic.signal;
			myOperator = myHarmonic.nonLinearOperator;
			if ~isempty(myOperator)
				myOperator = strrep(myOperator, 'x', 'mySignal.');
				mySignal = eval(myOperator);
                myHarmonic.signal = mySignal;
                mySignal = signalObj(myHarmonic.Fs, mySignal, myHarmonic.t);
                myHarmonic.components = mySignal;
                myHarmonic.componentsOrig = myHarmonic;
                myHarmonic.operations = [myHarmonic.operations sprintf('Non-linear operator %s was applied.', myOperator)];
            end
		end
		
		function myHarmonic = compModulations(myHarmonic)
		%Assign a carrier component and modulation signals by their indexes.
		%Also original (carr and mod) harmonics that need 2 be rest can be assigned.
		%Carrier and side freqs are the one component, that mb filtered and divided
		%on components, and it's mb only necessary include in the common signal.
			%myHarmonic.componentsOrig = myHarmonic.components;
			for i = 1:numel(myHarmonic.modulations)
				%Find necessary components and make a mod signal.
				carr = myHarmonic.modulations(i).carrierNum{1};
				%myHarmonic.includeVect(myHarmonic.modulations(i).carrierNum) = myHarmonic.modulations(i).restCarr;
				carr = getComponent(myHarmonic, carr);
				for j = 1:numel(myHarmonic.modulations(i).modNum)
					modComp = myHarmonic.modulations(i).modNum{j};
					modComp = getComponent(myHarmonic, modComp);
					%Variate a signal param - i.e. modulate.
					%modType - type of modulation (ampl/freq/phase), that points which parameter variate.
					param = carr.getParameters(myHarmonic.modulations(i).modType) + modComp.signal;
					carr = carr.setParameters(myHarmonic.modulations(i).modType, param);
					%myHarmonic.includeVect(myHarmonic.modulations(i).modNum(j)) = myHarmonic.modulations(i).restMod;
				end
				restComp = myHarmonic.modulations(i).restComp;  %Rest or not original components vector.
				compNumbers2del = ~myHarmonic.modulations(i).modNum{restComp};  %Numbers of modulations to rest.
                myHarmonic.componentsOrig = myHarmonic; nM = myHarmonic.modulations(i).carrierNum{1};
                myHarmonic.operations = [myHarmonic.operations sprintf('Component %d was added.', nM)];
                if ~myHarmonic.modulations(i).restCarr
                    compNumbers2del = [compNumbers2del myHarmonic.modulations(i).carrierNum{1}];  %Plus carrier.
                end
                compNumbers2del = find(compNumbers2del);
				%Delete unnecessary components.
				myHarmonic = deleteComponents(myHarmonic, compNumbers2del);
				%Add modulated signal.
				myHarmonic.components = [myHarmonic.components, carr];
			end
			%Sort components by F0.
            myHarmonic = sortComponents(myHarmonic);
            myHarmonic = compSignalVector(myHarmonic);
		end
		
		function myHarmonic = deleteComponents(myHarmonic, numbers2del)
			theCurrComponents = myHarmonic.components;
			compNums = 1:numel(theCurrComponents);
			compNums = (compNums == numbers2del); %Logical vector with found elems 2 del.
			compNums = find(~compNums);  %Numbers of elements 2 rest.
			myHarmonic.components = theCurrComponents(compNums);
		end
		
		function [myHarmonic, idxs] = sortComponents(myHarmonic)
% 			fieldNames = {'components', 'componentsOrig'};
% 			for i = 1:numel(fieldNames)
% 				theCurrComponents = myHarmonic.(fieldNames{i});
% 				centralFreqs = arrayfun(@(x) x.getCentralFreq, theCurrComponents);
% 				[~, idxs] = sort(centralFreqs);
% 				myHarmonic.(fieldNames{i}) = theCurrComponents(idxs);
%             end
            centralFreqs = arrayfun(@(x) x.getCentralFreq, myHarmonic.components);
            [~, idxs] = sort(centralFreqs);
            myHarmonic.components = myHarmonic.components(idxs);
		end
		
		function myComponent = getComponent(myHarmonic, numb)
			%numb is a string that assigns a components position in original components array
			%or in the current comp. arr., ect.
            if ~exist('numb', 'var'), numb = []; end
            if isempty(numb), numb = 1:numel(myHarmonic.components); end
			if isnumeric(numb)
				myComponent = myHarmonic.components(numb);
				return;
			end
			orig = strfind(numb, 'orig');
			orig = logical(find(orig));  %A flag of original vector.
			if orig
				myComps = myHarmonic.componentsOrig;
				numb = strrep(numb, 'orig', '');
			else
				myComps = myHarmonic.components;
			end
			numb = str2double(numb);
			myComponent = myComps(numb);
		end
		
		function myHarmonic = addModulations(myHarmonic, modStruct)
			%modStruct: harmonic numbers - carrier and mod components, type (ampl, freq, phase), index mod, delete from signal carrier or/and mod function.
			myHarmonic.modulations = [myHarmonic.modulations, modStruct];
        end
        
		function myHarmonic = setNonLinearOperator(myHarmonic, myOperator)
            myHarmonic.nonLinearOperator = myOperator;
        end

        function myHarmonic = clearOperations(myHarmonic)
			myHarmonic.modulations = [];
			myHarmonic.band = [];
			myHarmonic.nonLinearOperator = [];
        end
        
        function viewComponents(myHarmonic, indexes)
            myComponents = myHarmonic.components;
            if ~exist('indexes', 'var')
                indexes = 1:numel(myComponents);
            end
            if ~numel(myComponents)
                fprintf('Components of your signal are not setted.\n');
                return;
            end
            myComponents = myComponents(indexes);
            for i = 1:numel(myComponents)
                plotSignal(myComponents(i),  'amplitude');
            end
        end
        
        function viewFilterSet(myHarmonic, indexes)
            myFilterSet = myHarmonic.bandFilter;
            if ~exist('indexes', 'var')
                indexes = 1:numel(myFilterSet);
            end
            if ~numel(myFilterSet)
                fprintf('Filters on your signal are not setted.\n');
                return;
            end
            myFilterSet = myFilterSet(indexes);
            for i = 1:numel(myFilterSet)
                fvtool(myFilterSet{i});
            end
        end
    
        function [signalShifted, myHarmonic] = shiftSignal(myHarmonic, myAngle, numbers)
            if ~exist('myAngle', 'var'), myAngle = pi/2;  end
            if ~exist('numbers', 'var'), numbers = ''; end
            if isempty(numbers), numbers = 1:numel(myHarmonic.components); end
            myHarmonic.componentsOrig = myHarmonic;
            myHarmonic.operations = [myHarmonic.operations 'Harmonics shifting.'];
            for i = numbers
                [~, myHarmonic.components(i)] = shiftSignal(myHarmonic.components(i), myAngle);
            end
            myHarmonic = myHarmonic.compSignalVector; signalShifted = myHarmonic.getSignalVector;
        end
        
        function myFullPhase = compFullPhase(myHarmonic, mode)
			if ~exist('mode', 'var'), mode = 0; end
            if isstruct(mode), myFullPhase = compFullPhase@quasiHarmonic(myHarmonic, mode); return; end %4 quasiHarm-c functions.
            if ~isnumeric(mode)
                if strfind(mode, 'whole'), myFullPhase = compFullPhase@signalObj(myHarmonic, mode); return;
                    end %Compute phase of the whole signal.
                if nnz(strfind(mode, 'plot'))&&~nnz(strfind(mode, 'gcf')), figure; mode = [mode 'gcf'];
                    end %If plotting is assigned, plot at the one figure.
            end
            myFullPhase = arrayfun(@(x) compFullPhase(x, mode), myHarmonic.components, 'UniformOutput', false);
            if strfind(mode, 'numDiv'), myFullPhase = cellfun(@(x, y) x/y, myFullPhase, num2cell(1:numel(myFullPhase)), 'UniformOutput', false); end
        end
        
        function writeAudio(myHarmonic, fileName, indexes, mode)
            if ~exist('indexes', 'var'), indexes = []; end
            if ~exist('mode', 'var'), mode = ''; end
            if isempty(indexes), indexes = 1:numel(myHarmonic.components); end
            comps = arrayfun(@(x) x.getSignalVector, myHarmonic.components(indexes), 'UniformOutput', false);
            comps = cellfun(@(x) reshape(x, [], 1), comps, 'UniformOutput', false);
            if strfind(mode, 'sep')
                
            end
            y = horzcat(comps{:});
            if strfind(mode, 'sum'), y = sum(y, 2); end
            y = y/sum(max(abs(y))); %Normalize 2 avoid clipping.
            audiowrite(fileName, y, myHarmonic.Fs);
        end
        
        function myHarmonic = readAudio(myHarmonic, fileName, reset)
            [y, Fs] = audioread(fileName); mSO = signalObj(Fs, sum(y, 2)');
            myHarmonic = addSignalComponents(myHarmonic, mSO, reset);
        end
        
        function [signPolyH, params, S, exitflag] = approx(myHarmonic, alph, tolvalue, meth)
            if ~exist('alph', 'var'), alph = []; end
            if isempty(alph), alph = 1; end
            if ~exist('tolvalue', 'var'), tolvalue = []; end
            if isempty(tolvalue), tolvalue = 0.08; end
            if ~exist('meth', 'var'), meth = ''; end
            if isempty(meth), meth = 'approxFunMits'; end
            %=Mitsianok formula=
            %Get a and b.
            w = arrayfun(@(x) 2*pi*x.getCentralFreq, myHarmonic.components, 'UniformOutput', false);
            t = myHarmonic.getTimeVector; phi0 = cellfun(@(x) x*t, w, 'UniformOutput', false); %Base tone.
            phi = arrayfun(@(x) compFullPhase(x, 'unwrap') + pi/2, myHarmonic.components, 'UniformOutput', false);
            A = arrayfun(@(x) x.getAmplitude, myHarmonic.components, 'UniformOutput', false);
            switch meth
                case 'approxFunMits'
                    phi1 = cellfun(@(x, y) x - y, phi, phi0, 'UniformOutput', false);
                    A0 = cellfun(@(x) mean(x), A, 'UniformOutput', false);
                    A1 = cellfun(@(x, y) x - y, A, A0, 'UniformOutput', false);
                    a = cellfun(@(x, y, z) x.*cos(z) + y.*cos(z), A0, A1, phi1, 'UniformOutput', false);
                    b = cellfun(@(x, y, z) x.*sin(z) + y.*sin(z), A0, A1, phi1, 'UniformOutput', false);
%                     cosA = cellfun(@(x, y) x.*cos(y), b, phi0, 'UniformOutput', false);
%                     sinA = cellfun(@(x, y) x.*sin(y), a, phi0, 'UniformOutput', false);
                    varArr = horzcat(a{:}, b{:});
                case 'approxFunIntc'
                    varArr = horzcat(A{:}, phi{:});
            end
            options=optimset('tolX', tolvalue, 'MaxFunEvals', 3000*numel(varArr)); %, 'Display', 'iter', 'PlotFcns', @optimplotfval
            fun = @(x)approxFun(myHarmonic, x, alph, meth);
            [params, S, exitflag] = fminsearch(fun, varArr, options);
            %Synthesize approximation signal.
            n = length(myHarmonic);
            %b0 = varArr(1:n); 
            aVe = cell(size(myHarmonic.components)); bVe = aVe;
            for ai = 0:numel(myHarmonic.components)-1
                aVe{ai+1} = params(ai*n+1:(ai+1)*n); end
            for bi = ai+1:ai+numel(myHarmonic.components)
                bVe{bi-ai} = params(bi*n+1:(bi+1)*n); end
            w = arrayfun(@(x) 2*pi*x.getCentralFreq, myHarmonic.components, 'UniformOutput', false);
            sinComp = cellfun(@(x, y) x.*sin(y*t+pi/4), aVe, w, 'UniformOutput', false);
            cosComp = cellfun(@(x, y) x.*cos(y*t+pi/4), bVe, w, 'UniformOutput', false);
            cmps = cellfun(@(x, y) x+y, sinComp, cosComp, 'UniformOutput', false);
            signPolyH = myHarmonic.addSignalComponents(cmps, 1);
        end
        
        function S = approxFun(myHarmonic, varArr, alph, meth)
            %Function translates multivariable vector varArr into two cell
            %parameter arrays.
            if ~exist('meth', 'var'), meth = ''; end
            if isempty(meth), meth = 'approxFunMits'; end
            n = length(myHarmonic);
            %b0 = varArr(1:n); 
            aVe = cell(size(myHarmonic.components)); bVe = aVe;
            for ai = 0:numel(myHarmonic.components)-1
                aVe{ai+1} = varArr(ai*n+1:(ai+1)*n); end
            for bi = ai+1:ai+numel(myHarmonic.components)
                bVe{bi-ai} = varArr(bi*n+1:(bi+1)*n); end
            S = myHarmonic.(meth)(aVe, bVe, alph); %b0, 
        end
        
        function S = approxFunMits(myHarmonic, aVe, bVe, alph) %, b0
            %Function implements Mitsianok's objective function:
            %S = sum([y(ti)-yEst(ti)]^2) + alph*sum()...
            %b0 - count beginning; aVe - sin amplitude vector; bVe = cos amplitude vector.
            %aVe and bVe are cell arrays, containing vectors for each harmonic.
            bet = alph(end); if numel(alph) == 2, alph = alph(1); end %alpha = beta default - sin/cos and intercomponent coeffs resp.
            y = double(myHarmonic);
            %yEst - your signal estimation computed from gotten parameters and average freqz.
            t = myHarmonic.getTimeVector;
            w = arrayfun(@(x) 2*pi*x.getCentralFreq, myHarmonic.components, 'UniformOutput', false);
            sinComp = cellfun(@(x, y) x.*sin(y*t), aVe, w, 'UniformOutput', false);
            cosComp = cellfun(@(x, y) x.*cos(y*t), bVe, w, 'UniformOutput', false);
            yEst = sum(vertcat(sinComp{:}), 1) + sum(vertcat(cosComp{:}), 1); %b0 + 
            signMSE = sum((y-yEst).^2); %/length(y); %Error measure should not depend on signal length.
            %b0MSD = [0 alph*sum(diff(b0).^2)];
            aVeMSD = cellfun(@(x) alph*sum(diff(x).^2), aVe, 'UniformOutput', false);
            bVeMSD = cellfun(@(x) alph*sum(diff(x).^2), bVe, 'UniformOutput', false);
            %=Add intercomponent measure=
            phiVe = cellfun(@(x, y) atan(x./y), bVe, aVe, 'UniformOutput', false);
            intcptMSD = bet*sum(diff(phiVe{1}-phiVe{2}/2).^2);
            S = signMSE + sum(vertcat(aVeMSD{:}), 1) + sum(vertcat(bVeMSD{:}), 1) + sum(intcptMSD); % + b0MSD
        end
        
        function S = approxFunIntc(myHarmonic, aVe, phiVe, alph) %, b0
            %Function implements objective function considerig intercomponent relations:
            %S = sum([y(ti)-yEst(ti)]^2) + alph*sum()...
            %b0 - count beginning; aVe - amplitude vector; phiVe = phase vector.
            %aVe and bVe are cell arrays, containing vectors for each harmonic.
            y = double(myHarmonic);
            %yEst - your signal estimation computed from gotten parameters and average freqz.
            comps = cellfun(@(x, y) x.*sin(y), aVe, phiVe, 'UniformOutput', false);
            yEst = sum(vertcat(comps{:}), 1); %b0 + 
            signMSE = sum((y-yEst).^2)/length(y); %Error measure should not depend on signal length.
            %b0MSD = [0 alph*sum(diff(b0).^2)];
            aVeMSD = cellfun(@(x) alph*sum(diff(x).^2), aVe, 'UniformOutput', false);
            bVeMSD = cellfun(@(x) alph*sum(diff(x).^2), phiVe, 'UniformOutput', false); %Too great1!
            intcptMSD = alph*sum(diff(phiVe{1}-phiVe{2}/2).^2);
            S = signMSE + sum(vertcat(aVeMSD{:}), 1) + sum(vertcat(bVeMSD{:}), 1) + sum(intcptMSD); % + b0MSD
        end
        
        function [SNRest, pwrCln, pwrNoi] = SNRestClean(cleanSgn, noisySgn)
            %Function estimates each component and the whole SNR (the last element) using
            %formula: SNR = clean_pow/(noisy_pow - clean_pow);
            pwrCln = arrayfun(@(x) signPower(x), cleanSgn.getSignalComponents);
            pwrNoi = arrayfun(@(x) signPower(x), noisySgn.getSignalComponents);
            %pwrCln = [pwrCln signPower(cleanSgn)];
            %pwrNoi = [pwrNoi signPower(noisySgn)];
            SNRest = pwrCln./abs(pwrNoi - pwrCln);
        end
        
        function [pwr, compsPwr] = signPower(myHarmonic)
            compsPwr = arrayfun(@(x) signPower(x), myHarmonic.getSignalComponents);
            pwr = sum(compsPwr);
        end
        
        %Function implements non-linear harmonics production of assigned
        %base tone and harmonics number saving instantaneous amplitude and phase.
        function myHarmonic = modelHarmonicsNonLin(myHarmonic, harmoNumbs, harmoCoeffs, baseToneIdx)
            if ~exist('harmoNumbs', 'var'), harmoNumbs = []; end
            if isempty(harmoNumbs), harmoNumbs = 2; end
            if ~exist('harmoCoeffs', 'var'), harmoCoeffs = []; end
            if isempty(harmoCoeffs), harmoCoeffs = ones(size(harmoNumbs)); end
            if ~exist('baseToneIdx', 'var'), baseToneIdx = []; end
            if isempty(baseToneIdx), baseToneIdx = 1; end
            firstHarmc = polyHarmonic(myHarmonic.Fs, myHarmonic.components(baseToneIdx), 'time', myHarmonic.t);
            myOperator = arrayfun(@(x, y) sprintf('%f*x^%f', x, y), harmoCoeffs, harmoNumbs, 'UniformOutput', false);
            myOperator = strjoin(myOperator, '+');
            firstHarmc = setNonLinearOperator(firstHarmc, myOperator);
            firstHarmc = compNonLinearOperator(firstHarmc);
            myHarmonic.components(baseToneIdx) = myHarmonic.components(baseToneIdx).setSignalVector(double(firstHarmc));
            myHarmonic = compSignalVector(myHarmonic);
        end
        
        %Function implements non-linear harmonics production of assigned
        %base tone and harmonics number saving instantaneous amplitude and phase.
        function [myHarmonic, componentsAdded] = modelHarmonicsHlb(myHarmonic, harmoNumbs, harmoCoeffs, harmoAngles, baseToneIdx)
            if ~exist('harmoNumbs', 'var'), harmoNumbs = []; end
            if isempty(harmoNumbs), harmoNumbs = 2; end
            if ~exist('harmoCoeffs', 'var'), harmoCoeffs = []; end
            if isempty(harmoCoeffs), harmoCoeffs = ones(size(harmoNumbs)); end
            if ~exist('harmoAngles', 'var'), harmoAngles = []; end
            if isempty(harmoAngles), harmoAngles = zeros(size(harmoNumbs)); end
            if ~exist('baseToneIdx', 'var'), baseToneIdx = []; end
            if isempty(baseToneIdx), baseToneIdx = 1; end
            
            myAmplit = getAmplitude(myHarmonic.components(baseToneIdx));
            myAmplit = arrayfun(@(x) myAmplit*x, harmoCoeffs, 'UniformOutput', false);
            myPhs = compFullPhase(myHarmonic.components(baseToneIdx), 'unwrap');
            myPhs = arrayfun(@(x, y) myPhs*x+y, harmoNumbs, harmoAngles, 'UniformOutput', false); %+pi/2
            myCompts = cellfun(@(x, y) struct('A', x, 'f', 0, 'phi', y), myAmplit, myPhs, 'UniformOutput', false);
            myCompts = polyHarmonic(myHarmonic.Fs, myCompts, 'time', myHarmonic.t);
            if nargout == 2, componentsAdded = myCompts.getComponent(1:numel(harmoNumbs)); end
            %Measure harmonics amplitudes, compare with required, both normalized.
            [~, compsPwr] = signPower(myCompts);
            compsPwr = compsPwr/compsPwr(1); harmoCoeffs = harmoCoeffs/harmoCoeffs(1);
            correctionCoeffs = harmoCoeffs./compsPwr;
            
            myHarmonic.components(baseToneIdx) = myHarmonic.components(baseToneIdx).setSignalVector(double(myCompts));
            myHarmonic = compSignalVector(myHarmonic);
        end
        
    end
    


end