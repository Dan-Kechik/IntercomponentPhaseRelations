classdef quasiHarmonic < signalObj

properties(SetAccess = private, GetAccess = public)
%Original amplitude, frequency, phase struct.
signalParams

A
f
phi

noiseLevel
end

properties(SetAccess = protected, GetAccess = public)

    %Real signal parameters.
    amplitude
    fullPhase
    
end

methods(Access = public)
%Constructor.
	function myHarmonic = quasiHarmonic(myFs, mySignalParams, varargin)
        if ~exist('mySignalParams', 'var')
            mySignalParams = [];
        end
        if isempty(mySignalParams)
            mySignalParams = struct('A', 0, 'f', 1, 'phi', 0);
        end
        myHarmonic = myHarmonic@signalObj(myFs);
        %Struct with signal parameters - amplitude, frequency, phase, that
        %may be a constants, vectors (functions of time), a signals.
        %Struct also can contain variative params that sh.b. translated to signals.
        fNms = fieldnames(mySignalParams); slowlyArgs = cell(size(fNms)); iterations = zeros(size(fNms));
        for i = 1:numel(fNms)
            if numel(mySignalParams.(fNms{i})) == 1, continue; end
            if isa(mySignalParams.(fNms{i}), 'double'), continue; end
            slowlyArgs{i} = mySignalParams.(fNms{i})(2:end); %Slowly args are in multiple cells.
            mySignalParams.(fNms{i}) = mySignalParams.(fNms{i}){1}; iterations(i)=1; %Set average value and processing this field after.
        end
		myHarmonic.signalParams = mySignalParams;
        %varargin may also assign a signal parameters and set a time vector.
        if numel(varargin)
            %Set assigned params, including time.
            myHarmonic = setParameters(myHarmonic, varargin{:});
        end
        if ~sum(myHarmonic.t) && (numel(myHarmonic.t) == 1)  %Zero.   %isempty(myHarmonic.t)
            %If time is not assigned - set default time vector.
            myHarmonic = compTimeVector(myHarmonic, myFs/2, 'samples');
        end
        %Set slowly params when time vector was computed.
        for i = find(iterations') %Comp slowly comp of assigned fields.
            slowComp = signalObj(myFs, repmat(mySignalParams.(fNms{i}), size(myHarmonic.t)), myHarmonic.t);
            slowComp = slowComp.compSlowlyChangedComp(slowlyArgs{i}{:});
            mySignalParams.(fNms{i}) = slowComp;
        end
		myHarmonic.signalParams = mySignalParams;
        myHarmonic = resampleToTime(myHarmonic);
    end
    
	%Getters/setters...
	function myFullPhase = getFullPhase(myHarmonic)
		myFullPhase = myHarmonic.fullPhase;
    end
	
% 	function myAmplitude = getAmplitude(myHarmonic)
% 		myAmplitude = myHarmonic.amplitude;
%     end
    
	function mySignalParams = getOrigSignalParams(myHarmonic)
		mySignalParams = myHarmonic.signalParams;
        mySignalParams.t = myHarmonic.t;
    end
    
    %Function returns all original signal parameters in vectors.
	function mySignalParams = getSignalParams(myHarmonic)
		mySignalParams.A = myHarmonic.getParameters('A');
		mySignalParams.f = myHarmonic.getParameters('f');
		mySignalParams.phi = myHarmonic.getParameters('phi');
        mySignalParams.t = myHarmonic.t;
	end
    
	function myHarmonic = setSignalVector(myHarmonic, mySignal)
		myHarmonic.signal = mySignal;
        myHarmonic.fullPhase = compFullPhase(signalObj(myHarmonic.Fs, mySignal, myHarmonic.t));
        myHarmonic.amplitude = envelope(mySignal);
    end
    %...Getters/setters.
    
    %Functions of parameter computation and setting.
	function myHarmonic = setParameters(myHarmonic, varargin)
	%Amplitude, phase, frequency may be pointed like constants or functions.
	%Functions are vectors with length the same as time or only one (a few) periods,
	%which will be replicated on the whole time vector length, constants are too.
        myTime = [];
        myMode = '';
        for i = 1:2:numel(varargin)
            switch varargin{i}
                case 'w'
                    myHarmonic = setParam(myHarmonic, 'f', varargin{i+1}/(2*pi));
                case 'time'
                    myTime = reshape(double(varargin{i+1}), 1, []);
                case 'mode'
                    myMode = varargin{i+1};
                otherwise
                    %Amplitude, frequency, phase.
                    myHarmonic = setParam(myHarmonic, varargin{i}, varargin{i+1});
            end
        end
        if ~isempty(myTime)
            myHarmonic = compTimeVector(myHarmonic, myTime, myMode);
        end
	end
	
	function myVal = getParameters(myHarmonic, myParamName)
		if sum(strcmp(myParamName, {'A', 'f', 'phi'}))
			myVal = myHarmonic.signalParams.(myParamName);
            %Check, does a signal parameter assigned by some signal class.
            if ismethod(myVal, 'getSignalVector')
                myVal = myVal.getSignalVector;
            end
        elseif strcmp(myParamName, 'w')
			myVal = 2*pi*getParameters(myHarmonic, 'f');
        else
            myVal = myHarmonic.(myParamName);
		end
	end
	
	function myHarmonic = compSignalVector(myHarmonic)
        if isempty(myHarmonic.t)
            error('There is no time vector!')
        end
        if isempty(myHarmonic.fullPhase)
            myHarmonic.fullPhase = compFullPhase(myHarmonic);
        end
        if isempty(myHarmonic.amplitude)
            myHarmonic.amplitude = getParameters(myHarmonic, 'A');
        end
        myHarmonic.signal = myHarmonic.amplitude.*sin(myHarmonic.fullPhase);
	end

	function myHarmonic = compTimeVector(myHarmonic, myTime, mode)
		%Set time vector from Fs and pointed length (in seconds (def), in periods, in samples), if time is a scalar, or pointed time vector.
        if ~exist('mode', 'var')
            mode = [];
        end
        if strcmp(mode,  'periods')
            %Translate to seconds.
            myHarmonic = compTimeVector(myHarmonic, myHarmonic.Fs, 'samples'); %For period measurement.
            F0 =  getCentralFreq(myHarmonic);
            T = 1/F0;  %Seconds in the signal period.
            myTime = myTime*T; %Time in pointed periods number.
            mode = 'seconds';
        end
        myHarmonic = compTimeVector@signalObj(myHarmonic, myTime, mode);
        %myHarmonic = resampleToTime(myHarmonic);
	end
	
	function myHarmonic = setNoise(myHarmonic, myNoiseLevel)
		myHarmonic.noiseLevel = myNoiseLevel;
		myHarmonic.signal = awgn(myHarmonic.signal, myNoiseLevel, 'measured', 'linear');
    end
	
	function myFullPhase = compFullPhase(myHarmonic, mySignalParams)
        %Compute from our harmonic params or external struct.
        if ~exist('mySignalParams', 'var')
            mySignalParams = [];
        end
        if ischar(mySignalParams), myFullPhase = compFullPhase@signalObj(myHarmonic, mySignalParams); return; end
        if ~isstruct(mySignalParams) %It's mb need2unwrap flag.
            mySignalParams = getSignalParams(myHarmonic);
        end
        %Full phase mb computed for setted time vector.
        if ~isfield(mySignalParams, 't')
            mySignalParams.t = myHarmonic.t;
        end
        %Freq*time + phase - constant or vector with length the same like time vector.
        myFullPhase = 2*pi*mySignalParams.f.*mySignalParams.t + mySignalParams.phi;
%         %For comparison:
%         if numel(myHarmonic.signal) == numel(myFullPhase) %gag.
%             analytSign = hilbert(myHarmonic.signal);
%             myFullPhase = angle(analytSign);
%             myFullPhase = unwrap(myFullPhase);
%         end
    end
    
    function [signalShifted, myHarmonic] = shiftSignal(myHarmonic, myAngle)
        if ~exist('myAngle', 'var'), myAngle = pi/2;  end
        myHarmonic.fullPhase = []; myHarmonic.phi = myHarmonic.phi+myAngle;
        myHarmonic.signalParams.phi = myHarmonic.signalParams.phi+myAngle;
        myHarmonic = myHarmonic.compSignalVector; signalShifted = myHarmonic.getSignalVector;
    end
    
end

methods(Access = public, Sealed)
    
    function myAmplitude = getAmplitude(myHarmonic)
        myAmplitude = myHarmonic.A;
    end
    
    function myPhase = getPhase(myHarmonic)
        myPhase = myHarmonic.phi;
    end
    
    function myFrequency = getFrequency(myHarmonic)
        myFrequency = myHarmonic.f;
    end
    
    %Function returns all signal parameters in vectors, resampled to time vector.
	function mySignalParams = getHarmonicVectors(myHarmonic)
		mySignalParams.A = myHarmonic.A;
		mySignalParams.f = myHarmonic.f;
		mySignalParams.phi = myHarmonic.phi;
        mySignalParams.t = myHarmonic.t;
	end
    
end


methods(Access = private)
    
	function myHarmonic = setParam(myHarmonic, myParamName, myParamVal)
		if sum(strcmp({'f', 'A', 'phi'}, myParamName))
            myHarmonic.signalParams.(myParamName) = reshape(double(myParamVal), 1, []);
            myHarmonic = resampleToTime(myHarmonic);
		end
    end
    
    %Function implements resampling signal parameters to a new time vector.
    function myHarmonic = resampleToTime(myHarmonic)
			len = length(myHarmonic.t);
			if len  %If time vector is setted.
                mySignParams = getSignalParams(myHarmonic);
                paramNames = fieldnames(mySignParams);
                %Recalculate a signal parameters vectors from original data.
                %Fill a new struct of a signal params, where 
                for i = 1:numel(paramNames)
                    myParamName = paramNames{i};
                    myParamVal = mySignParams.(myParamName);
                    %Tiling by paramVal to the whole time length vector.
                    repeationNumber = floor(len/length(myParamVal));
                    myParamValResampled = repmat(myParamVal, 1, repeationNumber);
                    %The rest part of period.
                    restSamples = len - repeationNumber*length(myParamVal);
                    %myParamValResampled(end:end + restSamples) = myParamValResampled(1:1+restSamples);
                    myParamValResampled = [myParamValResampled myParamVal(1:restSamples)];
                    mySignParams.(myParamName) = myParamValResampled;
                    myHarmonic.(myParamName) = mySignParams.(myParamName);
                end
                myHarmonic.fullPhase = compFullPhase(myHarmonic, mySignParams);
                myHarmonic.amplitude = mySignParams.A;
                myHarmonic = compSignalVector(myHarmonic);
			end
    end
    
    
end
    
    methods (Access = public, Static = true)
        
        function myFullPhase = getSignalAngle(mySignal)
            hilbSign = hilbert(mySignal);
            myFullPhase = atan(hilbSign./mySignal);
        end
        
    end
    

end