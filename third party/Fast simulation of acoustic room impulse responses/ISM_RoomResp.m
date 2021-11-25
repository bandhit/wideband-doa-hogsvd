function [RIRvec] = ISM_RoomResp(Fs,beta,rt_type,rt_val,X_src,X_rcv,room,varargin)
%ISM_RoomResp  RIR based on Lehmann & Johansson's image-source method
%
% RIR = ISM_RoomResp(Fs,BETA,RT_TYPE,RT_VAL,SOURCE,SENSOR,ROOM)
% RIR = ISM_RoomResp( ... ,'arg1',val1,'arg2',val2,...)
% 
% This function generates the room impulse response (RIR) between a sound
% source and an acoustic sensor, based on various environmental parameters
% such as source and sensor positions, enclosure's dimensions and
% reflection coefficients, etc., according to Lehmann and Johansson's
% implementation of the image-source method (see below). The input
% parameters are defined as follows:
%
%       Fs: scalar, sampling frequency (in Hz). Eg: 8000.
%     BETA: 1-by-6 vector, corresponding to each wall's reflection 
%           coefficient: [x1 x2 y1 y2 z1 z2]. Index 1 indicates wall closest
%           to the origin. This function assumes strictly non-negative BETA
%           coefficients. Set to [0 0 0 0 0 0] to obtain anechoic response 
%           (direct path only), in which case the value of RT_VAL is 
%           discarded. E.g.: [0.75 0.75 0.8 0.25 0.3 0.9].
%  RT_TYPE: character string, measure of reverberation time used for the 
%           definition of the coefficients in BETA. Set to either 'T60' or
%           'T20'. 
%   RT_VAL: scalar, value of the reverberation time (in seconds) defined by
%           RT_TYPE. Set to 0 to obtain anechoic response (same effect as 
%           setting BETA to [0 0 0 0 0 0]), in which case the BETA 
%           coefficients are discarded. E.g.: 0.25.
%   SOURCE: 1-by-3 vector, indicating the location of the source in space 
%           (in m): [x y z]. E.g.: [1 1 1.5].
%   SENSOR: 1-by-3 vector, indicating the location of the microphone in 
%           space (in m): [x y z]. E.g.: [2 2 1.5].
%     ROOM: 1-by-3 vector, indicating the rectangular room dimensions 
%           (in m): [x_length y_length z_length]. E.g.: [4 4 3].
%
% In addition, a number of other (optional) parameters can be set using a 
% series of 'argument'--value pairs. The following parameters (arguments)
% can be used:
%
%          'c': scalar, speed of acoustic waves (in m/s). Defaults to 343.
%   'Delta_dB': scalar (in dB), parameter determining how much the resulting 
%               impulse response is cropped (i.e. RIR length): the impulse
%               response is computed until the time index where its overall
%               energy content has decreased by 'Delta_dB' decibels, after
%               which the computations stop. Not relevant if BETA=zeros(1,6). 
%               Defaults to 50.
%   'MaxDelay': scalar (in seconds), defines the desired length for the 
%               computed RIR. If defined as non-empty, this parameter 
%               overrides the setting of 'Delta_dB'. Use 'MaxDelay' if the
%               RIR length is known exactly prior to its computation. Not
%               relevant if BETA=zeros(1,6). Defaults to [].
% 'SilentFlag': set to 1 to disable all on-screen messages from this
%               function. Defaults to 0.
% 
% This function returns the time coefficients of the filter (transfer
% function) in the parameter RIR. The filter coefficients are real and
% non-normalised. The first value in the vector RIR, i.e., RIR(1),
% corresponds to time t=0. The number of coefficients returned is variable
% and results from the value of 'Delta_dB' defined by the user: the filter
% length will be as large as necessary to capture all the relevant
% highest-order reflections. 
% 
% This implementation uses Lehmann and Johansson's variant (see "Prediction
% of energy decay in room impulse responses simulated with an image-source
% model", J. Acoust. Soc. Am., vol. 124(1), pp. 269-277, July 2008) of
% Allen & Berkley's "Image Method for Efficiently Simulating Small-room
% Acoustics" (J. Acoust. Soc. Am., vol. 65(4), April 1979). This function
% implements a phase inversion for each sound reflection off the room's
% boundaries, which leads to more accurate room impulse responses (when
% compared to RIRs recorded in real acoustic environments). Also, the
% computations make use of fractional delay filters, which allow the
% representation of non-integer delays for all acoustic reflections.

% Release date: March 2012
% Author: Eric A. Lehmann, Perth, Australia (www.eric-lehmann.com)

% Explanations for the following code -------------------------------------
% This implementation of the image method principle has been speficically
% optimised for execution speed. The following code is based on the
% observation that for a particular dimension, the delays from the image
% sources to the receiver increases monotonically as the absolute value of
% the image index (m, n, or l) increases. Hence, all image sources whose
% indices are above or equal to a specific limit index (for which the
% received delay is above the relevant cut-off value) can be discarded. The
% following code checks, for each dimension, the delay of each received
% path and automatically determines when to stop, thus avoiding unnecessary
% computations (the amount of TF cropped depends on the 'Delta_dB' 
% parameter).
% The resulting number of considered image sources hence automatically
% results from environmental factors, such as the room dimensions, the
% source and sensor positions, and the walls' reflection coefficients. As a
% result, the length of the computed transfer function has an optimally
% minimum length (no extra padding with negligibly small values).
%--------------------------------------------------------------------------

VarList = {'SilentFlag'         0;       % set to 1 to disable on-screen messages
           'c'                  343;	 % sound propagation speed
           'Delta_dB'           50;      % attenuation limit
           'MaxDelay'           []};     % predefined RIR length, overrides value of 'Delta_dB'
eval(SetUserVars(VarList,varargin));

global RIRvec TimePoints           % not too pretty, but this avoids passing potentially large 
                                   % vectors to frequently called subfunctions...

%-=:=- Check user input:
if X_rcv(1)>=room(1) || X_rcv(2)>=room(2) || X_rcv(3)>=room(3) || X_rcv(1)<=0 || X_rcv(2)<=0 || X_rcv(3)<=0,
    error('Receiver must be within the room boundaries!');
elseif X_src(1)>=room(1) || X_src(2)>=room(2) || X_src(3)>=room(3) || X_src(1)<=0 || X_src(2)<=0 || X_src(3)<=0,
    error('Source must be within the room boundaries!');
elseif ~isempty(find(beta>=1,1)) || ~isempty(find(beta<0,1)),
    error('Parameter ''BETA'' must be in the range [0...1).');
end

beta = -abs(beta);      % implement phase inversion in Lehmann & Johansson's ISM implementation

X_src = X_src(:);       % Source location
X_rcv = X_rcv(:);       % Receiver location
beta = beta(:);         % Reflection coefficients
Rr = 2*room(:);         % Room dimensions

%-=:=- Calculate maximum time lag to consider in RIR -=:=-
if ~isequal(beta(:),zeros(6,1)) && rt_val~=0,      % non-anechoic case: compute RIR's decay time necessary to reach 
    if isempty(MaxDelay),          % Delta_dB (using Lehmann & Johansson's EDC approximation method)
        MaxDelay = ISM_RIR_DecayTime(Delta_dB,rt_type,rt_val,1-beta.^2,room,X_src,X_rcv,Fs,c);
    end
else                               % Anechoic case: allow for 5 times direct path in TF
    DPdel = norm(X_rcv - X_src)/c; % direct path delay in [s]
    MaxDelay = 5*DPdel;
	beta = zeros(6,1);		% in case rt_val=0 only
end
TForder = ceil(MaxDelay*Fs);       % total length of RIR [samp] to reach Delta_dB

TimePoints = ([0:TForder-1]/Fs).';
RIRvec = zeros(TForder,1);

%-=:=- Summation over room dimensions:
if ~SilentFlag, fprintf(' [ISM_RoomResp] Computing transfer function '); end;
for a = 0:1
    for b = 0:1
        for d = 0:1
            if ~SilentFlag, fprintf('.'); end;
            
            m = 1;              % Check delay values for m=1 and above
            FoundLValBelowLim = Check_lDim(a,b,d,m,X_rcv,X_src,Rr,c,MaxDelay,beta,Fs);
            while FoundLValBelowLim==1,
                m = m+1;
                FoundLValBelowLim = Check_lDim(a,b,d,m,X_rcv,X_src,Rr,c,MaxDelay,beta,Fs);
            end
            
            m = 0;              % Check delay values for m=0 and below
            FoundLValBelowLim = Check_lDim(a,b,d,m,X_rcv,X_src,Rr,c,MaxDelay,beta,Fs);
            while FoundLValBelowLim==1,
                m = m-1;
                FoundLValBelowLim = Check_lDim(a,b,d,m,X_rcv,X_src,Rr,c,MaxDelay,beta,Fs);
            end

        end
    end
end
if ~SilentFlag, fprintf('\n'); end;


%============
function [FoundLValBelowLim] = Check_lDim(a,b,d,m,X_rcv,X_src,Rr,c,MaxDelay,beta,Fs)

FoundLValBelowLim = 0;

l = 1;              % Check delay values for l=1 and above
FoundNValBelowLim = Check_nDim(a,b,d,l,m,X_rcv,X_src,Rr,c,MaxDelay,beta,Fs);
while FoundNValBelowLim==1,
    l = l+1;
    FoundNValBelowLim = Check_nDim(a,b,d,l,m,X_rcv,X_src,Rr,c,MaxDelay,beta,Fs);
end
if l~=1, FoundLValBelowLim = 1; end;

l = 0;              % Check delay values for l=0 and below
FoundNValBelowLim = Check_nDim(a,b,d,l,m,X_rcv,X_src,Rr,c,MaxDelay,beta,Fs);
while FoundNValBelowLim==1,
    l = l-1;
    FoundNValBelowLim = Check_nDim(a,b,d,l,m,X_rcv,X_src,Rr,c,MaxDelay,beta,Fs);
end
if l~=0, FoundLValBelowLim = 1; end;


%============
function [FoundNValBelowLim] = Check_nDim(a,b,d,l,m,X_rcv,X_src,Rr,c,MaxDelay,beta,Fs)

global RIRvec TimePoints

FoundNValBelowLim = 0;

n = 1;          % Check delay values for n=1 and above
dist = norm( [2*a-1; 2*b-1; 2*d-1].*X_src + X_rcv - Rr.*[n;l;m] );
foo_time = dist/c;
while foo_time<=MaxDelay,    % if delay is below TF length limit for n=1, check n=2,3,4...
    foo_amplitude = prod(beta.^abs([n-a; n; l-b; l; m-d; m])) / (4*pi*dist);
    RIRvec = RIRvec + foo_amplitude * sinc((TimePoints-foo_time)*Fs);
    n = n+1;
    dist = norm( [2*a-1; 2*b-1; 2*d-1].*X_src + X_rcv - Rr.*[n;l;m] );
    foo_time = dist/c;
end
if n~=1, FoundNValBelowLim = 1; end;

n = 0;          % Check delay values for n=0 and below
dist = norm( [2*a-1; 2*b-1; 2*d-1].*X_src + X_rcv - Rr.*[n;l;m] );
foo_time = dist/c;
while foo_time<=MaxDelay,    % if delay is below TF length for n=0, check n=-1,-2,-3...
    foo_amplitude = prod(beta.^abs([n-a; n; l-b; l; m-d; m])) / (4*pi*dist);
    RIRvec = RIRvec + foo_amplitude * sinc((TimePoints-foo_time)*Fs);
    n = n-1;
    dist = norm( [2*a-1; 2*b-1; 2*d-1].*X_src + X_rcv - Rr.*[n;l;m] );
    foo_time = dist/c;
end
if n~=0, FoundNValBelowLim = 1; end;
