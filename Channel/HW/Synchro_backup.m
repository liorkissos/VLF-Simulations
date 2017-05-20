function   Synchro( src,event,varargin )
%SYNCHRO Summary of this function goes here
%   Detailed explanation goes here

global BBB
global CCC
global peaks_actual
persistent nn
global timeout

if isempty(nn)
    % disp('1st iteration')
    nn=0;
end

prominence_ref=varargin{1};
Hpost_flt=varargin{2};
Pref=varargin{3};
Width=varargin{4};
h=varargin{5};


x=event.Data; % the data read from A/D

%% post sampling filtering

z=conv(Hpost_flt.numerator,x); % filtering out the aliased replicas that might harm the chirp detection

%% Chirp convolution
y=conv(h,z);

CCC=[CCC;y]; %  debugging variable


%%% NEW-start
Pseq=y'*y;
prominence=prominence_ref*sqrt(Pseq/Pref);

%prominence=prominence_ref; %%% TEMP!!!!!

% prominence1=prominence*sqrt(Pseq/1.5e8);
% prominence1=prominence*sqrt(Pseq/7.1e8);
% prominence1=prominence*sqrt(Pseq/1.3e9);
%prominence1=prominence*sqrt(Pseq/3.7e8);

%pks=findpeaks(y,'MinPeakProminence',prominence1,'MinPeakWidth',3);
%pks2=findpeaks(y,'MinPeakProminence',prominence,'MinPeakHeight',3.2e3);

pks=findpeaks(y,'MinPeakProminence',prominence,'MinPeakWidth',Width);
%%% NEW- end

%%% OLD-start

%pks=findpeaks(y,'MinPeakProminence',prominence,'MinPeakWidth',3);

%%% OLD- end

if numel(pks)==1
    nn=nn+1;
    if nn==1
        peaks_actual=pks;
    end
elseif numel(pks)>1
    error('too many peaks')
    src.stop()
end

%nn

if nn==1
    BBB = [BBB;x];
elseif nn==2
    BBB = [BBB;x];
    peaks_actual=[peaks_actual;pks];
    clear persistent nn
    % pause(2)
    src.stop()
    % pause(2)
    disp('packet identified. exiting')
    
end

if toc>timeout
    src.stop()
    disp('timeout expired. exiting')
end





