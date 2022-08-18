% Script to run models

%%%

%%author__ = Shanna Kulik & Prejaas Tewarie
%%contact__ = l.douw@amsterdamumc.nl
%%date__ = 2020
%%status__ = finished

%%%%%%%%%%%%%%%%%%%%
% Review History   %
%%%%%%%%%%%%%%%%%%%%

% Reviewed by Eduarda Centeno 2022/6

%%%%%%%%%%%%%%%%%%%%
% Requirements     %
%%%%%%%%%%%%%%%%%%%%

% Created and ran with Matlab R2018b

%%% Other m-files required:
% rangenorm.m
% Jansen_network_RK2.m
% fft_filt.m
% pli_matteo.m
% phaselockvalue_leak.m
% David Ferreira (2021). Quantiles (https://www.mathworks.com/matlabcentral/fileexchange/70279-quantiles), MATLAB Central File Exchange. Retrieved December 9, 2021.

%%% Toolboxes 
% Brain Connectivity Toolbox https://sites.google.com/site/bctnet/ 
% version: 2017-15-01

% CONTEST: A Controllable Test Matrix Toolbox for MATLAB. https://www.maths.ed.ac.uk/~dhigham/Publications/P83.pdf
% version: 1.2 | 20-May-2008

% octave-networks-toolbox: http://dx.doi.org/10.5281/zenodo.22398
% version: second release, as of August 2 2015

% EEGLAB https://eeglab.org/
% version: 14.0.0


%% Parameters

% MEG
sample = 1250;                              % sample frequency
f_low = 8;                                      % bandpass filter
f_high = 13;                                    % bandpass filter
N = 78;                                         % network size
win = floor(linspace(2,14,13).*sample);         % define window widths
sm = floor(sample.*0.25);                        % smoothing parameter

% simulation parameters NMM Jansen Rit
h = 0.0001;                               % integration time step
T = 20;                                   % observation time in sec
noise = 0.05;                             % noise parameter
P = 150;                                  % external input
Nsize = N;                                % network size
coupling = 0.1:0.012:0.3;

%% Get network

subs = {'asub'} 
for sub=1:numel(subs)
    namos = char(strcat(['/path/to/sub/connectome.aal.txt']));
    sc_unsorted = load(namos);
    sc = sc_unsorted(1:78,1:78);
    sc(:,:,sub) = (sc + sc').*~eye(78);
end

mean_sc = squeeze(mean(sc,3)); 

y=sort(mean_sc);
Q(1) = nanmedian(y(find(y<nanmedian(y))));
Q(3) = nanmedian(y(find(y>nanmedian(y))));
IQR = Q(3)-Q(1);
thres = Q(3) + 1.5*IQR;
index3 = find(mean_sc>thres);
mean_sc(index3) = Q(3) + 1.5*IQR;
sc_norm = rangenorm(mean_sc);

vel = 10;
load('/path/to/euclid_AAL.mat');
delay = euclid_AAL/vel;


%% simulations Jansen Rit model
params = zeros([0, 3]);   
npairs = size(params, 1);

tPLV_idx=zeros(npairs, Nsize,Nsize);
tPLI_idx=zeros(npairs, Nsize,Nsize);
tAEC_idx=zeros(npairs, Nsize,Nsize);
peak_freq=zeros(npairs);

parfor idx=1:npairs
    disp(['idx = ' num2str(idx)])
    subs_idx = params(idx,1);
    cp_idx = params(idx,3);
        
    PLIs=zeros(Nsize,Nsize);
    PLVs=zeros(Nsize,Nsize);
    AECs=zeros(Nsize,Nsize);

    % run NMM
    timeseries = Jansen_network_RK2(h,T,noise, P, cp_idx, Nsize,sc_norm,delay);
    
    % downsample
    down = 10;
    h2 = h * down;
    times = timeseries(:,1:down:end);
    data_filt_nmm = fft_filt(times',1/h2,f_low,f_high,N);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute functional connectivity
    PLIs(:,:) = pli_matteo(data_filt_nmm);
    PLVs(:,:) = phaselockvalue(data_filt_nmm);
    AECs(:,:) = corr(abs(hilbert(data_filt_nmm))).*~eye(N);
    
    % Power spectrum
    data_stat = timeseries(1,:);
    F = length(data_stat);
    Fs = 1/h;
    xdft = fft(data_stat);
    xdft = xdft(1:F/2+1);
    psdx = (1/(Fs*F)) * abs(xdft).^2;
    freq = 0:Fs/length(data_stat):Fs/2;
    freq(1:2)=[];
    psdx(1:2)=[];
    
%     Calculate peak frequency
    peak = max(psdx);
    peak_freq(idx) = freq(psdx == peak);
    
    tPLI_idx(idx,:,:)=PLIs;
    tPLV_idx(idx,:,:)=PLVs;
    tAEC_idx(idx,:,:)=AECs;

end
toc

%% final bookkeeping

PLI_JR = zeros(numel(runs),numel(coupling), Nsize,Nsize);
PLV_JR = zeros(numel(runs),numel(coupling), Nsize,Nsize);
AEC_JR = zeros(numel(runs),numel(coupling),Nsize,Nsize);
peak_freq_JR = zeros(numel(runs),numel(coupling)); 
pf = peak_freq(:,1);

idx = 1;

for sub=1
    for  it = 1:numel(runs)
        for k = 1:numel(coupling) 
            PLI_JR(it,k,:,:) = squeeze(tPLI_idx(idx,:,:));
            PLV_JR(it,k,:,:) = squeeze(tPLV_idx(idx,:,:));
            AEC_JR(it,k,:,:) = squeeze(tAEC_idx(idx,:,:));
            peak_freq_JR(it,k) = squeeze(pf(idx));
            
            idx = idx + 1;
        end
    end
end