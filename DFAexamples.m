function DFAexamples(n,duration,numOfTrials,fs)
% Some examples illustrating DFA vis-a-vis DFA + model selection...
%
% Ton & Daffertshofer, Model selection for identifying power-law scaling
% Neuroimage 136:215-26, 2016, doi:10.1016/j.neuroimage.2016.01.008
%
% See also FluctuationAnalysis, psdfgn
%
%                                               (c) marlow 2016/7
%                                     latest update March 5, 2017
%
% This file is released under the terms of the GNU General Public License,
% version 3. See http://www.gnu.org/licenses/gpl.html

if nargin<4 || isempty(fs), fs=1000; end              % def. fs = 1 kHz
if nargin<3 || isempty(fs), numOfTrials=10; end       % def. # trials = 10
if nargin<2 || isempty(duration), duration=7*60; end  % def. dur. = 7 min
if nargin<1 || isempty(n), n=1:4; end                 % run all examples

numOfSamples=duration*fs;

for i=n
        
    rng('shuffle'); % randomizing the random generator

    %% first generate a random signal and plot its power spectrum
    switch i

        case 1 % Gaussian white noise over the entire range
            
            Q=0.25;
            s=sprintf('Gaussian white noise (strength %g) from %g until %gs',Q,0,duration);
            y=sqrt(Q/fs)*randn(numOfTrials,numOfSamples)';

        case 2 % fractional Gaussian noise over the entire range

            H=0.9;
            
            s=sprintf('single fGn process with H=%g from %g until %gs',H,0,duration);
            y=psdfgn(numOfSamples,numOfTrials,H,0,fs);
            % as an alternative for psdfgn(...) you may use:
            % y=fftfgn(1,H,numOfTrials,numOfSamples,numOfSamples,1)';

        case 3 % two-partite correlations (short and long range)
               %
               % related reading;
               %
               % Collins & De Luca
               % Random walking during quiet standing
               % Phys. Rev. Lett. 73(5):764, 1994
               %
               % Newell, K et al. 
               % Stochastic processes in postural center-of-pressure profiles
               % Exp. Brain. Res. 113:158, 1997

            H=[0.7,0.05];
            bp=[0,1]; 
            
            bp=min(bp,duration); % breakpoints should be smaller than 1/2 of the total duration
            
            s=sprintf('two fGn processes with H=%g from %g to %gs followed by H=%g until %gs',H(1),bp(1),bp(2),H(2),duration);
            y=psdfgn(numOfSamples,numOfTrials,H,bp,fs);
            
        case 4 % fractional Gaussian noise over the entire range with confounding, additive oscillations
               %
               % related reading:
               %
               % Hu, Ivanov, Chen, Carpena, & Stanley
               % Effect of trends on detrended fluctuation analysis
               % Phys. Rev. E 64, 011114, 2001
               %
               % Chen, Ivanov, Hu, & Stanley
               % Effect of nonstationarities on detrended fluctuation analysis
               % Phys. Rev. E 65, 041107, 2002
               %
               % Chen, Hu, Carpena, Bernaola-Galvan, Stanley, & Ivanov
               % Effect of nonlinear filters on detrended fluctuation analysis
               % Phys. Rev. E 71, 011104, 2005

            H=0.7;
            f0=3; % estimate of the CoG natural frequency
            
            s=sprintf('single fGn process with H=%g from %g until %gs + %g Hz oscillations',H,0,duration,f0);
            y=psdfgn(numOfSamples,numOfTrials,H,0,fs);
            t=repmat((0:numOfSamples-1)'/fs,1,numOfTrials);
            y=y+5*sin(2*pi*f0*t+2*pi*repmat(rand(1,numOfTrials),numOfSamples,1)); % the phases of the additive oscillations are randomized over trials

    end
    
    fprintf('\nExample %d: %s\n',i,s);
    
    figure(i);
    set(gcf,'Name',s,'Position',[50*i,300-50*i,1024,256]);
    
    %% plot 'plain' power spectrum
    subplot(1,4,1);
    pspectrum(y,fs); % this is just for plotting the power spectrum
    title('power spectrum')

    %% convert fGn into fBm (or Gaussian noise into Brownian motion if H=0.5
    y=cumsum(y);

    %% define the total range of analysis
    maxNumOfSegments=400; % use values >= 1000 if you have sufficiently
                          % many data and a fast computer
    maxNumOfSegmentsPerRange=200;   % use values >= 500 if ... ;-)
    fullrange=[10/1000,duration/4]; % 10 ms until 1/4 of the recording time
    interval2evaluate=[fullrange,maxNumOfSegments];
    
    %% conventional DFA over the entire range - this suits examples 1 & 2
    subplot(1,4,2);
    % define the fitting range
    interval2fit=interval2evaluate;
    % run the analysis (including text & graphics report)
    fluctuationAnalysis(y,interval2evaluate,interval2fit,'DFA','verbose','graphics','fs',fs);
    title('DFA over full range'); set(gca,'xlim',fullrange); drawnow;

    %% conventional DFA over two distinct ranges - this suits examples 3 & 4
    subplot(1,4,3);
    % define the fitting ranges
    range1=[20/1000,800/1000];     % 20 ms until 800 ms
    range2=[1200/1000,duration/4]; % 1200 ms until 1/4 of the recording time
    interval2fit={ [range1,maxNumOfSegmentsPerRange], [range2,maxNumOfSegmentsPerRange] };
    % run the analysis (including text & graphics report)
    fluctuationAnalysis(y,interval2evaluate,interval2fit,'DFA','verbose','graphics','fs',fs);
    title('DFA over two ranges'); set(gca,'xlim',fullrange); drawnow;

    %% DFA including (Bayesian) non-parametric model selection - this suits all examples
    subplot(1,4,4);
    % define the fitting range
    interval2fit=interval2evaluate;
    % run the analysis (including text & graphics report)
    fluctuationAnalysis(y,interval2evaluate,interval2fit,'DFA+','verbose','graphics','fs',fs);
    title('DFA+ over full range'); set(gca,'xlim',fullrange); drawnow;

end

%% helper function for plotting the fft-based power spectrum
function [p,f]=pspectrum(y,fs)
% estimate the 'plain' power specturm using fft
N=size(y,1);
p=abs(fft(y/N)).^2;
nfft=size(y,1);
if rem(nfft,2) % nfft odd
    select = (1:(nfft+1)/2)';
else
    select = (1:nfft/2+1)'; % include DC AND Nyquist
end
% select the first half...
p=p(select,:);
f=(select-1)/(numel(select)-1)*(fs/2);

for k=1:size(p,2) % normalize every power spectrum
    p(:,k)=p(:,k)/trapz(f,p(:,k));
end
% compute the mean over trials
p=mean(p,2,'omitnan');

if nargout==0
    
    loglog(f,mean(p,2),'linewidth',2);

    xmin=floor(log10(min(f(f>0))));
    xmax=ceil(log10(max(f)));
    set(gca,'xlim',10.^[xmin,xmax],'xtick',10.^(xmin:1+((xmax-xmin)>5):xmax));

    grid on;
    if fs~=1, xlabel('f [Hz]'); else, xlabel('1/n'); end
    ylabel('p');

end

