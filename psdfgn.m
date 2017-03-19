function [y,p,f,py,fy]=psdfgn(N,M,H,bp,fs,verbose)
% x=psdfgn(N,M,H,fs) generating a zscored, real-valued fGn signal given
% a set of Hurst exponents (plus breakpoints)
%
% Input:
% - N = number of samples
% - M = number of trials
% - H = Hurst exponent(s); see below for examples
% - bp = breakpoint(s); see below for examples
% - fs = samping rate (optional, default = 1)
% - verbose = if true => plot power spectra (optional, default = false)
%
% Output:
% - y = real-valued NxM array
% - p = generating power spectrum (Nx1 array)
% - f = frequency axis corresponding to p
% - py = resulting (mean) power spectrum (Nx1 array)
% - fy = frequency axis corresponding to py
%
% Examples:
%
%   1 -- fractional Gaussian noise over the entire range
%
%     H=0.9; % persistent correlations
%     bp=0; % from 0s onwards
%     [y,p,f,py,fy]=psdfgn(numOfSamples,numOfTrials,H,bp,fs);
%
%   2 -- two-partite correlations (short and long range)
%
%     breakpoint=1; % 1s
%     H=[0.9,0.1]; % persistent correlations, followed by anti-persistent
%     bp=[0,1]; % H=0.9 from 0 to 1s, after that H=0.1 
%     [y,p,f,py,fy]=psdfgn(numOfSamples,numOfTrials,H,bp,fs);
%
%   3 -- three-partite correlations (short, middle, and long range)
%
%     H=[0.9,0.4,0.1];
%     bp=[0,1,10];
%     [y,p,f,py,fy]=psdfgn(numOfSamples,numOfTrials,H,bp,fs);
%
% NB: use y=psdfgn(...) and y=cumsum(y) to generate a fBm signal.
%
% See also psd2signal, randomizeFourierPhase, fftfgn
%
%                                                     (c) marlow 2016
%
% This file is released under the terms of the GNU General Public License,
% version 3. See http://www.gnu.org/licenses/gpl.html

if nargin<6, verbose=false; end 
if nargin<4, fs=1; end % set default sampling rate = 1

% set frequency axis
f=(0:N-1)'/(N-1)*(fs/2);

% extract breakpoints
if nargin<3 || isempty(bp)
    bp=zeros(size(H));
end

beta=nan(size(f));
p0=nan(size(f));
S0=1;
B0=0;

% construct the exponents of the power spectrum
for k=1:numel(bp)

    B=2*H(k)-1; % fGn process (H=1 -> 1/f;   H=1/2 -> p=const; H=0 -> p=f)
%   B=2*H(k)+1; % fBm process (H=1 -> 1/f^3; H=1/2 -> p=1/f^2; H=0 -> p=1/f)
    f0=1/bp(k);

    if ~isfinite(f0), f0=max(f); end
    i=f<=f0;
    beta(i)=B;
    p0(i)=S0*f0^(B-B0);
    S0=p0(i(1));
    B0=B;
    
end

p=p0.*f.^-beta; % make a 1/f power spectrum
p(1)=0; % remove DC (or NaNs))
p=p/trapz(f,p); % normalize the power

% generate a zscored, real-valued signal with power spectrum p
y=psd2signal(p,N,fs,'onesided');
% duplicate the signal M times
y=repmat(y,1,M);
% randomize the Fourier phase (keeps power intact)
y=randomizeFourierPhase(y);

if ~nargout || verbose || nargout>3 % this if for plotting or for extra output

    % estimate the 'plain' power specturm using fft
    py=abs(fft(y/N)).^2;
    nfft=size(y,1);
    if rem(nfft,2) % nfft odd
        select = (1:(nfft+1)/2)';
    else
        select = (1:nfft/2+1)'; % include DC AND Nyquist
    end
    % select the first half...
    py=py(select,:);
    fy=(select-1)/(numel(select)-1)*(fs/2);

    for k=1:size(py,2) % normalize every power spectrum
        py(:,k)=py(:,k)/trapz(fy,py(:,k));
    end
    % compute the mean over trials
    py=mean(py,2,'omitnan');

end

%% generate a plot if no output is wanted
if ~nargout || verbose

    loglog(f,p,fy,py,'linewidth',2);
    xmin=floor(log10(min([f(f>0);fy(fy>0)])));
    xmax=ceil(log10(max([f(:);fy(:)])));
    set(gca,'xlim',10.^[xmin,xmax],'xtick',10.^(xmin:1+((xmax-xmin)>5):xmax));
    grid on;
    legend('source power','generated power','Location','South');
    if fs~=1, xlabel('f [Hz]'); else, xlabel('1/n'); end
    ylabel('p');

end




