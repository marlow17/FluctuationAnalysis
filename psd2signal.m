function x=psd2signal(p,N,fs,range)
% x=psd2signal(p,N,fs,range) generating a zscored, real-valued signal given
% a certain power spectrum base on the Wiener-Khinchin theorem
%
% Input:
% - p = real-valued vector containing the power
% - N = number of samples (optional, default numel(p) if p is one-sided),
%                          numel(9)/2 otherwise)
% - fs = samping rate (optional, default = 1)
% - range = 'onesided' (default) or 'twosided'
%
% Output:
% - x = real-valued vector of length N
%
% See also ifft, psdfgn
%
%                                                     (c) marlow 2016
%
% This file is released under the terms of the GNU General Public License,
% version 3. See http://www.gnu.org/licenses/gpl.html

if nargin<4 || isempty(range), range='onesided'; end
if nargin<3 || isempty(fs), fs=1; end
if nargin<2 || isempty(N)
    if strcmp(range,'onesided')
        N=numel(p);
    else
        N=numel(p)/2;
    end
end

p=p(:);

% dublicate the power spectrum if only one side is given...
if strcmp(range,'onesided')
      p=[p;p(end-(1-rem(N,2)):-1:2)];
end

M=numel(p);
% inverse Fourier transforming the square root of the spectral power ...
x=ifft(sqrt(p*M*fs),M,'symmetric');
% remove dublicates
x=[x(1);(x(2:N)+x(end:-1:end-N+2))];

% z-scoring the output ...
x=zscore(x);

