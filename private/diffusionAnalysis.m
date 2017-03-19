function [Sn,conSn]=diffusionAnalysis(x,n,varargin)
%[Sn,conSn,n] = DiffusionAnalysis(x,n,varargin)
%
% Performs a diffusion analysis of a stochastic process
%
% Input: 
% - x = random time series (column)
% - n = vector containing sub-periods (column)
% - 'verbose' option flag to add a text report
%             (default = true if nargout==0, otherwise false = no report)
%
% Output:
% - Fn = mean squared displacement as a function of n
% - conFn = corresponding confidence intervals as a function of n presuming
%           that Fn are normally distributed
%
% See also fluctuationAnalysis, detrendedFluctuationAnalysis
%
%                                              (c) marlow 2012-16
%                                     latest update June 20, 2016
%
% This file is released under the terms of the GNU General Public License,
% version 3. See http://www.gnu.org/licenses/gpl.html

%% set defaults and check variable input
verbose=false; % default is no progress report
if numel(varargin)
    verbose=sum(strncmpi(varargin,'ver',3))~=0;
end

%% define the confidence interval
Sn=zeros(size(n));
conSn=squeeze(zeros([size(n),2]));

conf=95/100; % use a 95% confidence interval
if exist('icdf','file')
    conf=icdf('Normal',1-(1-conf)/2,0,1);
else
    conf=1.96;
end

%% loop over time lags
if verbose, fprintf('%s: looping %d segments ',mfilename,numel(n)); end

for i=1:numel(n)
    
    if verbose && mod(i-1,round(numel(n)/10))==0, fprintf('.'); end

    N=size(x,1)-n(i)+1;   
    sq=(x(1:N,:)-x((1:N)-1+n(i),:)).^2;    
    me=mean(sq(:),'omitnan');
    se=std(sq(:),'omitnan')/sqrt(numel(sq));
    
    % compute the mean diffusion & the confidence interval
    Sn(i)=me;
    conSn(i,:)=[me-conf*se,me+conf*se];
    
end

if verbose, fprintf('\n'); end

%% compute the confidence interval
conSn=max(conSn,0);
