function [Fn,conFn]=detrendedFluctuationAnalysis(x,n,varargin)
%[Fn,sFn,n] = detrendedFluctuationAnalysis(x,n,varargin)
%
% Performs the detrended fluctuation analysis of an fGn process
%
% Input: 
% - x = fGn time series (column)
% - n = vector containing sub-periods (column)
% - 'order'   followed by an integer: set order of polynomial used for
%             detrending (optional; default=1 = linear detrending)
% - 'verbose' flag to add a text report
%             (default = true if nargout==0, otherwise false = no report)
%
% Output:
% - Fn = divergence curve as a function of n
% - conFn = corresponding confidence intervals as a function of n presuming
%           that Fn are normally distributed
%
% See also fluctuationAnalysis, diffusionAnalysis
%
% Reference: Peng et al. (1994), Phys. Rev. E49:1685
%
%                                              (c) marlow 2012-15
%                                       latest update 11 jan 2015
%
% This file is released under the terms of the GNU General Public License,
% version 3. See http://www.gnu.org/licenses/gpl.html

%% set defaults and check variable input
order=1; % default is linear detrending
verbose=false; % default is no progress report
if numel(varargin)
    oi=find(strncmpi(varargin,'ord',3));
    if ~isempty(oi), order=varargin{oi+1}; end
    verbose=sum(strncmpi(varargin,'ver',3))~=0;
end

%% define the confidence interval
N=size(x,1);
Fn=zeros(size(n));
conFn=squeeze(zeros([size(n),2]));

conf=95/100; % use a 95% confidence interval
if exist('icdf','file')
    conf=icdf('Normal',1-(1-conf)/2,0,1);
else
    conf=1.96;
end

%% loop over segments
if verbose, fprintf('%s: looping %d segments ',mfilename,numel(n)); end

for i=1:numel(n)

    if verbose && mod(i-1,round(numel(n)/10))==0, fprintf('.'); end
    
    % get the max order for detrending
    o=min(order,n(i)-1);
    
    % get the number of segments per trila
    a=floor(N/n(i));
    sq=zeros(a,size(x,2));
    
    for k=1:size(x,2) % loop over trials
        
        % make the matrix of non-overlapping segments
        X=reshape(x(1:a*n(i),k),n(i),a);

        % detrend the segments
        switch o
            
            case 0 % DC-removal
                X=detrend(X,'constant');
            case 1 % linear detrending
                X=detrend(X,'linear');
            otherwise % nonlinear detrending
                ind=detrend(1:size(X,1),'constant')';
                for j=1:size(X,2)
                    c=polyfit(ind,X(:,j),o);
                    X(:,j)=X(:,j)-polyval(c,ind);
                end
                
        end
        
        % get the mean squares per segments
        sq(:,k)=mean(X.^2,1,'omitnan');
        
    end
    
    % average over segments & trials and get the standard error
    F=mean(sq(:),'omitnan');
    seF=std(sq(:),'omitnan')/sqrt(numel(sq));
    
    % store mean and confidence interval
    Fn(i)=F;
    conFn(i,:)=max([F-conf*seF,F+conf*seF],0);
    
end
if verbose, fprintf('\n'); end

%% compute the square root of the mean & confidence interval
Fn=sqrt(Fn);
conFn=sqrt(conFn);