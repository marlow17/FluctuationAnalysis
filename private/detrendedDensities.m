function [dens,xi,expect,conf]=detrendedDensities(X,n,N,varargin)
%[dens,xi,expect,conf]=detrendedDensities(X,n,N,order)
%
% Performs a detrending of an fGn process and estimates the probability
% density of the fluctuations as a function of segment length n
%
% Input:
% - X = fGn time series (column)
% - n = vector containing sub-periods (column)
% - N = number of bins when prob. densities (optional, default=100)
% - 'order'   followed by an integer: set order of polynomial used for
%             detrending (optional; default=1 = linear detrending)
% - 'verbose' flag to add a text report
%             (default = true if nargout==0, otherwise false = no report)
%
% Ton & Daffertshofer, Model selection for identifying power-law scaling
% Neuroimage 136:215-26, 2016, doi: 10.1016/j.neuroimage.2016.01.008
%
% See also ksdensity
%                                              (c) marlow 2012-16
%                                                 (c) robert 2015
%                                     latest update June 19, 2016
%
% This file is released under the terms of the GNU General Public License,
% version 3. See http://www.gnu.org/licenses/gpl.html

%% set defaults and check variable input
order=1; % default is linear detrending
verbose=false; % default is no progress report
if nargin<3, N=100; end % default number of bins
if numel(varargin)
    oi=find(strncmpi(varargin,'ord',3));
    if ~isempty(oi), order=varargin{oi+1}; end
    verbose=sum(strncmpi(varargin,'ver',3))~=0;
end

xi=nan(numel(n),N);
dens=xi;
U=nan(numel(n),1);
expect=U;

%% loop over segments
if verbose, fprintf('%s: estimating %d densities ',mfilename,numel(n)); end

for i=1:numel(n)
    
    if verbose && mod(i-1,round(numel(n)/10))==0, fprintf('.'); end
    
    % get the max order for detrending
    o=min(order,n(i)-1);
    
    % get the number of segments per trial
    a=floor(size(X,1)/n(i));
    sq=zeros(a,size(X,2));
    
    for k=1:size(X,2) % loop over trials
        
        % make the matrix of non-overlapping segments
        Y=reshape(X(1:a*n(i),k),n(i),a);
        % detrend the segments
        switch o
            
            case 0 % DC-removal
                Y=detrend(Y,'constant');
            case 1 % linear detrending
                Y=detrend(Y,'linear');
            otherwise % nonlinear detrending
                ind=detrend(1:size(Y,1),'constant')';
                for j=1:size(Y,2)
                    c=polyfit(ind,Y(:,j),o);
                    Y(:,j)=Y(:,j)-polyval(c,ind);
                end
                
        end
        % get the mean squares per segment
        sq(:,k)=mean(Y.^2);
        
    end % end loop over trials
    
    sq=sq(:); % consider consecutive segments & trials likewise independent
    
    % In case of zeros in sq(:), make them eps to avoid errors in ksdensity
    sq(sq==0)=eps;
    
    % estimate density by kernel density estimation and apply
    % logarithmically spaced positive support axis
    [~,Fdxi]=ksdensity(sq,'support','positive');
    
    % here we previously used an interpolation to get a log-spaced bin-axis
    % but the recent version of ksdensity supports 'oversampled' bins,
    % hence we just redo ksdensity with an adjust axis
    Fdxi=logspace(log10(min(Fdxi)),log10(max(Fdxi)),N);
    [dens(i,:),xi(i,:),U(i)]=ksdensity(sq,Fdxi,'support','positive');
    
end

if verbose, fprintf('\n'); end

%% rescale density for change in variable g(x) = sqrt(x)
xi=sqrt(xi);
dens=2*xi.*dens; %(2*sqrt(dens) = |d/dy(g^-1(y))|
U=sqrt(U);

%% extract the expectation of the density
for k=1:size(dens,1)
    expect(k)=trapz(xi(k,:),xi(k,:).*dens(k,:));
end
conf=max([expect-U,expect+U],0);
