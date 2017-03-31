function [alpha,details]=fluctuationAnalysis(fBm,n,nf,varargin)
%[alpha,details]=fluctuationAnalysis(fBm,n,nf,...)
%
% Main function for the fluctuation analysis of (a set of) fBm series
%
% Reference:
% Ton & Daffertshofer, Model selection for identifying power-law scaling
% Neuroimage 136:215-26, 2016, doi:10.1016/j.neuroimage.2016.01.008
%
% Example usages:
%
% >> x=randn(10,1000)'; % generate Gaussian white noise (10 trials with
%                       % 1000 samples each
%
%                       % default usage...
% >> alpha=fluctuationAnalysis(cumsum(x));
%                       % or if interested in the progress...
% >> alpha=fluctuationAnalysis(cumsum(x),[],[],'verbose');
%                       % or to visualize results...
% >> alpha=fluctuationAnalysis(cumsum(x),[],[],'graphics,'verbose');
%                       % this equals
% >> alpha=fluctuationAnalysis(cumsum(x),[],[],'graphics','verbose','DFA+');
%                       % all calls above include a non-parametric model
%                       % selection based on BIC (or AICc) as described in
%                       % the aforementioned paper
%
%                       % you can also use conventional DFA
% >> alpha=fluctuationAnalysis(cumsum(x),[],[],'graphics','verbose','DFA');
%
%                       % or conventional SDA
% >> alpha=fluctuationAnalysis(cumsum(x),[],[],'graphics','verbose','SDA');
%
%                       % or if interested in DFA with parametric model
%                       % selection (does not work in most examples tested
%                       % so far because the residuals of the model fit are
%                       % usually not normally distributed
% >> alpha=fluctuationAnalysis(cumsum(x),[],[],'graphics','verbose','DFA-');
%
%
% Prelude:
%
% If studying fractional Gaussian noise (fGn) one first has to compute
% the cumulative sum of the corresponding time series before estimating the
% exponent; recall diff(fBm) ~ fGn => H(fBm)=H(fGn)+1
%
%
% Input:
% - fBm       NxM matrix of M series with N samples each
% - n         segment sized to include estimating fluctuation strenghs
%             if numel(n)==3 
%                  n(1) = Nmin (default = min(10,length(fGn)))
%                  n(2) = Nmax (default = length(fGn)/4)
%                  n(3) = N = # of n-s (default = length(fGn)/4)
%             end
% - nf        indices to include in the linear fit default = n (see above)
%             if nf differs from n then densities will be
%             recomputed which will cost extra computation time
% - 'fs'      followed by a value: sampling rate; by using this option all
%             input/output value will be handled/given in [seconds] rather
%             than [samples] as dependent variable (default is samples)
%             detrending (optional; def.=1 for DFA, - for SDA)
% - 'order'   followed by an integer: set order of polynomial used for
%             detrending (optional; def.=1 for DFA, - for SDA)
% - 'bins'    followed by an integer: set the number of bins when
%             estimating the probability densities of F
% - 'DFA+'    flag to run DFA + non-parametric model comparison (default)
% - 'DFA-'    flag to run DFA + parametric model comp. (not recommended!)
% - 'DFA'     flag to run a traditional DFA instead of DFA+model comparison
% - 'SDA'     flag to run traditional SDA instead of DFA+model comparison
% - 'graphic' flag to add graphical output
%             (default = true if nargout==0, otherwise false = no graphics)
% - 'verbose' flag to add a text report
%             (default = true if nargout==0, otherwise false = no report)
%
% Output:
% - alpha     alpha = Hurst exponent if input is indeed an fGn and order=1
%             (will be set to NaN if model selection yield a min(BIC) for a
%             model other than the linear one (see ModelSelection)
%             Note that BIC (and AICc) strongly favor models with a small
%             numbe of parameters; you may also check the value of the
%             log-likelihood (i.e. details.loglikelihood, see below)
% - details   structure containing more details that depend on the the type
%             of analysis
%
% for the DFA + model comparison the details struct contains:
%
%   .F        divergence (diffusion) curve based on expectation values
%             Note that this may differ from conventional DFA as these 
%             values represent median values of the distributions of
%             fluctuation strenghts of consecutive segments rather than
%             'simple' mean values
%   .n        indices used in the estimate
%   .dens     kernel density estimates of F
%   .xi       logarithmically spaced support axis
%   .model    used models { function handle, [ parameters ] }
%   .BIC      BIC values per model
%   .AICc     AICc values per model
%   .loglikelihod ... per model
%
%    Comment: Support axis xi is logarithmically spaced to improve
%             accuracy. Both dens and xi correspond to the root mean
%             squared fluctuations, hence for modelSelection in double
%             logarithmic coordinates Fdens first has to be transformed to
%             double logarithmic coordinates according to the usual change
%             of variables
%
% for conventional DFA / SDA the details struct contains:
%
%   .F        divergence (diffusion) curve based on mean values
%   .n        indices used in the estimate
%   .R2       mean R^2 (coefficient of determination)
%   .conf     95% confidence interval of F
%   .P        coefficients of linear fit of the psd P=[beta,p_0]
%   .S        S contains fields for the triangular factor (R) from a QR
%             decomposition of the Vandermonde matrix of X, the degrees of
%             freedom (df), and the norm of the residuals (normr).
%
% See also detrendedFluctuationAnalysis, diffusionAnalysis, modelSelection,
% modelSelectionGauss (not recommended), reportDetails
%
%                                              (c) marlow 2012-17
%                                     latest update March 5, 2017
%
% This file is released under the terms of the GNU General Public License,
% version 3. See http://www.gnu.org/licenses/gpl.html

%% set default values
order=[];
% range over which fluctuations should be estimated
% default = 10 ... (signal length)/4 in (signal length)/4 steps
if nargin<2 || isempty(n)
    n=[min(10,size(fBm,1)),size(fBm,1)/4,round(size(fBm,1)/8)];
end
% range to fit the straight line (if appropriate) = range n (see above)
if nargin<3 || isempty(nf), nf=n; end

supportedAlgorithm={ 'DFA+','DFA-','DFA','SDA' };
algorithm=supportedAlgorithm{1};
graphics=false;
verbose=false;
fs=1;

%% check variable input (flags)
if numel(varargin)
    j=false(size(varargin));
    for i=1:numel(varargin)
        j(i)=ischar(varargin{i});
    end
    a=intersect(varargin(j),supportedAlgorithm);
    if ~isempty(a), algorithm=a{1}; end
    oi=find(strncmpi(varargin,'ord',3));
    if ~isempty(oi), order=varargin{oi+1}; end
    oi=find(strncmpi(varargin,'bin',3));
    if ~isempty(oi), M=varargin{oi+1}; end
    oi=find(strncmpi(varargin,'fs',2));
    if ~isempty(oi), fs=varargin{oi+1}; end
    graphics=sum(strncmpi(varargin,'gra',3))~=0;
    verbose=sum(strncmpi(varargin,'ver',3))~=0;
end
if verbose, verboseString='verbose'; else, verboseString=''; end

if ~exist('order','var') || isempty(order)
    order=1-strcmp(algorithm,'SDA');
end
if ~exist('M','var') || isempty(M)
    M=100;
end

% define indices over which the fluctuation should be estimated, if not
% provided explictly
n=range2index(algorithm,[max(order+1,round(n(1:2)*fs)),n(3)],order+2,size(fBm,1));
% if range of linear fits is a cell array, e.g. if there is more than one
% range of interest, then output will also be a cell array, otherwise not
if iscell(nf)
    notcell=false;
else
    notcell=true;
    nf={ nf };
end

for j=1:numel(nf)
    
    nf{j}=range2index(algorithm,[max(order+1,round(nf{j}(1:2)*fs)),...
        nf{j}(3)],order+2,size(fBm,1));
    
end

alpha=cell(size(nf));
details=cell(size(nf));

%% remove the DC (needed only when estimating the diffusion but it doesn't harm otherwise)
fBm=detrend(fBm,'constant');

%% run the analysis
if verbose, fprintf('\n%s: running %s ...\n',mfilename,algorithm); tic; end

switch algorithm
    
    case 'DFA+' % DFA via non-parametric model selection
        
        % first: get the fluctuation densities over then entire range
        [dens,xi,F,conf]=detrendedDensities(fBm,n,M,'order',order,verboseString);
        
        % second: get optimal model for every fit range
        for j=1:numel(nf)
            
            if verbose, fprintf('%s: range [%g ... %g] ...\n',...
                    mfilename,nf{j}(1)/fs,nf{j}(end)/fs); end

            details{j}.algorithm=algorithm;

            % store the full range results
            details{j}.fs=fs;
            details{j}.n_base=n/fs;
            details{j}.dens_base=dens;
            details{j}.xi_base=xi;
            details{j}.F_base=F;
            details{j}.conf_base=conf;

            if numel(n)==numel(nf{j}) && sum(n(:)-nf{j}(:))==0 % case interval2fit == interval2evaluate
                details{j}.dens=details{j}.dens_base;
                details{j}.xi=details{j}.xi_base;
                details{j}.F=details{j}.F_base;
            else
                [details{j}.dens,details{j}.xi,details{j}.F]=...
                    detrendedDensities(fBm,nf{j},M,'order',order,verboseString);
            end
            details{j}.xi=details{j}.xi;
            details{j}.n=nf{j}/fs;
            
            % model selection ... see doc modelSelection
            [...
                details{j}.AICc,...
                details{j}.BIC,...
                details{j}.model,...
                details{j}.loglikelihood]=...
                   modelSelection(...
                            log10(details{j}.n),...
                            log10(details{j}.xi),...
                            log(10)*details{j}.xi.*details{j}.dens,...
                            verboseString...
                            );
            
            [~,i]=min(details{j}.BIC);
            if i==1
                alpha{j}=details{j}.model{i,2}(2);
            else
                alpha{j}=nan;
            end
            
        end

    case 'DFA-' % DFA + model selection assuming Gaussianity
        
        % first: get the mean fluctuation values over then entire range
        [F,conf]=detrendedFluctuationAnalysis(fBm,n,'order',order,verboseString);
        
        % second: model selection for every fit range
        for j=1:numel(nf) % loop over to be fitted ranges

            if verbose, fprintf('%s: range [%g ... %g] ...\n',mfilename,nf{j}(1)/fs,nf{j}(end)/fs); end

            details{j}.algorithm=algorithm;

            % store the full range results
            details{j}.fs=fs;
            details{j}.n_base=n/fs;
            details{j}.F_base=F;
            details{j}.conf_base=conf;

            % rescaling to logarithmic axes
            details{j}.n=nf{j}/fs;
            details{j}.F=interp1(n/fs,F,nf{j}/fs,'pchip');

            % model selection under Gaussian assumption
            [...
                details{j}.AICc,...
                details{j}.BIC,...
                details{j}.model,...
                details{j}.loglikelihood]=...
                   modelSelectionGauss(...
                            log10(details{j}.n),...
                            log10(details{j}.F),...
                            verboseString...
                            );
            
            [~,i]=min(details{j}.BIC);
            if i==1
                alpha{j}=details{j}.model{i,2}(2);
            else
                alpha{j}=nan;
            end
            
        end
        
    case 'DFA' % conventional DFA
        
        % first: get the mean fluctuation values over then entire range
        [F,conf]=detrendedFluctuationAnalysis(fBm,n,'order',order,verboseString);

        % second: regression for every fit range
        for j=1:numel(nf) % loop over to be fitted ranges

            if verbose, fprintf('%s: range [%g ... %g] ...\n',mfilename,nf{j}(1)/fs,nf{j}(end)/fs); end

            details{j}.algorithm=algorithm;

            % store the full range results
            details{j}.fs=fs;
            details{j}.n_base=n/fs;
            details{j}.F_base=F;
            details{j}.conf_base=conf;

            % rescaling to logarithmic axes
            details{j}.n=nf{j}/fs;
            details{j}.F=interp1(n/fs,F,nf{j}/fs,'pchip');

            % fit the straight line
            [details{j}.P,details{j}.S]=polyfit(log10(details{j}.n),log10(details{j}.F),1);

            % get the slope and estimate R-squared
            alpha{j}=details{j}.P(1);
            details{j}.R2=(1-details{j}.S.normr^2/var(log10(details{j}.F))/...
                (numel(details{j}.n)-numel(details{j}.P)-1));
            
        end
        
    case 'SDA' % SDA
        
        % first: get the mean diffusion values over then entire range
        [F,conf]=diffusionAnalysis(fBm,n,verboseString);
        
        % second: regression for every fit range
        for j=1:numel(nf) % loop over to be fitted ranges
            
            if verbose, fprintf('%s: range [%g ... %g] ...\n',mfilename,nf{j}(1)/fs,nf{j}(end)/fs); end

            details{j}.algorithm=algorithm;

            % store the full range results
            details{j}.fs=fs;
            details{j}.n_base=n/fs;
            details{j}.F_base=F;
            details{j}.conf_base=conf;

            details{j}.n=nf{j}/fs;
            details{j}.F=interp1(n/fs,F,nf{j}/fs,'pchip');
            
            % fit the straight line
            [details{j}.P,details{j}.S]=polyfit(details{j}.n,details{j}.F,1);
            
            % get the slope and estimate R-squared
            alpha{j}=details{j}.P(1);
            details{j}.R2=(1-details{j}.S.normr^2/var(details{j}.F)/...
                (numel(details{j}.n)-numel(details{j}.P)-1));
            
        end
        
end

if verbose, fprintf('%s: ',mfilename); toc; end

%% verbose mode: create some output (plot and/or report)
if fs~=1, strg='\Delta{t}'; else, strg='\Delta{n}'; end
reportDetails(details,verbose,graphics,strg);

%% set output in as row vector of cells ...
if nargout && notcell
    alpha=alpha{:};
    details=details{:};
end


%% define index array over a range
function n=range2index(algorithm,n,nmin,nmax)
%
if n(1)>=n(2)
    error('range/index size %d ... %d is bogus',n(1),n(2));
end

if numel(n)==3
    switch algorithm
        
        case { 'DFA','DFA+','DFA-' }
            % create indices that are equidistant on a log scale
            % this allows for an unbiased estimate of the slope in log-log
            n=logspace(log10(n(1)),log10(n(2)),n(3))';
            
        case 'SDA'
            % for SDA a linear index range will do...
            n=linspace(n(1),n(2),n(3))';
            
    end
end

n=n(n>=nmin & n<=nmax); % we need to avoid NaNs and Infs
n=unique(round(n(:)));  % after rounding some indices might be double



