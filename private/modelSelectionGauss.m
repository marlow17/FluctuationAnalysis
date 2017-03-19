function [AICc,BIC,model,optimizedLogLikelihood]=modelSelectionGauss(X,Y,varargin)
% [AICc,BIC,model,optimizedLogLikelihood]=modelSelectionGauss(X,Y,varargin)
% function that uses mean values for model selection by assuming that the
% distribution of fit-residuals per model is Gaussian.
%
% NOTE!!!
% This function has been added for completeness though we cannot recommend
% using it as in the majority of examples we tested, the assumption of
% Gaussianity for the destribution of residuals appeared violated.
%
% Input:
%  - X: vector of windows; rescaled to logarithmic axis
%  - Y: matrix; densities of F; rescaled to logarithmic axis
%  - 'models' followed by a list for subselection models (opt., def. 1:11)
%  - 'verbose' if set a progress report is given
%
% Output:
%  - AICc: cell; AIC values corrected for small sample sizes
%  - BIC: cell; BIC values
%
%    The models considered can be found in 'defineModelCatalogue.m'
%    The parameters for optimal LL models can be found in model{:,2};
%    in case of linear model, here model #1 c(2) corresponds to the
%    scaling exponent;
%    note: imaginary parts of parameters are being ignored!
%
%  - optimzedLogLikelihood: vector; opt. log-likelihood val. for all models
%
% See also modelSelection, fminsearch
%
%                                              (c) marlow 2012-16
%                                      latest update July 04, 2016
%
% This file is released under the terms of the GNU General Public License,
% version 3. See http://www.gnu.org/licenses/gpl.html

%% optimization options
numberOfIterations=3;
options=optimset('TolX',1e-10,'Display','off');

%% set defaults and check variable input
model=defineModelCatalogue([],X);
MM=find(strncmpi(varargin,'mod',3));
if ~isempty(MM), MM=varargin{MM+1};
else
    MM=1:size(model,1);
end
model=model(MM,:);
verbose=sum(strncmpi(varargin,'ver',3))~=0;

optimizedLogLikelihood=zeros(size(model,1),1);
AICc=zeros(size(model,1),1);
BIC=zeros(size(model,1),1);

if verbose, fprintf('%s: optimizing models ',mfilename); end

for k=1:size(model,1)

    if verbose, fprintf('.'); end
        
    % initial guess on basis of least squares minimization of averaged values
    flsq=@(x,y,c,mu) sum((y-mu(x,c)).^2);
    initialValue=fminsearch(@(C) flsq(X,Y,C,model{k,1}),model{k,2},options);

    S=nan(numberOfIterations,1);
    B=repmat(initialValue,numberOfIterations,1);
    
    for l=1:numberOfIterations % loop for perturbing initial conditions in fminsearch
        
        [B(l,:),S(l)]=fminsearch(@(C) flsq(X,Y,C,model{k,1}),initialValue,options);
        
        if sum(abs(initialValue-B(l,:)))>eps % select new initial condition
            initialValue=B(l,:);
        else
            initialValue=initialValue+B(l,:).*(rand(size(B(l,:)))-0.5);
        end
        
    end
    
    N=numel(X);
    K=numel(model{k,2});
    
    % find maximum LL given model{k,1} and corresponding model{k,2}
    [~,n]=min(S);
    model{k,2}=B(n,:);
    
    LL=S(n)/N;
    optimizedLogLikelihood(k)=-N/2*log(LL);
    
    % compute the values of the information criteria
    AICc(k)=-2*optimizedLogLikelihood(k)+2*K+(2*K*(K+1))/(N-K-1);
    BIC(k)=-2*optimizedLogLikelihood(k)+K*log(N);
    
end

if verbose, fprintf('\n'); end

