function [AICc,BIC,model,optimizedLogLikelihood]=modelSelection(X,Xi,Y,varargin)
% [AICc,BIC,model,optimizedLogLikelihood]=modelSelection(X,Xi,Y,varargin)
% function that uses densities for model selection
%
% Input:
%  - X: vector of windows; rescaled to logarithmic axis
%  - Xi: matrix; support axis for Y; rescaled to logarithmic axis
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
% See also griddedinterpolant, fminsearch, lsqnonlin
%
% Ton & Daffertshofer, Model selection for identifying power-law scaling
% Neuroimage 136:215-26, 2016, doi: 10.1016/j.neuroimage.2016.01.008
%
%                                              (c) marlow 2012-17
%                                                 (c) robert 2016
%                                     latest update March 7, 2017
%
% This file is released under the terms of the GNU General Public License,
% version 3. See http://www.gnu.org/licenses/gpl.html

%% optimization options
numberOfIterations=3;
options=optimset('TolX',1e-10,'Display','off');
if exist('optimoptions','file')
    lsqOptions=optimoptions('lsqnonlin','Display','off');
else
    lsqOptions=[];
end

%% set defaults and check variable input
model=defineModelCatalogue([],X);    % default model list
MM=find(strncmpi(varargin,'mod',3)); % adjust according to input parameters
if ~isempty(MM), MM=varargin{MM+1};
else
    MM=1:size(model,1);
end
model=model(MM,:);
verbose=sum(strncmpi(varargin,'ver',3))~=0;

optimizedLogLikelihood=zeros(size(model,1),1);
AICc=zeros(size(model,1),1);
BIC=zeros(size(model,1),1);

Y0=nan(size(X)); % set the intial value to the expectation value;
for k=1:numel(Y0)
    Y0(k)=trapz(Xi(k,:),Xi(k,:).*Y(k,:));
end

Ym=cell(size(X));
for k=1:numel(Ym) % get interpolation grid to link data to the model(s)
    % added nearest neighbor extrapolation to avoid NaNs
    Ym{k}=griddedInterpolant(Xi(k,:),Y(k,:),'pchip','nearest');
end

if verbose, fprintf('%s: optimizing models',mfilename); end

for k=1:size(model,1)

    if verbose, fprintf(' %d',k); end
        
    % initial guess on basis of least squares minimization of averaged values
    if ~isempty(lsqOptions)       
        flsq=@(x,y,c,mu) y-mu(x,c);
        initialValue=lsqnonlin(@(C) flsq(X,Y0,C,model{k,1}),model{k,2},[],[],lsqOptions);
    else
        flsq=@(x,y,c,mu) sum((y-mu(x,c)).^2);
        initialValue=fminsearch(@(C) flsq(X,Y0,C,model{k,1}),model{k,2},options);
    end
    
    B=zeros(numberOfIterations,numel(model{k,2}));
    L=nan(numberOfIterations,1);
    L(1)=logLikelihood(X,Ym,initialValue,model{k,1}); 
    
    for l=1:numberOfIterations % loop for perturbing initial conditions in fminsearch
        
        [B(l,:),L(l)]=fminsearch(@(C)logLikelihood(X,Ym,C,model{k,1}),initialValue,options);
        
        if sum(abs(initialValue-B(l,:)))>eps % select new initial condition
            initialValue=B(l,:);
        else
            initialValue=initialValue+B(l,:).*(rand(size(B(l,:)))-0.5);
        end
        
    end
    
    N=numel(X);
    K=numel(model{k,2});

    % find maximum LL given model{k,1} and corresponding model{k,2}
    [optimizedLogLikelihood(k),n]=max(-L);
    model{k,2}=B(n,:);
    
    % compute the values of the information criteria
    AICc(k)=-2*optimizedLogLikelihood(k)+2*K+(2*K*(K+1))/(N-K-1);
    BIC(k)=-2*optimizedLogLikelihood(k)+K*log(N);
    
end

if verbose, fprintf('\n'); end

%% computing the logLikelihood ... used when optimizing the models
function lp=logLikelihood(x,y,c,mu)
% lp=logLikelihood(x,y,c,mu)
% input: x/y: independent/dependent variable.
% c: parameter values for the model
% mu: function handle for the function to be fitted.
% calculates loglikelihood based on fit and densities of fluctuations
yfit=mu(x,c);
yl=zeros(size(y));
for k=1:numel(yl)    
    yl(k)=y{k}(yfit(k));
end
% lp=-nansum(log(yl)); %-log likelihood % use this for old Matlab versions
lp=-sum(log(yl),'omitnan'); %-log likelihood
% the omission of NaN should no longer be necessary given the added nearest
% neighbor extrapolation in griddedInterpolatant (see above) but it does
% not cost much extra time, i.e. we keep it for the time being
