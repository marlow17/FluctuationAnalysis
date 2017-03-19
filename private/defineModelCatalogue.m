function model=defineModelCatalogue(n,X)
% model=setModelCatalogue(n) define a set of models to use for model
% comparison
%
% The models considered are
%
%    model={... % function handle                  % initial values for c
%        @(x,c) (c(1)+c(2)*x),                     [0,0]; ...
%        @(x,c) (c(1)+c(2)*x.^2),                  [0,0];...
%        @(x,c) (c(1)+c(2)*x+c(3)*x.^2),           [0,0,0]; ...
%        @(x,c) (c(1)+c(2)*x.^3),                  [0,0];...
%        @(x,c) (c(1)+c(2)*x+c(3)*x.^3),           [0,0,0];...
%        @(x,c) (c(1)+c(2)*x.^2+c(3)*x.^3),        [0,0,0];...
%        @(x,c) (c(1)+c(2)*x+c(3)*x.^2+c(4)*x.^3), [0,0,0,0];...
%        @(x,c) (c(1)+c(2)*exp(c(3)*x)),           [0,0,0]; ...
%        @(x,c) real((1/log(10))*log(c(1)*(1-exp(-c(2)*...
%                         exp(log(10)*x)))+c(3))), [0.1,0.1,0.1];...
%        @(x,c) ((x<c(4)).*(c(1)+c(2)*x))+...
%               ((x>=c(4)).*((c(1)+(c(2)-c(3))*c(4))+c(3)*x)), ...
%                                                  [0,0,0,mean(X)];...
%        @(x,c) ((x<c(5)).*(c(1)+c(2)*x))+...
%               ((x>=c(5)&x<c(6)).*((c(1)+(c(2)-c(3))*c(5))+c(3)*x))+...
%               ((x>=c(6)).*((c(1)+(c(2)-c(3))*c(4))+...
%                       (c(3)-c(4))*c(6)+c(4)*x)), [0,0,0,0,mean(X)*2/3,...
%                                                      mean(X)*4/3]
%      };
%
% the parameters for the models can be found in model{:,2};
% in case of linear model, here model #1 c(2) corresponds to the
% scaling exponent;
%
% Input:
%  - n: number(s) of models to be included (optional; default = 1:11)
%  - X: data to fit (helps to set initial values; optional; default = 1)
%
% OUTput:
%   model: defined models { function handle, [ parameters ] }
%
% See also modelSelection
%
% Ton & Daffertshofer, Model selection for identifying power-law scaling
% Neuroimage 136:215-26, 2016, doi: 10.1016/j.neuroimage.2016.01.008
%
%                                              (c) marlow 2012-16
%                                                 (c) robert 2016
%                                      latest update jun 18, 2016
%
% This file is released under the terms of the GNU General Public License,
% version 3. See http://www.gnu.org/licenses/gpl.html

if nargin<2 || isempty(X)
    X=1;
end

model={...
    @(x,c) (c(1)+c(2)*x),                     [0,0]; ...
    @(x,c) (c(1)+c(2)*x.^2),                  [0,0];...
    @(x,c) (c(1)+c(2)*x+c(3)*x.^2),           [0,0,0]; ...
    @(x,c) (c(1)+c(2)*x.^3),                  [0,0];...
    @(x,c) (c(1)+c(2)*x+c(3)*x.^3),           [0,0,0];...
    @(x,c) (c(1)+c(2)*x.^2+c(3)*x.^3),        [0,0,0];...
    @(x,c) (c(1)+c(2)*x+c(3)*x.^2+c(4)*x.^3), [0,0,0,0];...
    @(x,c) (c(1)+c(2)*exp(c(3)*x)),           [0,0,0]; ...
    @(x,c) real((1/log(10))*log(c(1)*(1-exp(-c(2)*exp(log(10)*x)))+c(3))), [0.1,0.1,0.1];...
    @(x,c) (x<c(4)).*(c(1)+c(2)*x)+(x>=c(4)).*(c(1)+(c(2)-c(3))*c(4)+c(3)*x), [0,0,0,mean(X)];...
    @(x,c) (x<c(5)).*(c(1)+c(2)*x)+(x>=c(5)&x<c(6)).*(c(1)+(c(2)-c(3))*c(5)+c(3)*x)+(x>=c(6)).*(c(1)+(c(2)-c(3))*c(4)+(c(3)-c(4))*c(6)+c(4)*x), [0,0,0,0,mean(X),mean(X)*2/3,mean(X)*4/3];...
    };

if nargin && ~isempty(n)
    model=model(n,:);
end