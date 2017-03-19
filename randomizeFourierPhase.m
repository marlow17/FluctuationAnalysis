function y=randomizeFourierPhase(x)
% y=randomizeFourierPhase(x) randomizes the Fourier phases of x
% cf Kantz & Schreiber, Nonlinear time series analysis, Cambridge Univ.
% Press, 1997, p98
%
% See also fft, ifft
%                                                     (c) marlow 2016
%
% This file is released under the terms of the GNU General Public License,
% version 3. See http://www.gnu.org/licenses/gpl.html

y=nan(size(x));
for k=1:size(x,2)
    
    phi=2*pi*rand(size(x,1),1);
    if rem(size(x,1),2)
        phi((end-1)/2+2:end,:)=-phi((end-1)/2+1:-1:2,:);
        phi(1,:)=0;
    else
        phi(end/2+2:end,:)=-phi(end/2:-1:2,:);
        phi(1,:)=0; phi(end/2+1,:)=0;
    end
    
    y(:,k)=fft(x(:,k));
    y(:,k)=ifft(y(:,k).*exp(sqrt(-1)*phi),size(y,1),'symmetric');
    
end
%
