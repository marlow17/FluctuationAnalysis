function reportDetails(details,report,graphics,strg)
% reportDetails(details,report,graphics,strg) helper function for reporting
% details of DFA+model selection, conventional DFA, or SDA
%
% input
% - details = struct output of FluctuationAnalysis
% - report = boolean, if true a table is printed in the command window
% - graphics = boolean, if true a plot is generated
%
% see also FluctuationAnalysis
%
%                                                     (c) marlow 2016
%
% This file is released under the terms of the GNU General Public License,
% version 3. See http://www.gnu.org/licenses/gpl.html

% set some defaults
if nargin<4, strg='\Delta{n}'; end
if nargin<3, graphics=true; end
if nargin<2, report=true; graphics=true; end

if usejava('jvm') && ~feature('ShowFigureWindows')
    graphics=false;
end

if ~iscell(details)
    details={ details };
end

if report
    printDetails(details);
end

if graphics 
    plotDetails(details,strg);
end

%%
function printDetails(details)

tablesize=60;

fprintf('\ntotal range of assessment: epoch sizes [%g ... %g] \n',...
    details{1}.n_base(1),details{1}.n_base(end));

switch details{1}.algorithm
    
    case 'DFA' 
        
        fprintf('%s\n',repmat('-',tablesize,1));
        
        for j=1:numel(details)
            
            fprintf('fit results over epoch sizes [%g ... %g]\n',...
                details{j}.n(1),details{j}.n(end))
            fprintf('-> model: %s',fcn2str({'c(2)+c(1)*x',details{j}.P}));
            fprintf(' (R2=%f)\n',details{j}.R2);
            fprintf('%s\n',repmat('-',tablesize,1));

        end

    case 'SDA' 

        fprintf('%s\n',repmat('-',tablesize,1));
        
        for j=1:numel(details)
            
            fprintf('fit results over epoch sizes [%g ... %g]\n',...
                details{j}.n(1),details{j}.n(end))
            fprintf('-> model: %s',fcn2str({'c(2)+c(1)*x',details{j}.P}));
            fprintf(' (R2=%f)\n',details{j}.R2);
            fprintf('%s\n',repmat('-',tablesize,1));

        end
        
    case { 'DFA+', 'DFA-' }
        
        fprintf('number of included models:\t%d\n',size(details{1}.model,1));
        fprintf('%s\n',repmat('-',tablesize,1));
        
        for j=1:numel(details)
            
            fprintf('fit results over epoch sizes [%g ... %g]\n',...
                details{j}.n(1),details{j}.n(end))
            fprintf('   | log-lh |   AICc |    BIC | model\n');
            fprintf('%s\n',repmat('-',tablesize,1));
            
            for i=1:size(details{j}.model,1)
                fprintf('%2d | %6d | %6d | %6d | %s\n',i, ...
                    round(details{j}.loglikelihood(i)), ...
                    round(details{j}.AICc(i)), ...
                    round(details{j}.BIC(i)), fcn2str(details{j}.model(i,1)) ...
                    );
            end
            fprintf('%s\n',repmat('-',tablesize,1));
            [~,k]=max(details{j}.loglikelihood);
            fprintf('max. log-likelihood\t-> model %2d: %s\n',...
                k,fcn2str(details{j}.model(k,:)));
            [~,k]=min(details{j}.AICc);
            fprintf('minimum AICc\t\t-> model %2d: %s\n',...
                k,fcn2str(details{j}.model(k,:)));
            [~,k]=min(details{j}.BIC);
            fprintf('minimum BIC\t\t-> model %2d: %s\n',...
                k,fcn2str(details{j}.model(k,:)));
            fprintf('%s\n',repmat('-',tablesize,1));
            
        end
end

%%
function plotDetails(details,xstrg)

h=zeros(size(details));
strg=cell(size(details));

% plot confidence interval / density function over the entire range
switch details{1}.algorithm

    case 'DFA'
        
        % plot mean divergence plus confidence interval
        cla; patch([details{1}.n_base;details{1}.n_base(end:-1:1)],...
            [details{1}.conf_base(:,1);details{1}.conf_base(end:-1:1,2)],...
            ones(1,3)*.7,'edgecolor','none','facealpha',0.5); hold on;
        hold on;
        plot(details{1}.n_base,details{1}.F_base,'k-','linewidth',2);
        
        co=get(gca,'colororder');

        for j=1:numel(details)

            [pp,delta]=polyval(details{j}.P,log10(details{j}.n),details{j}.S);
            delta=[pp-delta,pp+delta];
            delta=10.^delta;
            pp=10.^pp;
            patch([details{j}.n;details{j}.n(end:-1:1)],[delta(:,1);delta(end:-1:1,2)],...
                co(j,:)*.7,'edgecolor','none','facealpha',0.5);
            h(j)=plot(details{j}.n,pp,'color',co(j,:),'linewidth',2);
            strg{j}=fcn2str({'c(2)+c(1)*x',details{j}.P});

        end
        
        legend(h,strg,'location','southeast');
        set(gca,'xscale','log','yscale','log');
        ylabel(sprintf('F(%s)',xstrg));
        xmin=floor(log10(min(details{1}.n_base(details{1}.n_base>0))));
        xmax=ceil(log10(max(details{1}.n_base(details{1}.n_base>0))));
        set(gca,'xlim',10.^[xmin,xmax],'xtick',10.^(xmin:1+((xmax-xmin)>5):xmax));

    case 'SDA'
        
        % plot mean diffusion plus confidence interval
        cla; patch([details{1}.n_base;details{1}.n_base(end:-1:1)],...
            [details{1}.conf_base(:,1);details{1}.conf_base(end:-1:1,2)],...
            ones(1,3)*.7,'edgecolor','none','facealpha',0.5); hold on;
        hold on;
        plot(details{1}.n_base,details{1}.F_base,'k-','linewidth',2);
        
        co=get(gca,'colororder');

        for j=1:numel(details)
            
            [pp,delta]=polyval(details{j}.P,details{j}.n,details{j}.S);
            delta=[pp-delta,pp+delta];
            patch([details{j}.n;details{j}.n(end:-1:1)],[delta(:,1);delta(end:-1:1,2)],...
                co(j,:)*.7,'edgecolor','none','facealpha',0.5);
            h(j)=plot(details{j}.n,pp,'color',co(j,:),'linewidth',2);
            strg{j}=fcn2str({'c(2)+c(1)*x',details{j}.P});

        end
        
        legend(h,strg,'location','southeast');
        set(gca,'xscale','lin','yscale','lin');
        ylabel(sprintf('\\Delta{X}^2(%s)',xstrg));

    case 'DFA+'
            
        % plot mean divergence plus density functions
        dens_rs=log(10)*details{1}.xi_base.*details{1}.dens_base;
        xi_rs=log10(details{1}.xi_base);
        Faxis=linspace(min(xi_rs(:)),max(xi_rs(:)),2000);
        d=nan(numel(Faxis),numel(details{1}.n_base));
        for i=1:numel(details{1}.n_base)
            d(:,i)=interp1(xi_rs(i,:),dens_rs(i,:),Faxis,'pchip',nan)';
        end
        cla; colormap(1-gray(256));
        surf(details{1}.n_base,10.^Faxis,d.^(1/3),...
            'edgecolor','none','facecolor','interp','facealpha',0.5);
        view(0,90); hold on;
        plot3(details{1}.n_base,details{1}.F_base,...
            (max(d(:)).^(1/3))*ones(size(details{1}.F_base)),...
            'k-','linewidth',2);
        co=get(gca,'colororder');
        
        for j=1:numel(details)
            
            [~,bestmodel]=min(details{j}.BIC);
            pp=details{j}.model{bestmodel,1}(log10(details{j}.n),details{j}.model{bestmodel,2});
            h(j)=plot3(details{j}.n,10.^pp,(max(d(:)).^(1/3))*ones(size(pp)),...
                'color',co(j,:),'linewidth',2);
            strg{j}=fcn2str(details{j}.model(bestmodel,:));

        end
                
        legend(h,strg,'location','southeast');
        set(gca,'xscale','log','yscale','log');
        ylabel(sprintf('F(%s)',xstrg));
        xmin=floor(log10(min(details{1}.n_base(details{1}.n_base>0))));
        xmax=ceil(log10(max(details{1}.n_base(details{1}.n_base>0))));
        set(gca,'xlim',10.^[xmin,xmax],'xtick',10.^(xmin:1+((xmax-xmin)>5):xmax));
        
    case 'DFA-'
            
        % plot mean divergence 
        cla; plot(details{1}.n_base,details{1}.F_base,'k-','linewidth',2);
        hold on;
        co=get(gca,'colororder');
        
        for j=1:numel(details)
            
            ni=details{j}.n;
            [~,bestmodel]=min(details{j}.BIC);
            pp=details{j}.model{bestmodel,1}(log10(details{j}.n),details{j}.model{bestmodel,2});
            h(j)=plot(ni,10.^pp,...
                'color',co(j,:),'linewidth',2);
            strg{j}=fcn2str(details{j}.model(bestmodel,:));

        end
            
        legend(h,strg,'location','southeast');
        set(gca,'xscale','log','yscale','log');
        ylabel(sprintf('F(%s)',xstrg));
        xmin=floor(log10(min(details{1}.n_base(details{1}.n_base>0))));
        xmax=ceil(log10(max(details{1}.n_base(details{1}.n_base>0))));
        set(gca,'xlim',10.^[xmin,xmax],'xtick',10.^(xmin:1+((xmax-xmin)>5):xmax));

end

hold off;
title(details{1}.algorithm);
set(gca,'ylim',[min(details{1}.F_base),max(details{1}.F_base)],...
    'xgrid','on','ygrid','on');
if strcmp(xstrg,'\Delta{t}')
    xlabel([xstrg ' [s]']);
else
    xlabel(xstrg);
end

function strg=fcn2str(model,varname)

if ischar(model{1})
    strg=model{1};
else
    strg=func2str(model{1});
end
i=strfind(strg,'@(x,c)');
if ~isempty(i)
     strg=strg(i+6:end);
end


if numel(model)>1

    specials={ ... % this is a quite ugly solution... to be updated
        'c(1)+(c(2)-c(3))*c(4)',...
        'model{2}(1)+(model{2}(2)-model{2}(3))*model{2}(4)';...
        'c(1)+(c(2)-c(3))*c(5)',...
        'model{2}(1)+(model{2}(2)-model{2}(3))*model{2}(5)'; ...
        'c(1)+(c(2)-c(3))*c(4)+(c(3)-c(4))*c(6)',...
        'model{2}(1)+(model{2}(2)-model{2}(3))*model{2}(4)+(model{2}(3)-model{2}(4))*model{2}(6)'; ...
        };
    for j=1:size(specials,1)
            
        i=strfind(strg,specials{j,1});
        if isempty(i), continue; end
        s=eval(specials{j,2});
        strg=[strg(1:i-1) sprintf('%.2f',s) strg(i+numel(specials{j,1}):end)];
        
    end
    
    for k=1:numel(model{2})
                    
        se=sprintf('c(%d)',k);
        while ~isempty(strfind(strg,se))
            i=strfind(strg,se);
            strg=[strg(1:i-1) ...
                sprintf('%.2f',model{2}(k)) strg(i+numel(se):end)];
            i=strfind(strg,'+-');
            if ~isempty(i)
                strg=strg([1:i-1,i+1:end]);
            end
        end
        
    end

end

s2remove={'.*','*','./','.^'};

for k=1:numel(s2remove)
    while ~isempty(strfind(strg,s2remove{k}))
        i=strfind(strg,s2remove{k});
        strg=strg([1:i-1,i+1:end]);
    end
end

if nargin>1 && ~isempty(varname)
    while ~isempty(strfind(strg,'x'))
        i=strfind(strg,'x');
        strg=[strg(1:i-1) '\Delta{n}' strg(i+1:end)];
    end
end


