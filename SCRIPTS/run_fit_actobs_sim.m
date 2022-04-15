% Rouault et al. (2022) BioRxiv. Controllability boosts neural and
% cognitive correlates of changes-of-mind in uncertain environments.
% This script reproduced the main behavioral and psychometric analyses and
% statistics reported in the paper.

% task 1 = uncontrollable C- condition
% task 2 = controllable C+ condition

% This is a two-part script:
% 1) Run the first section to fit the model
% 2) Run the second section to analyse fitting output

% You first need to install Bayesian Adaptive Direct Search:
% https://github.com/lacerbi/bads


%% 1) Run the first section to fit the model

% clear workspace
clear;
close all;
clc;

% set analysis parameters
expename = 'Experiment_1';

% add fitting toolboxes to path
addpath('./bads/'); % Bayesian Adaptive Direct Search

% define list of subjects
switch expename
    case 'Experiment_1'
        subjlist = setdiff(01:17,[04,16]);
    case 'Experiment_2A'
        subjlist = setdiff(01:20,[06,15]);
    otherwise
        error('Undefined experiment!');
end
nsubj = numel(subjlist);

for isubj = 1:nsubj
    
    for itask = 1:2
        
        fprintf('fitting sub%02d task%d\n\n\n',isubj,itask);
        
        % load experiment data
        switch expename
            case 'Experiment_1'
                fname = sprintf('./DATA/%s_S%02d_task%d_expdata.mat','ACTOBS_C',subjlist(isubj),itask);
            case 'Experiment_2A'
                fname = sprintf('./DATA/%s_rule1_S%02d_task%d_expdata.mat','ACTOBS_D_',subjlist(isubj),itask);
        end
        load(fname,'dat');
        
        % recreate additional information
        blkind = 0; % current block index
        seqdir = 0; % current "true" sequence direction (1 or 2)
        epinum = 0; % episode number in current block, reset at break
        seqpos = 0; % sequence position in current episode, reset at reversal
        dat.epinum = nan(size(dat.seqdir));
        dat.seqpos = nan(size(dat.seqdir));
        for iseq = 1:numel(dat.seqdir)
            if dat.blkind(iseq) ~= blkind
                blkind = dat.blkind(iseq);
                seqdir = dat.seqdir(iseq);
                epinum = 1;
                seqpos = 1;
            elseif dat.seqdir(iseq) ~= seqdir
                seqdir = dat.seqdir(iseq);
                epinum = epinum+1;
                seqpos = 1;
            end
            dat.epinum(iseq) = epinum;
            dat.seqpos(iseq) = seqpos;
            seqpos = seqpos+1;
        end
        dat.seqlen = cellfun(@length,dat.smpang);
        dat.seqllr = cellfun(@sum,dat.smpllr);
        
        % store two times participant data: "c" here and "csub" later on
        c = getc(dat.seqind',dat.seqpos',dat.seqdir',dat.seqllr',dat.rbef',dat.raft',dat.cbef',dat.caft');
        out_fit.c = c ;
        
        % configure model fitting
        cfg         = [];
        cfg.seqind  = dat.seqind; % sequence index in current block
        cfg.seqpos  = dat.seqpos; % sequence position in current episode
        cfg.seqdir  = dat.seqdir; % sequence direction
        cfg.seqllr  = cellfun(@sum,dat.smpllr); % sequence evidence
        cfg.seqlen  = cellfun(@length,dat.smpllr); % sequence length
        cfg.rbef    = dat.rbef; % response before sequence
        cfg.raft    = dat.raft; % response after sequence
        cfg.cbef    = dat.cbef; % confidence before sequence
        cfg.caft    = dat.caft; % confidence after sequence
        cfg.nsmp    = 1e3; % number of samples used by particle filter
        cfg.nres    = 1e3; % number of bootstrap resamples (vbmc)
        cfg.nval    = 1e2; % number of validation samples (bads)
        cfg.nrun    = 5; % number of random starting points (bads)
        cfg.verbose = 2; % display level
        cfg.ssel    = 0; % no selection noise
        
        cfg.fitrev  = false;
        cfg.fitrep  = true;
        cfg.fitcnf  = true;
        cfg.resamp  = false;
        
        cfg.gcnf = []; % "1" to force the confidence gain on switches to be nominal
        cfg.scnf = []; % "1" to force the confidence noise to be zero
        
        % fit model
        out_fit = fit_actobs_sim(cfg);
        
        % save results
        switch expename
            case 'Experiment_1'
                fname_out = sprintf('%s_sub%02d_task%d_fit_bads.mat',expename,subjlist(isubj),itask);
            case 'Experiment_2A'
                fname_out = sprintf('%s_rule1_sub%02d_task%d_fit_bads.mat',expename,subjlist(isubj),itask);
        end
        save(fname_out,'out_fit');
        
    end
    
end

out_fit




%% 2) Run the second section to analyse fitting output

expename = 'Experiment_1';
% define list of subjects
switch expename
    case 'Experiment_1'
        subjlist = setdiff(01:17,[04,16]);
    case 'Experiment_2A'
        subjlist = setdiff(01:20,[06,15]);
    otherwise
        error('Undefined experiment!');
end
nsubj = numel(subjlist);

all_rrev = NaN(nsubj,2,8);
all_crev = NaN(nsubj,2,8);
all_rrev_fit = NaN(nsubj,2,8);
all_crev_fit = NaN(nsubj,2,8);

pconf = NaN(nsubj,2);
pconf_fit = NaN(nsubj,2);

all_rrep = NaN(nsubj,2,8);
all_crep = NaN(nsubj,2,8);
all_rrep_fit = NaN(nsubj,2,8);
all_crep_fit = NaN(nsubj,2,8);

for isubj = 1:nsubj
    for itask = 1:2
        
        switch expename
            case 'Experiment_1'
                fname_out = sprintf('%s_sub%02d_task%d_fit_bads.mat',expename,subjlist(isubj),itask);
            case 'Experiment_2A'
                fname_out = sprintf('%s_rule1_sub%02d_task%d_fit_bads.mat',expename,subjlist(isubj),itask);
        end
        load(fname_out)
        
        % retrieve parameters
        h(isubj,itask)    = out_fit.h;
        sinf(isubj,itask) = out_fit.sinf;
        scnf(isubj,itask) = out_fit.scnf;
        tcnf(isubj,itask) = out_fit.tcnf;
        gcnf(isubj,itask) = out_fit.gcnf;
        
        % retrieve reversal curves
        all_rrev(isubj,itask,:) = out_fit.csub.rrev;
        all_crev(isubj,itask,:) = out_fit.csub.crev;
        
        all_rrev_fit(isubj,itask,:) = mean(out_fit.cfit.rrev,2);
        all_crev_fit(isubj,itask,:) = mean(out_fit.cfit.crev,2);
        
        % retrieve repetition curves
        all_rrep(isubj,itask,:) = out_fit.csub.rrep;
        all_crep(isubj,itask,:) = out_fit.csub.crep;
        
        all_rrep_fit(isubj,itask,:) = mean(out_fit.cfit.rrep,2);
        all_crep_fit(isubj,itask,:) = mean(out_fit.cfit.crep,2);
        
        % load experiment data to get pconf for error bars subplot4
        switch expename
            case 'Experiment_1'
                fname = sprintf('./DATA/%s_S%02d_task%d_expdata.mat','ACTOBS_C',subjlist(isubj),itask);
            case 'Experiment_2A'
                fname = sprintf('./DATA/%s_rule1_S%02d_task%d_expdata.mat','ACTOBS_D_',subjlist(isubj),itask);
        end
        load(fname,'dat');
        
        pconf(isubj,itask) = mean(dat.caft == 2);
        pconf_fit = pconf;
    end
end


% colors
rgb = [[0,0.5,1];[1,0,0.5]];
rgb = 0.75*rgb+0.25;


% Plot response and confidence reversal curves
figure(1)

% 1/ for responses
% fit
yavghat = squeeze(mean(all_rrev_fit,1)); % mean
yerrhat = squeeze(std(all_rrev_fit,[],1)/sqrt(nsubj)); % S.E.M.
% data
yavg = squeeze(mean(all_rrev,1));
yerr = squeeze(std(all_rrev,[],1)/sqrt(nsubj));


pbar = 4/3;
xbin = -3.5:+3.5;
xvec = xbin;

subplot(2,2,1)
hold on;
xlim([-4,+4]);
ylim([0,1]);
for itask = 1:2
    patch([xvec,fliplr(xvec)],[yavghat(itask,:)+yerrhat(itask,:),fliplr(yavghat(itask,:)-yerrhat(itask,:))], ...
        0.5*(rgb(itask,:)+1),'EdgeColor','none');%fit
end
for itask = 1:2
    plot(xvec,yavghat(itask,:),'Color',rgb(itask,:),'LineWidth',1);%fit
end
for itask = 1:2
    errorbar(xvec,yavg(itask,:),yerr(itask,:),'k.','Marker','o',...
        'MarkerEdgeColor','k','LineStyle','none','LineWidth',1,...
        'MarkerFaceColor',rgb(itask,:),'MarkerSize',6); %data
end
plot([0,0],ylim,'k-');
plot(xlim,0.5*[1,1],'k--');
hold off
set(gca,'Layer','top','Box','off','TickDir','out','PlotBoxAspectRatio',[pbar,1,1]);
set(gca,'TickLength',[1,1]*0.02/max(pbar,1));
set(gca,'FontName','Helvetica','FontSize',8);
set(gca,'XTick',[-3.5:-0.5,+0.5:+3.5],'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'});
set(gca,'YTick',0:0.2:1);
xlabel('trial position from reversal','FontSize',8);
ylabel('fraction selecting correct','FontSize',8);
axes = findobj(gcf, 'type', 'axes');
for a = 1:length(axes)
    if axes(a).YColor <= [1 1 1]
        axes(a).YColor = [0 0 0];
    end
    if axes(a).XColor <= [1 1 1]
        axes(a).XColor = [0 0 0];
    end
end



% 2/ for confidence judgments
% fit
y = all_crev_fit;
yavghat = squeeze(mean(y,1)); % mean
% compute within-subject error bars (relevant for statistics)
y = bsxfun(@minus,y,squeeze(reshape(mean(pconf_fit,2),[1,1,nsubj])));
yerrhat = squeeze(std(y,[],1)/sqrt(nsubj)); % S.E.M.

% data
y = all_crev;
yavg = squeeze(mean(y,1));
% compute within-subject error bars (relevant for statistics)
y = bsxfun(@minus,y,squeeze(reshape(mean(pconf,2),[1,1,nsubj])));
yerr = squeeze(std(y,[],1)/sqrt(nsubj));


pbar = 4/3;
xbin = -3.5:+3.5;
xvec = xbin;

subplot(2,2,2)
hold on;
xlim([-4,+4]);
ylim([0.2,0.9]);
for itask = 1:2
    patch([xvec,fliplr(xvec)],[yavghat(itask,:)+yerrhat(itask,:),fliplr(yavghat(itask,:)-yerrhat(itask,:))], ...
        0.5*(rgb(itask,:)+1),'EdgeColor','none');%fit
end
for itask = 1:2
    plot(xvec,yavghat(itask,:),'Color',rgb(itask,:),'LineWidth',1);%fit
end
for itask = 1:2
    errorbar(xvec,yavg(itask,:),yerr(itask,:),'k.','Marker','o',...
        'MarkerEdgeColor','k','LineStyle','none','LineWidth',1,...
        'MarkerFaceColor',rgb(itask,:),'MarkerSize',6); %data
end
plot([0,0],ylim,'k-');
plot(xlim,mean(pconf(:))*[1,1],'k--');
hold off
set(gca,'Layer','top','Box','off','TickDir','out','PlotBoxAspectRatio',[pbar,1,1]);
set(gca,'TickLength',[1,1]*0.02/max(pbar,1));
set(gca,'FontName','Helvetica','FontSize',8);
set(gca,'XTick',[-3.5:-0.5,+0.5:+3.5],'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'});
set(gca,'YTick',-1:0.1:1);
xlabel('trial position from reversal','FontSize',8);
ylabel('fraction high confident','FontSize',8);
axes = findobj(gcf, 'type', 'axes');
for a = 1:length(axes)
    if axes(a).YColor <= [1 1 1]
        axes(a).YColor = [0 0 0];
    end
    if axes(a).XColor <= [1 1 1]
        axes(a).XColor = [0 0 0];
    end
end





% Plot response and confidence repetition curves

% 1/ for responses
% fit
yavghat = squeeze(mean(all_rrep_fit,1)); % mean
yerrhat = squeeze(std(all_rrep_fit,[],1)/sqrt(nsubj)); % S.E.M.
% data
yavg = squeeze(mean(all_rrep,1));
yerr = squeeze(std(all_rrep,[],1)/sqrt(nsubj));

pbar = 4/3;
xbin = -3.5:+3.5;
xvec = xbin;


subplot(2,2,3)
hold on;
xlim([-4,+4]);
ylim([0,1]);
for itask = 1:2
    patch([xvec,fliplr(xvec)],[yavghat(itask,:)+yerrhat(itask,:),fliplr(yavghat(itask,:)-yerrhat(itask,:))], ...
        0.5*(rgb(itask,:)+1),'EdgeColor','none');%fit
end
for itask = 1:2
    plot(xvec,yavghat(itask,:),'Color',rgb(itask,:),'LineWidth',1);%fit
end
for itask = 1:2
    errorbar(xvec,yavg(itask,:),yerr(itask,:),'k.','Marker','o',...
        'MarkerEdgeColor','k','LineStyle','none','LineWidth',1,...
        'MarkerFaceColor',rgb(itask,:),'MarkerSize',6); %data
end
plot([0,0],ylim,'k-');
plot(xlim,0.5*[1,1],'k--');
hold off
set(gca,'Layer','top','Box','off','TickDir','out','PlotBoxAspectRatio',[pbar,1,1]);
set(gca,'TickLength',[1,1]*0.02/max(pbar,1));
set(gca,'FontName','Helvetica','FontSize',8);
set(gca,'XTick',-4:2:+4,'XTickLabel',{'-4.0','-2.0','0','2.0','4.0'});
set(gca,'YTick',0:0.2:1);
xlabel('evidence direction (logLR)','FontSize',8);
ylabel('fraction repeat','FontSize',8);
axes = findobj(gcf, 'type', 'axes');
for a = 1:length(axes)
    if axes(a).YColor <= [1 1 1]
        axes(a).YColor = [0 0 0];
    end
    if axes(a).XColor <= [1 1 1]
        axes(a).XColor = [0 0 0];
    end
end



% 2/ for confidence
% fit
y = all_crep_fit;
yavghat = squeeze(mean(y,1)); % mean
% compute within-subject error bars (relevant for statistics)
y = bsxfun(@minus,y,squeeze(reshape(mean(pconf_fit,2),[1,1,nsubj])));
yerrhat = squeeze(std(y,[],1)/sqrt(nsubj)); % S.E.M.

% data
y = all_crep;
yavg = squeeze(mean(y,1));
% compute within-subject error bars (relevant for statistics)
y = bsxfun(@minus,y,squeeze(reshape(mean(pconf,2),[1,1,nsubj])));
yerr = squeeze(std(y,[],1)/sqrt(nsubj));

subplot(2,2,4)
hold on;
xlim([-4,+4]);
ylim([0.2,1]);%ylim([0,1]);
for itask = 1:2
    patch([xvec,fliplr(xvec)],[yavghat(itask,:)+yerrhat(itask,:),fliplr(yavghat(itask,:)-yerrhat(itask,:))], ...
        0.5*(rgb(itask,:)+1),'EdgeColor','none');%fit
end
for itask = 1:2
    plot(xvec,yavghat(itask,:),'Color',rgb(itask,:),'LineWidth',1);%fit
end
for itask = 1:2
    errorbar(xvec,yavg(itask,:),yerr(itask,:),'k.','Marker','o',...
        'MarkerEdgeColor','k','LineStyle','none','LineWidth',1,...
        'MarkerFaceColor',rgb(itask,:),'MarkerSize',6); %data
end
plot([0,0],ylim,'k-');
plot(xlim,0.5*[1,1],'k--');
hold off
set(gca,'Layer','top','Box','off','TickDir','out','PlotBoxAspectRatio',[pbar,1,1]);
set(gca,'TickLength',[1,1]*0.02/max(pbar,1));
set(gca,'FontName','Helvetica','FontSize',8);
set(gca,'XTick',-4:2:+4,'XTickLabel',{'-4.0','-2.0','0','2.0','4.0'});
set(gca,'YTick',0:0.2:1);
xlabel('evidence direction (logLR)','FontSize',8);
ylabel('fraction high confident','FontSize',8);
axes = findobj(gcf, 'type', 'axes');
for a = 1:length(axes)
    if axes(a).YColor <= [1 1 1]
        axes(a).YColor = [0 0 0];
    end
    if axes(a).XColor <= [1 1 1]
        axes(a).XColor = [0 0 0];
    end
end


figure(1);
fname = './Figure_model_validation';
print(fname,'-painters','-dpdf');



% Analysis of best-fitting parameters

% case without ssel fitted, comparison across conditions:
[~,p_h]    = ttest(h(:,1),h(:,2)) ;
[~,p_sinf] = ttest(sinf(:,1),sinf(:,2)) ;
[~,p_tcnf] = ttest(tcnf(:,1),tcnf(:,2)) ;
[~,p_scnf] = ttest(scnf(:,1),scnf(:,2)) ;
[~,p_gcnf] = ttest(gcnf(:,1),gcnf(:,2)) ;

figure(2)
% Note. selection noise is typically fixed at zero.

subplot(2,3,1)
hold on;
xlim([0.5,2.5]);
ylim([0,.5]);
errorbar(1:2,[mean(h(:,1)) 0],[std(h(:,1)) 0]/sqrt(nsubj), ...
    'Color',rgb(1,:),'LineStyle','none','LineWidth',3);
errorbar(1:2,[0 mean(h(:,2))],[0 std(h(:,2))]/sqrt(nsubj), ...
    'Color',rgb(2,:),'LineStyle','none','LineWidth',3);
title(['p = ',num2str(p_h)]);
ylabel('hazard rate');
set(gca,'FontName','Helvetica','FontSize',8);
set(gca,'XTick',1:2,'XTickLabel',{'C-','C+'});
xlabel('condition');
hold off

subplot(2,3,2)
hold on;
xlim([0.5,2.5]);
ylim([0,1]);
errorbar(1:2,[mean(sinf(:,1)) 0],[std(sinf(:,1)) 0]/sqrt(nsubj), ...
    'Color',rgb(1,:),'LineStyle','none','LineWidth',3);
errorbar(1:2,[0 mean(sinf(:,2))],[0 std(sinf(:,2))]/sqrt(nsubj), ...
    'Color',rgb(2,:),'LineStyle','none','LineWidth',3);
set(gca,'FontName','Helvetica','FontSize',8);
set(gca,'XTick',1:2,'XTickLabel',{'C-','C+'});
title(['p = ',num2str(round(p_sinf*100)/100)]);
ylabel('inference noise');
hold off

subplot(2,3,3)
hold on;
xlim([0.5,2.5]);
ylim([0,2.5]);
errorbar(1:2,[mean(tcnf(:,1)) 0],[std(tcnf(:,1)) 0]/sqrt(nsubj), ...
    'Color',rgb(1,:),'LineStyle','none','LineWidth',3);
errorbar(1:2,[0 mean(tcnf(:,2))],[0 std(tcnf(:,2))]/sqrt(nsubj), ...
    'Color',rgb(2,:),'LineStyle','none','LineWidth',3);
set(gca,'FontName','Helvetica','FontSize',8);
set(gca,'XTick',1:2,'XTickLabel',{'C-','C+'});
title(['p = ',num2str(p_tcnf)]);
ylabel('confidence threshold');
hold off

subplot(2,3,4)
hold on;
xlim([0.5,2.5]);
errorbar(1:2,[mean(scnf(:,1)) 0],[std(scnf(:,1)) 0]/sqrt(nsubj), ...
    'Color',rgb(1,:),'LineStyle','none','LineWidth',3);
errorbar(1:2,[0 mean(scnf(:,2))],[0 std(scnf(:,2))]/sqrt(nsubj), ...
    'Color',rgb(2,:),'LineStyle','none','LineWidth',3);
title(['p = ',num2str(round(p_scnf*100)/100)]);
ylabel('metacognitive noise');
set(gca,'FontName','Helvetica','FontSize',8);
set(gca,'XTick',1:2,'XTickLabel',{'C-','C+'});
xlabel('condition');
hold off

subplot(2,3,5)
hold on;
xlim([0.5,2.5]);
errorbar(1:2,[mean(gcnf(:,1)) 0],[std(gcnf(:,1)) 0]/sqrt(nsubj), ...
    'Color',rgb(1,:),'LineStyle','none','LineWidth',3);
errorbar(1:2,[0 mean(gcnf(:,2))],[0 std(gcnf(:,2))]/sqrt(nsubj), ...
    'Color',rgb(2,:),'LineStyle','none','LineWidth',3);
title(['p = ',num2str(round(p_gcnf*100)/100)]);
ylabel('confidence gain');
set(gca,'FontName','Helvetica','FontSize',8);
set(gca,'XTick',1:2,'XTickLabel',{'C-','C+'});
xlabel('condition');
hold off

figure(2);
fname = './Figure_bestfitting_params';
print(fname,'-painters','-dpdf');
