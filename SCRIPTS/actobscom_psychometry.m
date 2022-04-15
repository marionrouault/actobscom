% Rouault et al. (2022) BioRxiv. Controllability boosts neural and
% cognitive correlates of changes-of-mind in uncertain environments.
% This script reproduced the main behavioral and psychometric analyses and
% statistics reported in the paper.

% 4 main psychometric analyses:
% - Choice reversal curves
% - Confidence reversal curve
% - Choice repetition curve
% - Confidence repetition curve

% Fit of the choice-PSE on the next trial

% Performance, confidence, and discrimination

% task 1 = uncontrollable C- condition
% task 2 = controllable C+ condition


% set experiment referring to the Experiments presented in the paper:
expename = 'Experiment_1_2A'; % see the three possibilities below

% set list of subjects
switch expename
    case 'Experiment_1'
        % 1 session = 576 trials = 288 task 1 + 288 task 2
        subjlist = setdiff(01:17,[04,16]);
    case 'Experiment_2A'
        subjlist = setdiff(01:20,[06,15]);
    case 'Experiment_1_2A' % pooling Exp1 and Exp2A
        subjlistC = setdiff(01:17,[04,16]);
        subjlistD = setdiff(01:20,[06,15]);
        
        subjlist = [subjlistC subjlistD];
    otherwise
        error('undefined experiment name!');
end

nsubj = numel(subjlist);



% sequence position offsets for response reversal curve
ioff = -4:+3; % from first trial following reversal
npos = numel(ioff);

% evidence direction bins for response repetition curve
emin = [-inf,  -3,  -2,  -1,   0,  +1,  +2,  +3];
emax = [  -3,  -2,  -1,   0,  +1,  +2,  +3,+inf];
nbin = numel(emax);

paccu = nan(nsubj,2); % mean response accuracy
pconf = nan(nsubj,2); % mean fraction high confident
discrim = nan(nsubj,2); % mean confidence on correct vs. error responses

prev_resp = cell(nsubj,2); % response reversal curve
prev_conf = cell(nsubj,2); % confidence associated with response reversal curve
prep_resp = cell(nsubj,2); % response repetition curve
prep_conf = cell(nsubj,2); % confidence associated with response repetition curve

tcst_resp = nan(nsubj,2); % response reversal time constant
pmax_resp = nan(nsubj,2); % asymptotic response accuracy following reversal
tcst_conf = nan(nsubj,2); % confidence reversal time constant
pdrp_conf = nan(nsubj,2); % confidence drop following reversal

pse_resp = nan(nsubj,2); % point of subjective equivalence for response repetitions
woe_resp = nan(nsubj,2); % weight of evidence for responses
lap_resp = nan(nsubj,2); % proportion of blind response repetitions (lapses)

% same analysis as repetition curve for after repeat/switch trials only:
pse_resp_af_rep = nan(nsubj,2);
woe_resp_af_rep = nan(nsubj,2);
pse_resp_af_swi = nan(nsubj,2);
woe_resp_af_swi = nan(nsubj,2);

% same analysis as repetition curve for after high/low conf trials only:
pse_resp_af_hi = nan(nsubj,2);
woe_resp_af_hi = nan(nsubj,2);
pse_resp_af_lo = nan(nsubj,2);
woe_resp_af_lo = nan(nsubj,2);



for isubj = 1:nsubj
    
    disp(['Processing subject ',num2str(subjlist(isubj))])
    
    for itask = 1:2
        
        % load experiment data
        switch expename
            case 'Experiment_1'
                fname = sprintf('%s/%s_S%02d_task%d_expdata.mat','./DATA','ACTOBS_C',subjlist(isubj),itask);
                load(fname,'dat');
            case 'Experiment_2A'
                fname = sprintf('%s/%s_S%02d_task%d_expdata.mat','./DATA','ACTOBS_D_rule1',subjlist(isubj),itask);
                load(fname,'dat');
            case 'Experiment_1_2A'
                if isubj <= length(subjlistC)
                    fname = sprintf('%s/%s_S%02d_task%d_expdata.mat','./DATA','ACTOBS_C',subjlist(isubj),itask);
                else
                    fname = sprintf('%s/%s_S%02d_task%d_expdata.mat','./DATA','ACTOBS_D_rule1',subjlist(isubj),itask);
                end
                load(fname,'dat');
        end
        
        
        % recreate additional information
        blkind = 0; % current block index
        seqdir = 0; % current "true" sequence direction (1 or 2)
        epinum = 0; % episode number in current block, reset at break
        seqpos = 0; % sequence position in current episode, reset at reversal
        dat.epinum = nan(size(dat.seqdir));
        dat.seqpos = nan(size(dat.seqdir));
        dat.seqlen = nan(size(dat.seqdir));
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
            dat.seqlen(iseq) = length(dat.smpang{iseq}) ;
        end
        
        % perform simulations or retrieve data:
        cfg = [];
        cfg.seqind = dat.seqind;
        cfg.seqdir = dat.seqdir;
        cfg.seqlen = dat.seqlen;
        cfg.epinum = dat.epinum;
        cfg.seqpos = dat.seqpos;
        cfg.itask = itask;
        % sum the samples within each sequence to get the whole llr:
        cfg.seqllr = cellfun(@sum,dat.smpllr);
        
        
        % compute response accuracy and confidence
        paccu(isubj,itask) = mean(dat.raft == dat.seqdir); % fraction correct
        pconf(isubj,itask) = mean(dat.caft == 2); % fraction confident [1 low conf, 2 high conf]
        iscor = dat.raft == dat.seqdir ;
        conf_cor = dat.caft(iscor == 1) ;
        conf_err = dat.caft(iscor == 0) ;
        discrim(isubj,itask) = mean(conf_cor==2)-mean(conf_err==2) ;
        
        
        
        
        % ========= compute reversal curves =========
        
        i1st = find(dat.taskid == itask & dat.seqind > 1 & dat.seqpos == 1);
        for ipos = 1:npos
            prev_resp{isubj,itask}(ipos) = mean(dat.raft(i1st+ioff(ipos)) == dat.seqdir(i1st));
            prev_conf{isubj,itask}(ipos) = mean(dat.caft(i1st+ioff(ipos)) == 2);
        end
        
        % Fit exponential learning model to responses and confidence
        cfg.raft = dat.raft; % response provided after the current sequence
        cfg.caft = dat.caft;
        
        % Fit responses
        out = fit_model_prev(cfg);
        tcst_resp(isubj,itask) = out.tcst; % time constant
        pmax_resp(isubj,itask) = out.pmax; % asymptotic accuracy
        rhat = out.pr; % psychometric predictions
        
        
        % Fit confidence judgments
        out = fit_model_crev(cfg);
        % NB: "crev" fits one more param than "prev": "pmin", for computing confidence drop as extra variable
        tcst_conf(isubj,itask) = out.tcst; % time constant
        pdrp_conf(isubj,itask) = out.pdrp; % confidence drop
        chat = out.pr(:,1); % psychometric predictions
        
        % compute psychometric predictions for responses and
        % confidence judgments
        for ipos = 1:npos
            i1 = find(dat.seqdir(i1st) == 1);
            prev_resp_hat{isubj,itask}(ipos) = mean(rhat(sub2ind(size(rhat),i1st+ioff(ipos),dat.seqdir(i1st))));
            prev_conf_hat{isubj,itask}(ipos) = mean(chat(i1st+ioff(ipos)));
        end
        
        
        
        
        % ========= compute repetition curves =========
        
        ifilt = find(dat.seqind > 1);
        raft = dat.raft(ifilt);
        rbef = dat.rbef(ifilt); % previous response: gives you repeat when compared w raft
        caft = dat.caft(ifilt);
        slen = cellfun(@length,dat.smpllr(ifilt)); % llr of each sample within a sequence
        edir = cellfun(@sum,dat.smpllr(ifilt)).*(3-2*rbef);
        rdir = raft == rbef;
        cdir = caft == 2;
        
        % fraction confident
        pconfhat(isubj,itask) = mean(dat.caft==2);
        
        
        for ibin = 1:nbin
            ifilt = edir > emin(ibin) & edir < emax(ibin);
            prep_resp{isubj,itask}(ibin) = mean(rdir(ifilt));
            prep_conf{isubj,itask}(ibin) = mean(cdir(ifilt));
        end
        
        % fit response repetition curve with 2 or 3-param sigmoid function)
        % (subplot3):
        x = edir(:);
        y = rdir(:);
        
        [b,~,stat] = glmfit(x,y,'binomial','link','logit');
        if glmval(b,-10,'logit') > 0.5
            pse_resp(isubj,itask) = -inf;
        elseif glmval(b,+10,'logit') < 0.5
            pse_resp(isubj,itask) = +inf;
        else
            pse_resp(isubj,itask) = -fzero(@(L)glmval(b,L,'logit')-0.5,[-10,+10]);
        end
        woe_resp(isubj,itask) = b(2);
        lap_resp(isubj,itask) = 0;
        
        
        
        pr_hat = glmval(b,x,'logit');
        for ibin = 1:nbin
            ifilt = edir > emin(ibin) & edir < emax(ibin);
            prep_resp_hat{isubj,itask}(ibin) = mean(pr_hat(ifilt));
        end
        
        
        % Here, also fit pse and woe for subsamples of trials:
        rep_filt = zeros(1,length(rdir));    % 1 after a repeat trial, 0 after a switch trial
        hi_filt = zeros(1,length(rdir));     % 1 after a high conf trial, 0 after a low conf trial
        rep_filt_hi = zeros(1,length(rdir)); % 1 after a repeat trial with high conf, 0 after a repeat trial with low conf
        swi_filt_hi = zeros(1,length(rdir)); % 1 after a switch trial with high conf, 0 after a switch trial with low conf
        
        for irep = 2:length(rdir)
            if rdir(irep-1) == 1
                rep_filt(irep) = 1;
                if caft(irep-1) == 2
                    rep_filt_hi(irep) = 1;
                    hi_filt(irep) = 1;
                elseif caft(irep-1) == 1
                    rep_filt_hi(irep) = 0;
                    hi_filt(irep) = 0;
                end
            elseif rdir(irep-1) == 0
                rep_filt(irep) = 0;
                if caft(irep-1) == 2
                    swi_filt_hi(irep) = 1;
                    hi_filt(irep) = 1;
                elseif caft(irep-1) == 1
                    swi_filt_hi(irep) = 0;
                    hi_filt(irep) = 0;
                end
            end
        end
        
        % PSE fitted separately for after a repeat and after a switch:
        x = edir(rep_filt==1);
        y = rdir(rep_filt==1)';
        [b,~,stat] = glmfit(x,y,'binomial','link','logit');
        if glmval(b,-10,'logit') > 0.5
            pse_resp_af_rep(isubj,itask) = -inf;
        elseif glmval(b,+10,'logit') < 0.5
            pse_resp_af_rep(isubj,itask) = +inf;
        else
            pse_resp_af_rep(isubj,itask) = -fzero(@(L)glmval(b,L,'logit')-0.5,[-10,+10]);
        end
        woe_resp_af_rep(isubj,itask) = b(2);
        
        x = edir(rep_filt==0);
        y = rdir(rep_filt==0)';
        [b,~,stat] = glmfit(x,y,'binomial','link','logit');
        if glmval(b,-10,'logit') > 0.5
            pse_resp_af_swi(isubj,itask) = -inf;
        elseif glmval(b,+10,'logit') < 0.5
            pse_resp_af_swi(isubj,itask) = +inf;
        else
            pse_resp_af_swi(isubj,itask) = -fzero(@(L)glmval(b,L,'logit')-0.5,[-10,+10]);
        end
        woe_resp_af_swi(isubj,itask) = b(2);
        
        
        % PSE fitted separately for after a high / a low confidence trial:
        x = edir(hi_filt==1);
        y = rdir(hi_filt==1)';
        [b,~,stat] = glmfit(x,y,'binomial','link','logit');
        if glmval(b,-10,'logit') > 0.5
            pse_resp_af_hi(isubj,itask) = -inf;
        elseif glmval(b,+10,'logit') < 0.5
            pse_resp_af_hi(isubj,itask) = +inf;
        else
            pse_resp_af_hi(isubj,itask) = -fzero(@(L)glmval(b,L,'logit')-0.5,[-10,+10]);
        end
        woe_resp_af_hi(isubj,itask) = b(2);
        
        x = edir(hi_filt==0);
        y = rdir(hi_filt==0)';
        [b,~,stat] = glmfit(x,y,'binomial','link','logit');
        if glmval(b,-10,'logit') > 0.5
            pse_resp_af_lo(isubj,itask) = -inf;
        elseif glmval(b,+10,'logit') < 0.5
            pse_resp_af_lo(isubj,itask) = +inf;
        else
            pse_resp_af_lo(isubj,itask) = -fzero(@(L)glmval(b,L,'logit')-0.5,[-10,+10]);
        end
        woe_resp_af_lo(isubj,itask) = b(2);
        
        
        
        
        % Predicting fitted confidence associated with response repetition
        % curve by mixture of two 2-param sigmoid functions
        % weighted by the choice repetition curve:
        
        
        % (1/4) repeat sigmoid
        
        % select repeat trials
        x = edir(rdir==1)';
        y = cdir(rdir==1)';
        
        % perform logistic regression
        [b_rep,~,~] = glmfit(x,y,'binomial','link','logit');
        if glmval(b_rep,-10,'logit') > 0.5
            pse_rep(isubj,itask) = -inf;
        elseif glmval(b_rep,+10,'logit') < 0.5
            pse_rep(isubj,itask) = +inf;
        else
            pse_rep(isubj,itask) = -fzero(@(L)glmval(b_rep,L,'logit')-0.5,[-10,+10]);
        end
        woe_rep(isubj,itask) = b_rep(2);
        
        % re-apply but on all trials, not on repeat trials only:
        pr_rep = glmval(b_rep,edir(:),'logit');
        
        
        % (2/4) switch sigmoid
        
        % select switch trials
        x = edir(rdir==0)';
        y = cdir(rdir==0)';
        
        % perform logistic regression
        [b_swi,~,~] = glmfit(x,y,'binomial','link','logit');
        if glmval(b_swi,-10,'logit') < 0.5
            pse_swi(isubj,itask) = -inf;
        elseif glmval(b_swi,+10,'logit') > 0.5
            pse_swi(isubj,itask) = +inf;
        else
            pse_swi(isubj,itask) = -fzero(@(L)glmval(b_swi,L,'logit')-0.5,[-10,+10]);
        end
        woe_swi(isubj,itask) = b_swi(2);
        
        % re-apply but on all trials, not on switch trials only:
        pr_swi = glmval(b_swi,edir(:),'logit');
        
        
        % (3/4) mixture of rep sigmoid and swi sigmoid weighted by
        % repetition curve, trial by trial:
        conf_pr_hat = pr_hat.*pr_rep + (1-pr_hat).*pr_swi;
        
        
        % (4/4) binning for the plot by evidence levels (subplot4):
        for ibin = 1:nbin
            ifilt = edir > emin(ibin) & edir < emax(ibin);
            prep_conf_hat{isubj,itask}(ibin) = mean(conf_pr_hat(ifilt));
        end
        
        
        
        
    end % end of loop over itask (i.e. condition C-/C+)
    
end % end of loop over participants





%%%% ------- Statistical tests for psychometric parameters ------- %%%%


% Choice reversal curves
disp(['- Choice reversal curves:'])

[~,pval,~,stats] = ttest(tcst_resp(:,1),tcst_resp(:,2));
disp(['choice: time constant between C- and C+ t=',num2str(stats.tstat), ...
    ' df=',num2str(stats.df),' p=',num2str(pval)])

[~,pval,~,stats] = ttest(pmax_resp(:,1),pmax_resp(:,2));
disp(['choice: asymptotic accuracy between C- and C+ t=',num2str(stats.tstat), ...
    ' df=',num2str(stats.df),' p=',num2str(pval)])


% Confidence reversal curves
disp(['- Confidence reversal curves:'])

[~,pval,~,stats] = ttest(tcst_conf(:,1),tcst_conf(:,2));
disp(['confidence: time constant between C- and C+ t=',num2str(stats.tstat), ...
    ' df=',num2str(stats.df),' p=',num2str(pval)])

[~,pval,~,stats] = ttest(pdrp_conf(:,1),pdrp_conf(:,2));
disp(['confidence: confidence drop between C- and C+ t=',num2str(stats.tstat), ...
    ' df=',num2str(stats.df),' p=',num2str(pval)])


% Choice repetition curves
disp(['- Choice repetition curves:'])

[~,pval,~,stats] = ttest(pse_resp(:,1),pse_resp(:,2));
disp(['choice: PSE between C- and C+ t=',num2str(stats.tstat), ...
    ' df=',num2str(stats.df),' p=',num2str(pval)])

[~,pval,~,stats] = ttest(woe_resp(:,1),woe_resp(:,2));
disp(['choice: weight of evidrnce between C- and C+ t=',num2str(stats.tstat), ...
    ' df=',num2str(stats.df),' p=',num2str(pval)])


% Confidence repetition curves, each sigmoid has 2 params:
disp(['- Confidence repetition curves:'])

[~,pval,~,stats] = ttest(pse_rep(:,1),pse_rep(:,2));
disp(['confidence: pse_rep between C- and C+ t=',num2str(stats.tstat), ...
    ' df=',num2str(stats.df),' p=',num2str(pval)])

[~,pval,~,stats] = ttest(woe_rep(:,1),woe_rep(:,2));
disp(['confidence: woe_rep between C- and C+ t=',num2str(stats.tstat), ...
    ' df=',num2str(stats.df),' p=',num2str(pval)])

[~,pval,~,stats] = ttest(pse_swi(:,1),pse_swi(:,2));
disp(['confidence: pse_swi between C- and C+ t=',num2str(stats.tstat), ...
    ' df=',num2str(stats.df),' p=',num2str(pval)])
% => these will be performed using jackknifing

[~,pval,~,stats] = ttest(woe_swi(:,1),woe_swi(:,2));
disp(['confidence: woe_swi between C- and C+ t=',num2str(stats.tstat), ...
    ' df=',num2str(stats.df),' p=',num2str(pval)])

disp 'Note: jackknifed statistics will be necessary for some of the parameters'



%% plot reversal and repetition curves

rgb = [[0,0.5,1];[1,0,0.5]];
rgb = 0.75*rgb+0.25;

% offset individual data points on x-axis
sss = 0.17 ;
ttt = 0.27 ;



% -------------- plot reversal curves --------------


% 1/ for responses

% fit
yhat = permute(reshape(cat(2,prev_resp_hat{:}),[8,nsubj,2]),[3,1,2]);
yavghat = mean(yhat,3); % mean
yerrhat = std(yhat,[],3)/sqrt(nsubj); % S.E.M.
% data
y = permute(reshape(cat(2,prev_resp{:}),[8,nsubj,2]),[3,1,2]);
yavg = mean(y,3); % mean
yerr = std(y,[],3)/sqrt(nsubj); % S.E.M.

pbar = 4/3;
xbin = -3.5:+3.5;
xvec = xbin;


figure(1)
subplot(2,2,1)
hold on;
xlim([-4,+4]);
ylim([0,1]);
for itask = 1:2
    patch([xvec,fliplr(xvec)],[yavghat(itask,:)+yerrhat(itask,:),fliplr(yavghat(itask,:)-yerrhat(itask,:))], ...
        0.5*(rgb(itask,:)+1),'EdgeColor','none'); %fit
end
for itask = 1:2
    plot(xvec,yavghat(itask,:),'Color',rgb(itask,:),'LineWidth',1.5); %fit
    
    errorbar(xvec,yavg(itask,:),yerr(itask,:),'k.','Marker','o',...
        'MarkerEdgeColor','k','LineStyle','none','LineWidth',1.5,...
        'MarkerFaceColor',rgb(itask,:),'MarkerSize',6); %data
end
plot([0,0],ylim,'k-');
plot(xlim,0.5*[1,1],'k--');
hold off
set(gca,'Layer','top','Box','off','TickDir','out','PlotBoxAspectRatio',[pbar,1,1]);
set(gca,'TickLength',[1,1]*0.02/max(pbar,1));
set(gca,'FontName','Helvetica','FontSize',12);
set(gca,'XTick',[-3.5:-0.5,+0.5:+3.5],'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'});
set(gca,'YTick',0:0.2:1);
xlabel('trial position from reversal','FontSize',12);
ylabel('fraction reported correct hidden state','FontSize',12);
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

%fit
yhat = permute(reshape(cat(2,prev_conf_hat{:}),[8,nsubj,2]),[3,1,2]);
yavghat = mean(yhat,3); % mean
yhat = bsxfun(@minus,yhat,reshape(mean(pconfhat,2),[1,1,nsubj]));
yerrhat = std(yhat,[],3)/sqrt(nsubj); % S.E.M.

%data
y = permute(reshape(cat(2,prev_conf_hat{:}),[8,nsubj,2]),[3,1,2]);
yavg = mean(y,3);
% compute within-subject error bars (relevant for statistics)
y = bsxfun(@minus,y,reshape(mean(pconf,2),[1,1,nsubj]));
yerr = std(y,[],3)/sqrt(nsubj);

pbar = 4/3;
xbin = -3.5:+3.5;
xvec = xbin;

figure(1)
subplot(2,2,2)
hold on;
xlim([-4,+4]);
ylim([0.4,0.7]);
for itask = 1:2
    patch([xvec,fliplr(xvec)],[yavghat(itask,:)+yerrhat(itask,:),fliplr(yavghat(itask,:)-yerrhat(itask,:))], ...
        0.5*(rgb(itask,:)+1),'EdgeColor','none'); %fit
end
for itask = 1:2
    plot(xvec,yavghat(itask,:),'Color',rgb(itask,:),'LineWidth',1.5); %fit
    errorbar(xvec,yavg(itask,:),yerr(itask,:),'k.','Marker','o',...
        'MarkerEdgeColor','k','LineStyle','none','LineWidth',1.5,...
        'MarkerFaceColor',rgb(itask,:),'MarkerSize',6); %data
end
plot([0,0],ylim,'k-');
plot(xlim,mean(pconf(:))*[1,1],'k--');
hold off
set(gca,'Layer','top','Box','off','TickDir','out','PlotBoxAspectRatio',[pbar,1,1]);
set(gca,'TickLength',[1,1]*0.02/max(pbar,1));
set(gca,'FontName','Helvetica','FontSize',12);
set(gca,'XTick',[-3.5:-0.5,+0.5:+3.5],'XTickLabel',{'-4','-3','-2','-1','1','2','3','4'});
set(gca,'YTick',-1:0.1:1);
xlabel('trial position from reversal','FontSize',12);
ylabel('fraction high confident','FontSize',12);
axes = findobj(gcf, 'type', 'axes');
for a = 1:length(axes)
    if axes(a).YColor <= [1 1 1]
        axes(a).YColor = [0 0 0];
    end
    if axes(a).XColor <= [1 1 1]
        axes(a).XColor = [0 0 0];
    end
end





% -------------- plot repetition curves --------------


% 1/ for responses
% fit
yhat = permute(reshape(cat(2,prep_resp_hat{:}),[8,nsubj,2]),[3,1,2]);
yavghat = mean(yhat,3); % mean
yerrhat = std(yhat,[],3)/sqrt(nsubj); % S.E.M.
% data
y = permute(reshape(cat(2,prep_resp{:}),[8,nsubj,2]),[3,1,2]);
yavg = mean(y,3);
yerr = std(y,[],3)/sqrt(nsubj);

pbar = 4/3;
xbin = -3.5:+3.5;
xvec = xbin;
figure(1)
subplot(2,2,3)
hold on;
xlim([-4,+4]);
ylim([0,1]);
for itask = 1:2
    patch([xvec,fliplr(xvec)],[yavghat(itask,:)+yerrhat(itask,:),fliplr(yavghat(itask,:)-yerrhat(itask,:))], ...
        0.5*(rgb(itask,:)+1),'EdgeColor','none');%fit
end
for itask = 1:2
    plot(xvec,yavghat(itask,:),'Color',rgb(itask,:),'LineWidth',2);%fit
    errorbar(xvec,yavg(itask,:),yerr(itask,:),'k.','Marker','o',...
        'MarkerEdgeColor','k','LineStyle','none','LineWidth',1.5,...
        'MarkerFaceColor',rgb(itask,:),'MarkerSize',6); %data
end
plot([0,0],ylim,'k-');
plot(xlim,0.5*[1,1],'k--');
hold off
set(gca,'Layer','top','Box','off','TickDir','out','PlotBoxAspectRatio',[pbar,1,1]);
set(gca,'TickLength',[1,1]*0.02/max(pbar,1));
set(gca,'FontName','Helvetica','FontSize',12);
set(gca,'XTick',-4:2:+4,'XTickLabel',{'-4.0','-2.0','0','2.0','4.0'});
set(gca,'YTick',0:0.2:1);
xlabel('evidence direction (logLR)','FontSize',12);
ylabel('fraction repeat','FontSize',12);
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

%fit
yhat = permute(reshape(cat(2,prep_conf_hat{:}),[8,nsubj,2]),[3,1,2]);
yavghat = mean(yhat,3);
% compute error bars
yerrhat = std(yhat,[],3)/sqrt(nsubj); % S.E.M.

%data
y = permute(reshape(cat(2,prep_conf{:}),[8,nsubj,2]),[3,1,2]);
yavg = mean(y,3);
% compute within-subject error bars (relevant for statistics)
y = bsxfun(@minus,y,reshape(mean(pconf,2),[1,1,nsubj]));
yerr = std(y,[],3)/sqrt(nsubj);

pbar = 4/3;
xbin = -3.5:+3.5;
xvec = xbin;
figure(1)
subplot(2,2,4)
hold on;
xlim([-4,+4]);
ylim([0.2,1]);
for itask = 1:2
    patch([xvec,fliplr(xvec)],[yavghat(itask,:)+yerrhat(itask,:),fliplr(yavghat(itask,:)-yerrhat(itask,:))], ...
        0.5*(rgb(itask,:)+1),'EdgeColor','none');%fit
end
for itask = 1:2
    plot(xvec,yavghat(itask,:),'Color',rgb(itask,:),'LineWidth',2);%fit
    errorbar(xvec,yavg(itask,:),yerr(itask,:),'k.','Marker','o',...
        'MarkerEdgeColor','k','LineStyle','none','LineWidth',1.5,...
        'MarkerFaceColor',rgb(itask,:),'MarkerSize',6); %data
end
plot([0,0],ylim,'k-');
plot(xlim,mean(pconf(:))*[1,1],'k--');
hold off
set(gca,'Layer','top','Box','off','TickDir','out','PlotBoxAspectRatio',[pbar,1,1]);
set(gca,'TickLength',[1,1]*0.02/max(pbar,1));
set(gca,'FontName','Helvetica','FontSize',12);
set(gca,'XTick',-4:2:+4,'XTickLabel',{'-4.0','-2.0','0','2.0','4.0'});
set(gca,'YTick',0:0.2:1);
xlabel('evidence direction (logLR)','FontSize',12);
ylabel('fraction high confident','FontSize',12);
axes = findobj(gcf, 'type', 'axes');
for a = 1:length(axes)
    if axes(a).YColor <= [1 1 1]
        axes(a).YColor = [0 0 0];
    end
    if axes(a).XColor <= [1 1 1]
        axes(a).XColor = [0 0 0];
    end
end







%% Model-free statistics for performance, confidence, and discrimination

disp 'Statistics for basic behaviour:'

[~,p_paccu,~,stats_accu] = ttest(paccu(:,1),paccu(:,2)) ;
disp(['mean accuracy between C- and C+ t=',num2str(stats_accu.tstat), ...
    ' df=',num2str(stats_accu.df),' p=',num2str(p_paccu)])

[~,p_pconf,~,stats_pconf] = ttest(pconf(:,1),pconf(:,2)) ;
disp(['mean confidence between C- and C+ t=',num2str(stats_pconf.tstat), ...
    ' df=',num2str(stats_pconf.df),' p=',num2str(p_pconf)])

[~,p_discrim,~,stats_discrim] = ttest(discrim(:,1),discrim(:,2)) ;
disp(['discrimination capacity between C- and C+ t=',num2str(stats_discrim.tstat), ...
    ' df=',num2str(stats_discrim.df),' p=',num2str(p_discrim)])



figure(2)

subplot(2,2,1)
hold on;
xlim([0.5,2.5]);
ylim([0.5,1]);
errorbar(1:2,[mean(paccu(:,1)) 0],[std(paccu(:,1)) 0]/sqrt(nsubj), ...
    'Color',rgb(1,:),'LineStyle','none','LineWidth',3);
errorbar(1:2,[0 mean(paccu(:,2))],[0 std(paccu(:,2))]/sqrt(nsubj), ...
    'Color',rgb(2,:),'LineStyle','none','LineWidth',3);
set(gca,'FontName','Helvetica','FontSize',12,'LineWidth',1.5);
set(gca,'XTick',1:2,'XTickLabel',{'C-','C+'});
xlabel('condition');
title(['p=',num2str(p_paccu)]);
ylabel('performance');
hold off

subplot(2,2,2)
hold on;
xlim([0.5,2.5]);
ylim([0,1]);
errorbar(1:2,[mean(pconf(:,1)) 0],[std(pconf(:,1)) 0]/sqrt(nsubj), ...
    'Color',rgb(1,:),'LineStyle','none','LineWidth',3);
errorbar(1:2,[0 mean(pconf(:,2))],[0 std(pconf(:,2))]/sqrt(nsubj), ...
    'Color',rgb(2,:),'LineStyle','none','LineWidth',3);
set(gca,'FontName','Helvetica','FontSize',12,'LineWidth',1.5);
set(gca,'XTick',1:2,'XTickLabel',{'C-','C+'});
xlabel('condition');
title(['p=',num2str(p_pconf)]);
ylabel('proportion high confidence');
hold off

subplot(2,2,3)
hold on;
xlim([0.5,2.5]);
ylim([0,1]);
errorbar(1:2,[mean(discrim(:,1)) 0],[std(discrim(:,1)) 0]/sqrt(nsubj), ...
    'Color',rgb(1,:),'LineStyle','none','LineWidth',3);
errorbar(1:2,[0 mean(discrim(:,2))],[0 std(discrim(:,2))]/sqrt(nsubj), ...
    'Color',rgb(2,:),'LineStyle','none','LineWidth',3);
set(gca,'FontName','Helvetica','FontSize',12,'LineWidth',1.5);
set(gca,'XTick',1:2,'XTickLabel',{'C-','C+'});
xlabel('condition');
title(['p=',num2str(p_discrim)]);
ylabel('discrimination');
hold off





%% pse and woe repetition curves: parameters for subsamples of trials
% (related to the main choice repetition curve)

figure(3)

% after a repeat vs. after a switch trial
subplot(2,2,1)
hold on;
bar(1, mean(pse_resp_af_rep(:,1)),'FaceColor',rgb(1,:),'EdgeColor','k','LineWidth',2,'LineStyle','-')
hold on
bar(2, mean(pse_resp_af_swi(:,1)),'FaceColor',rgb(1,:),'EdgeColor','k','LineWidth',2,'LineStyle','-')
hold on
bar(3, mean(pse_resp_af_rep(:,2)),'FaceColor',rgb(2,:),'EdgeColor','k','LineWidth',2,'LineStyle','-')
hold on
bar(4, mean(pse_resp_af_swi(:,2)),'FaceColor',rgb(2,:),'EdgeColor','k','LineWidth',2,'LineStyle','-')
hold on
errorbar(1:4,mean([pse_resp_af_rep(:,1) pse_resp_af_swi(:,1) pse_resp_af_rep(:,2) pse_resp_af_swi(:,2)]), ...
    std([pse_resp_af_rep(:,1) pse_resp_af_swi(:,1) pse_resp_af_rep(:,2) pse_resp_af_swi(:,2)])/sqrt(nsubj), ...
    'Color','k','LineStyle','none','LineWidth',3);
for o = 1:nsubj
    hold on;
    plot([1+sss 1+ttt],[pse_resp_af_rep(o,1) pse_resp_af_rep(o,1)],'k','LineWidth',2) ;
    plot([2+sss 2+ttt],[pse_resp_af_swi(o,1) pse_resp_af_swi(o,1)],'k','LineWidth',2) ;
    plot([3+sss 3+ttt],[pse_resp_af_rep(o,2) pse_resp_af_rep(o,2)],'k','LineWidth',2) ;
    plot([4+sss 4+ttt],[pse_resp_af_swi(o,2) pse_resp_af_swi(o,2)],'k','LineWidth',2) ;
end
set(gca,'FontName','Helvetica','FontSize',12,'LineWidth',1.5, ...
    'XTick',1:4,'XTickLabel',{'after','after','after','after'})
xlabel('repeat switch repeat switch','fontsize',12)
ylabel('PSE','fontsize',12)
axis([0 5 -.3 3.7])
hold off


subplot(2,2,3)
hold on;
bar(1, mean(woe_resp_af_rep(:,1)),'FaceColor',rgb(1,:),'EdgeColor','k','LineWidth',2,'LineStyle','-')
hold on
bar(2, mean(woe_resp_af_swi(:,1)),'FaceColor',rgb(1,:),'EdgeColor','k','LineWidth',2,'LineStyle','-')
hold on
bar(3, mean(woe_resp_af_rep(:,2)),'FaceColor',rgb(2,:),'EdgeColor','k','LineWidth',2,'LineStyle','-')
hold on
bar(4, mean(woe_resp_af_swi(:,2)),'FaceColor',rgb(2,:),'EdgeColor','k','LineWidth',2,'LineStyle','-')
hold on

errorbar(1:4,mean([woe_resp_af_rep(:,1) woe_resp_af_swi(:,1) woe_resp_af_rep(:,2) woe_resp_af_swi(:,2)]), ...
    std([woe_resp_af_rep(:,1) woe_resp_af_swi(:,1) woe_resp_af_rep(:,2) woe_resp_af_swi(:,2)])/sqrt(nsubj), ...
    'Color','k','LineStyle','none','LineWidth',3);
for o = 1:nsubj
    hold on;
    plot([1+sss 1+ttt],[woe_resp_af_rep(o,1) woe_resp_af_rep(o,1)],'k','LineWidth',2) ;
    plot([2+sss 2+ttt],[woe_resp_af_swi(o,1) woe_resp_af_swi(o,1)],'k','LineWidth',2) ;
    plot([3+sss 3+ttt],[woe_resp_af_rep(o,2) woe_resp_af_rep(o,2)],'k','LineWidth',2) ;
    plot([4+sss 4+ttt],[woe_resp_af_swi(o,2) woe_resp_af_swi(o,2)],'k','LineWidth',2) ;
end
set(gca,'FontName','Helvetica','FontSize',12,'LineWidth',1.5, ...
    'XTick',1:4,'XTickLabel',{'after','after','after','after'})
xlabel('repeat switch repeat switch','fontsize',12)
ylabel('sensitivity to evidence','fontsize',12)
axis([0 5 .2 7.3])
hold off

% after a high vs. a low confidence trial:
subplot(2,2,2)
hold on;
bar(1, mean(pse_resp_af_hi(:,1)),'FaceColor',rgb(1,:),'EdgeColor','k','LineWidth',2,'LineStyle','-')
hold on
bar(2, mean(pse_resp_af_lo(:,1)),'FaceColor',rgb(1,:),'EdgeColor','k','LineWidth',2,'LineStyle','-')
hold on
bar(3, mean(pse_resp_af_hi(:,2)),'FaceColor',rgb(2,:),'EdgeColor','k','LineWidth',2,'LineStyle','-')
hold on
bar(4, mean(pse_resp_af_lo(:,2)),'FaceColor',rgb(2,:),'EdgeColor','k','LineWidth',2,'LineStyle','-')
hold on
errorbar(1:4,mean([pse_resp_af_hi(:,1) pse_resp_af_lo(:,1) pse_resp_af_hi(:,2) pse_resp_af_lo(:,2)]), ...
    std([pse_resp_af_hi(:,1) pse_resp_af_lo(:,1) pse_resp_af_hi(:,2) pse_resp_af_lo(:,2)])/sqrt(nsubj), ...
    'Color','k','LineStyle','none','LineWidth',3);
for o = 1:nsubj
    hold on;
    plot([1+sss 1+ttt],[pse_resp_af_hi(o,1) pse_resp_af_hi(o,1)],'k','LineWidth',2) ;
    plot([2+sss 2+ttt],[pse_resp_af_lo(o,1) pse_resp_af_lo(o,1)],'k','LineWidth',2) ;
    plot([3+sss 3+ttt],[pse_resp_af_hi(o,2) pse_resp_af_hi(o,2)],'k','LineWidth',2) ;
    plot([4+sss 4+ttt],[pse_resp_af_lo(o,2) pse_resp_af_lo(o,2)],'k','LineWidth',2) ;
end
set(gca,'FontName','Helvetica','FontSize',12,'LineWidth',1.5, ...
    'XTick',1:4,'XTickLabel',{'after','after','after','after'})
xlabel('high low high low','fontsize',12)
ylabel('PSE','fontsize',12)
axis([0 5 -.3 3.7])
hold off

subplot(2,2,4)
hold on;
bar(1, mean(woe_resp_af_hi(:,1)),'FaceColor',rgb(1,:),'EdgeColor','k','LineWidth',2,'LineStyle','-')
hold on
bar(2, mean(woe_resp_af_lo(:,1)),'FaceColor',rgb(1,:),'EdgeColor','k','LineWidth',2,'LineStyle','-')
hold on
bar(3, mean(woe_resp_af_hi(:,2)),'FaceColor',rgb(2,:),'EdgeColor','k','LineWidth',2,'LineStyle','-')
hold on
bar(4, mean(woe_resp_af_lo(:,2)),'FaceColor',rgb(2,:),'EdgeColor','k','LineWidth',2,'LineStyle','-')
hold on
errorbar(1:4,mean([woe_resp_af_hi(:,1) woe_resp_af_lo(:,1) woe_resp_af_hi(:,2) woe_resp_af_lo(:,2)]), ...
    std([woe_resp_af_hi(:,1) woe_resp_af_lo(:,1) woe_resp_af_hi(:,2) woe_resp_af_lo(:,2)])/sqrt(nsubj), ...
    'Color','k','LineStyle','none','LineWidth',3);
for o = 1:nsubj
    hold on;
    plot([1+sss 1+ttt],[woe_resp_af_hi(o,1) woe_resp_af_hi(o,1)],'k','LineWidth',2) ;
    plot([2+sss 2+ttt],[woe_resp_af_lo(o,1) woe_resp_af_lo(o,1)],'k','LineWidth',2) ;
    plot([3+sss 3+ttt],[woe_resp_af_hi(o,2) woe_resp_af_hi(o,2)],'k','LineWidth',2) ;
    plot([4+sss 4+ttt],[woe_resp_af_lo(o,2) woe_resp_af_lo(o,2)],'k','LineWidth',2) ;
end
set(gca,'FontName','Helvetica','FontSize',12,'LineWidth',1.5, ...
    'XTick',1:4,'XTickLabel',{'after','after','after','after'})
xlabel('high low high low','fontsize',12)
ylabel('sensitivity to evidence','fontsize',12)
axis([0 5 .2 7.3])
hold off





