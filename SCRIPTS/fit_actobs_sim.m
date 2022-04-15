function [out] = fit_actobs_sim(cfg)

% check configuration structure
if ~all(isfield(cfg,{'seqind','seqpos','seqdir','seqllr','seqlen'}))
    error('Missing experiment data!');
end
if ~all(isfield(cfg,{'rbef','raft'}))
    error('Missing response data!');
end
if ~all(isfield(cfg,{'cbef','caft'}))
    if cfg.fitcnf
        error('Missing confidence data!');
    end
    cfg.cbef(:) = nan;
    cfg.caft(:) = nan;
end
if ~isfield(cfg,'fitrev')
    cfg.fitrev = false;
end
if ~isfield(cfg,'fitrep')
    cfg.fitrep = false;
end
if ~isfield(cfg,'fitcnf')
    cfg.fitcnf = ...
        all(isfield(cfg,{'cbef','caft'})) && ...
        ~any(isnan(cfg.cbef)) && ~any(isnan(cfg.caft));
end
if ~isfield(cfg,'resamp')
    cfg.resamp = false;
end
if ~isfield(cfg,'nsmp')
    cfg.nsmp = 1e3;
end
if ~isfield(cfg,'nres')
    cfg.nres = 1e3;
end
if ~isfield(cfg,'nval')
    cfg.nval = 1e2;
end
if ~isfield(cfg,'nrun')
    cfg.nrun = 10;
end
if ~isfield(cfg,'verbose')
    cfg.verbose = 0;
end

% get experiment data
seqind = cfg.seqind(:); % sequence index in current block
seqpos = cfg.seqpos(:); % sequence position in current episode
seqdir = cfg.seqdir(:); % sequence direction
seqllr = cfg.seqllr(:); % sequence log-likelihood ratio
seqlen = cfg.seqlen(:); % sequence length

% get number of sequences
nseq = numel(seqind);

% get fitting configuration
fitrev  = cfg.fitrev; % fit reversal curves?
fitrep  = cfg.fitrep; % fit repetition curves?
fitcnf  = cfg.fitcnf; % fit confidence?
resamp  = cfg.resamp; % resample log-belief conditioned upon response?
nsmp    = cfg.nsmp; % number of samples used in simulations (fitalgo = vbmc)
nval    = cfg.nval; % number of validation samples (fitalgo = bads)
nrun    = cfg.nrun; % number of random starting points (fitalgo = bads)
verbose = cfg.verbose; % fitting display level

% set additional fitting parameters
lmin = 0.5/nsmp;

% do not fit confidence parameters
if ~fitcnf
    cfg.scnf = nan;
    cfg.tcnf = nan;
    cfg.gcnf = nan;
end

% get subject reversal and repetition curves
csub = getc(seqind,seqpos,seqdir,seqllr, ...
    cfg.rbef(:),cfg.raft(:),cfg.cbef(:),cfg.caft(:));

% define model parameters
pnam = {}; % name
pmin = []; % minimum value
pmax = []; % maximum value
pfun = {}; % log-prior function
pini = []; % initial value
pplb = []; % plausible lower bound
ppub = []; % plausible upper bound
% 1/ perceived hazard rate
pnam{1,1} = 'h';
pmin(1,1) = 1e-6;
pmax(1,1) = 1-(1e-6);
pfun{1,1} = @(x)betapdf(x,1,7);
pini(1,1) = 0.125;
pplb(1,1) = betainv(0.1587,1,7);
ppub(1,1) = betainv(0.8413,1,7);
% 2/ inference noise
pnam{1,2} = 'sinf';
pmin(1,2) = 0;
pmax(1,2) = 10;
pfun{1,2} = @(x)gampdf(x,1,0.5);
pini(1,2) = 0.5;
pplb(1,2) = gaminv(0.1587,1,0.5);
ppub(1,2) = gaminv(0.8413,1,0.5);
% 3/ selection noise
pnam{1,3} = 'ssel';
pmin(1,3) = 0;
pmax(1,3) = 10;
pfun{1,3} = @(x)gampdf(x,1,1);
pini(1,3) = 1;
pplb(1,3) = gaminv(0.1587,1,1);
ppub(1,3) = gaminv(0.8413,1,1);
% 4/ confidence noise
pnam{1,4} = 'scnf';
pmin(1,4) = 0;
pmax(1,4) = 10;
pfun{1,4} = @(x)gampdf(x,1,1);
pini(1,4) = 1;
pplb(1,4) = gaminv(0.1587,1,1);
ppub(1,4) = gaminv(0.8413,1,1);
% 5/ confidence threshold
pnam{1,5} = 'tcnf';
pmin(1,5) = -10;
pmax(1,5) = +10;
pfun{1,5} = @(x)normpdf(x,0,1);
pini(1,5) = 0;
pplb(1,5) = norminv(0.1587,0,1);
ppub(1,5) = norminv(0.8413,0,1);
% 6/ confidence gain during switches
pnam{1,6} = 'gcnf';
pmin(1,6) = 0;
pmax(1,6) = 10;
pfun{1,6} = @(x)gampdf(x,2,1);
pini(1,6) = 1;
pplb(1,6) = gaminv(0.1587,2,1);
ppub(1,6) = gaminv(0.8413,2,1);

% set number of parameters
npar = numel(pnam);

% define fixed parameters
pfix = cell(1,npar);
for i = 1:npar
    if isfield(cfg,pnam{i}) && ~isempty(cfg.(pnam{i}))
        pfix{i} = cfg.(pnam{i});
    end
end

% define free parameters
ifit = cell(1,npar);
pfit_ini = [];
pfit_min = [];
pfit_max = [];
pfit_plb = [];
pfit_pub = [];
n = 1;
for i = 1:npar
    if isempty(pfix{i}) % free parameter
        ifit{i} = n;
        pfit_ini = cat(2,pfit_ini,pini(i));
        pfit_min = cat(2,pfit_min,pmin(i));
        pfit_max = cat(2,pfit_max,pmax(i));
        pfit_plb = cat(2,pfit_plb,pplb(i));
        pfit_pub = cat(2,pfit_pub,ppub(i));
        n = n+1;
    end
end

ntrl = nseq; % number of trials
nfit = length(pfit_ini); % number of fitted parameters

if nfit > 0
    
    % fit model using Bayesian Adaptive Direct Search
    if ~exist('bads','file')
        error('BADS missing from path!');
    end
    
    % configure BADS
    options = bads('defaults');
    options.UncertaintyHandling = true; % noisy objective function
    options.NoiseFinalSamples = nval; % number of samples
    switch verbose % display level
        case 0, options.Display = 'none';
        case 1, options.Display = 'final';
        case 2, options.Display = 'iter';
    end
    
    % fit model using multiple random starting points
    fval   = nan(1,nrun);
    xhat   = cell(1,nrun);
    output = cell(1,nrun);
    for irun = 1:nrun
        done = false;
        while ~done
            % set random starting point
            n = 1;
            for i = 1:npar
                if isempty(pfix{i}) % free parameter
                    % sample starting point uniformly within plausible bounds
                    pfit_ini(n) = unifrnd(pplb(i),ppub(i));
                    n = n+1;
                end
            end
            % fit model using BADS
            [xhat{irun},fval(irun),exitflag,output{irun}] = ...
                bads(@(x)getnll(x), ...
                pfit_ini,pfit_min,pfit_max,pfit_plb,pfit_pub,[],options);
            if exitflag > 0
                done = true;
            end
        end
    end
    % find best fit among random starting points
    [fval,irun] = min(fval);
    xhat   = xhat{irun};
    output = output{irun};
    
    % get full parameter set with best-fitting values
    phat = getpval(xhat);
    
    % create output structure with best-fitting values
    out = cell2struct(phat(:),pnam(:));
    out = rmmissp(out); % remove missing parameters
    
    % store fitting information
    out.nsmp = nsmp; % number of samples used by particle filter
    out.nval = nval; % number of validation samples
    out.nrun = nrun; % number of random starting points
    out.ntrl = ntrl; % number of trials
    out.nfit = nfit; % number of fitted parameters
    
    % get maximum log-likelihood
    out.ll = -output.fval; % estimated log-likelihood
    out.ll_sd = output.fsd; % estimated s.d. of log-likelihood
    
    % get complexity-penalized fitting metrics
    out.aic = -2*out.ll+2*nfit+2*nfit*(nfit+1)/(ntrl-nfit+1); % AIC
    out.bic = -2*out.ll+nfit*log(ntrl); % BIC
    
    % get parameter values
    out.xnam = pnam(cellfun(@isempty,pfix));
    out.xhat = xhat;
    
    % store additional output from BADS
    out.output = output;
    
else
    
    % get full parameter set with fixed values
    phat = getpval([]);
    
    % create output structure with fixed parameter values
    out = cell2struct(phat(:),pnam(:));
    out = rmmissp(out); % remove missing parameters
    
end

% store configuration structure
out.cfg = cfg;

if resamp
    % store fraction of propagated samples
    [~,out.psmp] = sim_c(phat{:});
end

% store reversal and repetition curves
resamp = false;
out.csub = csub; % subject
out.cfit = sim_c(phat{:}); % best-fitting predictions

    function [pval] = getpval(p)
        % get parameter values
        pval = cell(1,npar);
        for k = 1:npar
            if isempty(pfix{k}) % free parameter
                pval{k} = p(ifit{k});
            else % fixed parameter
                pval{k} = pfix{k};
            end
        end
    end


    function [nll] = getnll(p)
        % get parameter values
        pval = getpval(p);
        % get negative log-likelihood
        if fitrev || fitrep
            % fit reversal and/or repetition curves
            nll = -getll_c(pval{:});
        else
            % fit responses
            nll = -getll_r(pval{:});
        end
    end



    function [ll,ll_sd] = getll_c(varargin)
        % simulate reversal and repetition curves
        c = sim_c(varargin{:});
        % compute log-likelihood
        ll = 0;
        if fitrev % fit reversal curve
            rrev_hat = mean(c.rrev,2);
            rrev_hat = (1-lmin*2)*rrev_hat+lmin;
            ll = ll+ ...
                sum((csub.rrev.*csub.nrev).*log(rrev_hat))+ ...
                sum(((1-csub.rrev).*csub.nrev).*log(1-rrev_hat));
        end
        if fitrep % fit repetition curve
            rrep_hat = mean(c.rrep,2);
            rrep_hat = (1-lmin*2)*rrep_hat+lmin;
            ll = ll+ ...
                sum((csub.rrep.*csub.nrep).*log(rrep_hat))+ ...
                sum(((1-csub.rrep).*csub.nrep).*log(1-rrep_hat));
        end
        if fitcnf % fit confidence
            if fitrev % fit reversal curve
                crev_hat = mean(c.crev,2);
                crev_hat = (1-lmin*2)*crev_hat+lmin;
                ll = ll+ ...
                    sum((csub.crev.*csub.nrev).*log(crev_hat))+ ...
                    sum(((1-csub.crev).*csub.nrev).*log(1-crev_hat));
            end
            if fitrep % fit repetition curve
                crep_hat = mean(c.crep,2);
                crep_hat = (1-lmin*2)*crep_hat+lmin;
                ll = ll+ ...
                    sum((csub.crep.*csub.nrep).*log(crep_hat))+ ...
                    sum(((1-csub.crep).*csub.nrep).*log(1-crep_hat));
            end
        end
        if nargout > 1
            % compute bootstrapped log-likelihood s.d.
            ll_res = zeros(nres,1);
            for ires = 1:nres
                jres = randsample(nsmp,nsmp,true);
                if fitrev % fit reversal curve
                    rrev_hat = mean(c.rrev(:,jres),2);
                    rrev_hat = (1-lmin*2)*rrev_hat+lmin;
                    ll_res(ires) = ll_res(ires)+ ...
                        sum((csub.rrev.*csub.nrev).*log(rrev_hat))+ ...
                        sum(((1-csub.rrev).*csub.nrev).*log(1-rrev_hat));
                end
                if fitrep % fit repetition curve
                    rrep_hat = mean(c.rrep(:,jres),2);
                    rrep_hat = (1-lmin*2)*rrep_hat+lmin;
                    ll_res(ires) = ll_res(ires)+ ...
                        sum((csub.rrep.*csub.nrep).*log(rrep_hat))+ ...
                        sum(((1-csub.rrep).*csub.nrep).*log(1-rrep_hat));
                end
                if fitcnf % fit confidence
                    if fitrev % fit reversal curve
                        crev_hat = mean(c.crev(:,jres),2);
                        crev_hat = (1-lmin*2)*crev_hat+lmin;
                        ll_res(ires) = ll_res(ires)+ ...
                            sum((csub.crev.*csub.nrev).*log(crev_hat))+ ...
                            sum(((1-csub.crev).*csub.nrev).*log(1-crev_hat));
                    end
                    if fitrep % fit repetition curve
                        crep_hat = mean(c.crep(:,jres),2);
                        crep_hat = (1-lmin*2)*crep_hat+lmin;
                        ll_res(ires) = ll_res(ires)+ ...
                            sum((csub.crep.*csub.nrep).*log(crep_hat))+ ...
                            sum(((1-csub.crep).*csub.nrep).*log(1-crep_hat));
                    end
                end
            end
            ll_sd = max(std(ll_res),1e-6);
        end
    end

    function [ll,ll_sd] = getll_r(varargin)
        % simulate responses
        [rbef,raft,cbef,caft] = sim_r(varargin{:});
        % compute log-likelihood
        ll = 0;
        r_hat = mean(raft == 1,2);
        r_hat = (1-lmin*2)*r_hat+lmin;
        r_hat(cfg.raft == 2) = 1-r_hat(cfg.raft == 2);
        ll = ll+sum(log(r_hat));
        if fitcnf
            c_hat = mean(caft == 1,2);
            c_hat = (1-lmin*2)*c_hat+lmin;
            c_hat(cfg.caft == 2) = 1-c_hat(cfg.caft == 2);
            ll = ll+sum(log(c_hat));
        end
        if nargout > 1
            % compute bootstrapped log-likelihood s.d.
            ll_res = zeros(nres,1);
            for ires = 1:nres
                jres = randsample(nsmp,nsmp,true);
                r_hat = mean(raft(:,jres) == 1,2);
                r_hat = (1-lmin*2)*r_hat+lmin;
                r_hat(cfg.raft == 2) = 1-r_hat(cfg.raft == 2);
                ll_res(ires) = ll_res(ires)+sum(log(r_hat));
                if fitcnf
                    c_hat = mean(caft(:,jres) == 1,2);
                    c_hat = (1-lmin*2)*c_hat+lmin;
                    c_hat(cfg.caft == 2) = 1-c_hat(cfg.caft == 2);
                    ll_res(ires) = ll_res(ires)+sum(log(c_hat));
                end
            end
            ll_sd = max(std(ll_res),1e-6);
        end
    end

    function [c,psmp] = sim_c(varargin)
        % simulate responses
        [rbef,raft,cbef,caft,psmp] = sim_r(varargin{:});
        % get reversal and repetition curves
        c = getc(seqind,seqpos,seqdir,seqllr,rbef,raft,cbef,caft);
    end

    function [rbef,raft,cbef,caft,psmp] = sim_r(h,sinf,ssel,scnf,tcnf,gcnf)
        % simulate responses
        xt = zeros(nseq,nsmp); % log-belief
        rbef = nan(nseq,nsmp); % response before sequence
        raft = nan(nseq,nsmp); % response after sequence
        cbef = nan(nseq,nsmp); % confidence before sequence
        caft = nan(nseq,nsmp); % confidence after sequence
        psmp = ones(nseq,1); % fraction of propagated samples
        for iseq = 1:nseq
            % set state before sequence
            if seqind(iseq) > 1
                % set response before sequence
                rbef(iseq,:) = raft(iseq-1,:);
                if fitcnf
                    % set confidence before sequence
                    cbef(iseq,:) = caft(iseq-1,:);
                end
            else % 1st trial of each block
                % use same response as subject
                rbef(iseq,:) = cfg.rbef(iseq);
                if fitcnf
                    % use same confidence as subject
                    cbef(iseq,:) = cfg.cbef(iseq);
                end
            end
            % update log-belief
            if seqind(iseq) > 1
                xt(iseq,:) = upfun(xt(iseq-1,:),h);
            end
            xt(iseq,:) = normrnd(xt(iseq,:)+seqllr(iseq),sqrt(seqlen(iseq))*sinf);
            % apply selection noise
            xr = normrnd(xt(iseq,:),ssel);
            % compute response
            raft(iseq,:) = 1+(xr < 0);
            if fitcnf
                % compute log-belief in favor of response
                xc = xt(iseq,:).*(3-2*raft(iseq,:));
                % apply confidence gain during switches
                iswi = raft(iseq,:) ~= rbef(iseq,:);
                xc(iswi) = xc(iswi)*gcnf;
                % apply confidence noise
                xc = normrnd(xc,scnf);
                % compute confidence
                caft(iseq,:) = 1+(xc > tcnf);
                if resamp % resample log-belief
                    ires = ...
                        raft(iseq,:) == cfg.raft(iseq) & ...
                        caft(iseq,:) == cfg.caft(iseq);
                    psmp(iseq) = mean(ires);
                    if any(ires)
                        % bootstrap resample
                        xt(iseq,:) = xt(iseq,randsample(find(ires),nsmp,true));
                    else
                        xt(iseq,:) = 0;
                    end
                end
            else
                if resamp % resample log-belief
                    ires = raft(iseq,:) == cfg.raft(iseq);
                    psmp(iseq) = mean(ires);
                    if any(ires)
                        % bootstrap resample
                        xt(iseq,:) = xt(iseq,randsample(find(ires),nsmp,true));
                    else
                        xt(iseq,:) = 0;
                    end
                end
            end
        end
    end

    function [o] = rmmissp(o)
        % remove missing parameters from output structure
        for k = 1:npar
            if isfield(o,pnam{k}) && isnan(o.(pnam{k}))
                o = rmfield(o,pnam{k});
            end
        end
    end

end

function [x] = upfun(x,h)
% update log-belief
x = x+log((1-h)./h+exp(-x))-log((1-h)./h+exp(+x));
end