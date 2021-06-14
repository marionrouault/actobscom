function [out] = fit_model_crev(cfg)

seqind = cfg.seqind(:);
seqdir = cfg.seqdir(:);
epinum = cfg.epinum(:);
seqpos = cfg.seqpos(:);
caft   = cfg.caft(:);

nseq = numel(seqind);
indx = sub2ind([nseq,2],(1:nseq)',3-caft);

pfit_ini = [ 1;0.75;0.25];
pfit_min = [ 0;0.00;0.00];
pfit_max = [10;1.00;1.00];

pval = fmincon(@fmin,pfit_ini,[],[],[],[],pfit_min,pfit_max,[], ...
    optimset('Display','notify','FunValCheck','on','Algorithm','interior-point','TolX',1e-20,'MaxFunEvals',1e6, ...
    'MaxIter',1e4));

out      = [];
out.tcst = pval(1);
out.pmax = pval(2);
out.pmin = pval(3);
out.pr   = get_pr(out.tcst,out.pmax,out.pmin);

% compute confidence drop as extra variable
out.pdrp = out.pmin+(out.pmax-out.pmin)*(1-exp(-1/out.tcst));
out.pdrp = out.pmax-out.pdrp;

    function [f] = fmin(p)
    f = -get_llh(p(1),p(2),p(3));
    end

    function [llh] = get_llh(tcst,pmax,pmin)
    pr = get_pr(tcst,pmax,pmin);
    llh = sum(log(max(pr(indx),1e-20)));
    end

    function [pr] = get_pr(tcst,pmax,pmin)
    pr = nan(nseq,1);
    for iseq = 1:nseq
        if seqpos(iseq) == 1
            pini = pmin;
            pend = pmax;
        end
        pcur = pini+(pend-pini)*(1-exp(-seqpos(iseq)/tcst));
        pr(iseq) = pcur;
    end
    pr = [pr,1-pr];
    end

end