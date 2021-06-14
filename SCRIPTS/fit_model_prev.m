function [out] = fit_model_prev(cfg)

seqind = cfg.seqind(:);
seqdir = cfg.seqdir(:);
epinum = cfg.epinum(:);
seqpos = cfg.seqpos(:);
raft   = cfg.raft(:);

nseq = numel(seqind);
indx = sub2ind([nseq,2],(1:nseq)',raft);

pfit_ini = [ 1;0.75];
pfit_min = [ 0;0.50];
pfit_max = [10;1.00];

pval = fmincon(@fmin,pfit_ini,[],[],[],[],pfit_min,pfit_max,[], ...
    optimset('Display','notify','FunValCheck','on','Algorithm','interior-point','TolX',1e-20,'MaxFunEvals',1e6));


out      = [];
out.tcst = pval(1);
out.pmax = pval(2);
out.pr   = get_pr(out.tcst,out.pmax);

    function [f] = fmin(p)
    f = -get_llh(p(1),p(2));
    end

    function [llh] = get_llh(tcst,pmax)
    pr = get_pr(tcst,pmax);
    llh = sum(log(max(pr(indx),1e-20)));
    end

    function [pr] = get_pr(tcst,pmax)
    pr = nan(nseq,1);
    for iseq = 1:nseq
        if seqind(iseq) == 1
            pini = 0.5;
            if seqdir(iseq) == 1
                pend = pmax;
            else
                pend = 1-pmax;
            end
        elseif seqpos(iseq) == 1
            pini = pcur;
            if seqdir(iseq) == 1
                pend = pmax;
            else
                pend = 1-pmax;
            end
        end
        pcur = pini+(pend-pini)*(1-exp(-seqpos(iseq)/tcst));
        pr(iseq) = pcur;
    end
    pr = [pr,1-pr];
    end

end