function [c] = getc(seqind,seqpos,seqdir,seqllr,rbef,raft,cbef,caft)

if nargin < 6
    error('Missing input arguments!');
end

% has confidence?
hasc = nargin == 8 && ~isempty(cbef) && ~isempty(caft);

% set trial offsets for reversal curves
ioff = -4:+3;
npos = numel(ioff);

% set bin limits for repetition curves
blim = cat(1, ...
    [-inf,  -3,  -2,  -1,   0,  +1,  +2,  +3], ...
    [  -3,  -2,  -1,   0,  +1,  +2,  +3,+inf]);
nbin = size(blim,2);

% initialize output array
c = [];

% get number of datasets
ndat = size(raft,2);

% get reversal curves
ifilt = find(seqind > 1 & seqpos == 1);
nrev = nan(npos,ndat);
rrev = nan(npos,ndat);
for idat = 1:ndat
    for ipos = 1:npos
        nrev(ipos,idat) = numel(ifilt);
        rrev(ipos,idat) = mean(raft(ifilt+ioff(ipos),idat) == seqdir(ifilt));
    end
end
c.nrev = nrev;
c.rrev = rrev;
if hasc % has confidence
    crev = nan(npos,ndat);
    for idat = 1:ndat
        for ipos = 1:npos
            crev(ipos,idat) = mean(caft(ifilt+ioff(ipos),idat) == 2);
        end
    end
    c.crev = crev;
end

% get repetition curves
ifilt = find(seqind > 1);
nrep = nan(nbin,ndat);
rrep = nan(nbin,ndat);
for idat = 1:ndat
    xi = seqllr(ifilt).*(3-2*rbef(ifilt,idat));
    ri = raft(ifilt,idat) == rbef(ifilt,idat);
    for ibin = 1:nbin
        jfilt = xi >= blim(1,ibin) & xi < blim(2,ibin);
        nrep(ibin,idat) = nnz(jfilt);
        rrep(ibin,idat) = mean(ri(jfilt));
    end
end
c.nrep = nrep;
c.rrep = rrep;
if hasc % has confidence
    crep = nan(nbin,ndat);
    for idat = 1:ndat
        xi = seqllr(ifilt).*(3-2*rbef(ifilt,idat));
        ci = caft(ifilt,idat) == 2;
        for ibin = 1:nbin
            jfilt = xi >= blim(1,ibin) & xi < blim(2,ibin);
            crep(ibin,idat) = mean(ci(jfilt));
        end
    end
    c.crep = crep;
    
    % fraction confident for within subject error bars
    c.pconf = mean(caft == 2);

end

end