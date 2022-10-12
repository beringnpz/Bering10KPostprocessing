function Idx = btemp2indices(btemp, gmask, varargin)
%BTEMP2INDICES Calculate mean temp and frac-below indices
%
% btemp2indices(btemp, gmask)
% btemp2indices(btemp, gmask, w, )
% btemp2indices(btemp, gmask, w, thresh)
%
% Input variables:
%
%   btemp:  bottom temperature array, xi x eta x time
%
%   gmask:  mask, xi x eta, calculate indices for gmask=1 cells only
%
% Optional input variables:
%
%   w:      grid weighting [ones(size(gmask)]
%
%   thresh: temperature threshholds for cold pool indices.  [0 1 2]
%
% Output variables:
%
%   Idx:    table with columns:
%           t:      datetime
%           btemp:  mean bottom temp in masked region
%           cpoolx: fraction of region with temp less than x

[nxi, neta] = size(btemp);

p = inputParser;
p.addParameter('w', true(nxi,neta));
p.addParameter('thresh', [0 1 2]);

p.parse(varargin{:});
Opt = p.Results;

% Set up # masks and thresholds

if ~isequal([size(btemp,1) size(btemp,2)], [size(gmask,1) size(gmask,2)])
    error('First two dimensions of btemp and gmask must match');
end

nmask = size(gmask,3);
nt = size(btemp,3);
nthresh = length(Opt.thresh);

% Calculate indices

Idx.thresh = Opt.thresh;
Idx.btemp = nan(nt,nmask);
Idx.cpool = nan(nt,nmask,nthresh);


for im = 1:nmask
    Idx.btemp(:,im) = local(btemp, gmask(:,:,im), 'weight', Opt.w, 'omitnan');

    for it = 1:nthresh
        Idx.cpool(:,im,it) = local(btemp<Opt.thresh(it), gmask(:,:,im), 'weight', Opt.w, 'omitnan');
    end
end

% % Threshhold column names
% 
% cpcol = cellstr(num2str(Opt.thresh(:)));
% cpcol = strtrim(strrep(cpcol, '.', 'p'));
% 
% Idx = struct('t', t);
% Idx.btemp = local(btemp, gmask, 'weight', Opt.w, 'omitnan');
% for ii = 1:length(Opt.thresh)
%     Idx.(['cpool' cpcol{ii}]) = local(btemp<Opt.thresh(ii), gmask, 'weight', Opt.w, 'omitnan');
% end
% 
% Idx = struct2table(Idx);