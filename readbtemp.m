function [btemp, tmod, mask] = readbtemp(varargin)
%READJUL1BTEMP Reads bottom temperature from simulation
%
% readbtemp
% readbtemp(param1, val1, ...)
%
% Optional input arguments (passed as parameter/value pairs):
%
%   sim:    name of simulation, expected to be found under roms-for-public
%           folder ['B10K-K20_CORECFS']
%
%   dates:  dates to read.  Can be either:
%           NaT:  all model times read 
%           2-element datetime array: read all times between dates
%           n-element datetime array: read specific times closest to each
%           [datetime(1970:present, 07, 01)]
%
%   mask:   mask determining horizontal hyperslab to read.  Can be 2D or
%           3D, with first 2 dimensions matching the size of the ROMS grid.  
%           Data will be read within the bounding box of any(mask,3).  If
%           empty, entire domain is read. [[]]
%
% Output:
%
%   btemp:  nxi_trim x neta_trim x nyr array of bottom temp values
%
%   t:      model dates matching requested.
%
%   mask:   copy of mask, trimmed to the bounding box of any(mask,3) 

% Copyright 2022 Kelly Kearney


% Parse input

p = inputParser;

p.addParameter('sim', 'B10K-K20_CORECFS',  @(x) validateattributes(x, {'char'}, {'scalartext'}));
yrnow = year(datetime('today'));
if datetime('today') > datetime(yrnow,7,1)
    yrmax = yrnow;
else
    yrmax = yrnow - 1;
end
p.addParameter('dates', datetime(1970:yrmax,7,1),  @(x) validateattributes(x, {'datetime'}, {}));
p.addParameter('mask', []);

p.parse(varargin{:});
Opt = p.Results;


% Bottom temperature files associated with simulation

F = dir(fullfile(moxdir, 'roms_for_public', Opt.sim, 'Level2', '*average_temp_bottom5m.nc'));
fname = fullfile({F.folder}, {F.name});

% Check mask size

if ~isempty(Opt.mask)
    I = ncinfo(fname{1});
    [~,loc] = ismember({'xi_rho', 'eta_rho'}, {I.Dimensions.Name});
    nxi = I.Dimensions(loc(1)).Length;
    neta = I.Dimensions(loc(2)).Length;
    validateattributes(Opt.mask, {'logical'}, {'size', [nxi neta NaN]});
end

% Read times from relevant files

if ~isempty(Opt.dates)
    tlim = ncdatelim(fname, 'ocean_time');
    isin = ~(tlim(:,2)<min(Opt.dates) | tlim(:,1)>max(Opt.dates));
    fname = fname(isin);
end

% I = cellfun(@(x) ncinfo(x, 'ocean_time'), fname);
% fidx = arrayfun(@(a,b) ones(1,a)*b, [I.Size], 1:length(fname), 'uni', 0);
% tidx = arrayfun(@(a) 1:a, [I.Size], 'uni', 0);
% fidx = cat(2, fidx{:});
% tidx = cat(2, tidx{:});
t = ncdateread(fname, 'ocean_time');

% Set up time hyperslabs to minimize number of reads

if isnat(Opt.dates)
    tmod = t;
    tscs = [1 Inf 1];
elseif length(Opt.dates) == 2
    isin = t >= Opt.dates(1) & t <= Opt.dates(2);
    tmod = t(isin);
    tscs = [find(isin,1) sum(isin) 1];
else
    islice = interp1(t, 1:length(t), Opt.dates, 'nearest');
    tmod = t(islice);
    
    tmp = [1 diff(islice)] ~= 1;
    
    sidx = islice(union(1, find(tmp)));
    eidx = islice(union(find(tmp)-1, length(islice)));
    
    tscs = [sidx; eidx-sidx+1; ones(size(sidx))]';
end

% If mask provided, determine horizontal hyperslab

if isempty(Opt.mask)
    hscs = {};
else
    [xi,eta] = ind2sub([nxi neta], find(any(Opt.mask, 3)));
    xilim = minmax(xi);
    etalim = minmax(eta);
    mask = Opt.mask(xilim(1):xilim(2), etalim(1):etalim(2), :);
    hscs = {'xi_rho', [xilim(1) diff(xilim)+1 1], 'eta_rho', [etalim(1) diff(etalim)+1 1]};
end
    
% Read data

nread = size(tscs,1);
btemp = cell(nread,1);

for ii = 1:nread
    Tmp = ncstruct(fname, 'temp', struct('ocean_time', tscs(ii,:), hscs{:}));
    btemp{ii} = Tmp.temp;
end
btemp = cat(3, btemp{:});
