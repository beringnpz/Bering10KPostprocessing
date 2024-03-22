function addjul1btemp(sim)
%ADDJUL1BTEMP Add July 1 model-based bottom temp to cold pool index file
%
% addjul1btemp(sim)
% addjul1btemp(sim, p1, v1, ...)
%
% This function calculates cold pool indices using data from a ROMS
% simulation extracted on Jul 1 of each year of the simualtion.  Indices
% are calculated across a few variants of the masking polygons for the
% eastern Bering Sea shelf.  Data is saved to a preexisting file for the
% simulation (see buildcoldpoolncfile.m for details).
%
% Input variables:
%
%   sim:    Name of simulation in the post-processed dataset.
%
% Optional input variables (passed as parameter/value pairs) [default]:
%
%   years:      vector of years for which indices will be extracted (if not
%               already in the file).  Any years not currently in the file
%               will be appended to the time dimension; years should be
%               chronological though this is not enforced by this script.
%               [1970:present] 
%
%   grdfile:    path to ROMS grid file associated with the simulation
%               ['[moxdir]/roms_for_public/Bering10K_extended_grid.nc']

% Copyright 2024 Kelly Kearney
 
% Parse input

p = inputParser;
p.addParameter('years', 1970:year(datetime('today')), @(x) validateattributes({'numeric'}, {'vector', 'integer', '>=', 1970, '<=', year(datetime('today'))}));
p.addParameter('grdfile', fullfile(moxdir, '/roms_for_public/Bering10K_extended_grid.nc'));

p.parse(varargin{:});
Opt = p.Results;

% File setup

outname  = fullfile(moxdir, 'roms_for_public', sim, 'Level3', sprintf('%s_coldpool.nc', sim));

if ~exist(outname, 'file')
    error('Cold pool index file (%s) not found', outname);
end

% Masks

Grd = ncstruct(Opt.grdfile);

M = cpindexmasks(Grd);

nmask = size(M.mask,3);

% Grid weighting

gmask = any(M.mask,3);

extractblock = @(x) x(any(gmask,2),any(gmask,1));

w = extractblock(Grd.area_feast);

% Thresholds

thresh = ncread(outname, 'threshold');
nthresh = length(thresh);

% Read data and calculate indices

ttarget = datetime(Opt.years, 7, 1);

tcomplete = ncdateread(outname, 'time');
tf = ismember(ttarget, tcomplete);
ttarget = ttarget(~tf);


if ~isempty(ttarget)
    fprintf('Calculating July1 indices\n');
    [btemp, ~, mask] = readbtemp('sim', sim, 'dates', ttarget, 'mask', gmask);

    % Calculate indices

    Idx = btemp2indices(btemp, mask, 'w', w, 'thresh', thresh);

    % Add data to file

    fprintf('Adding July 1 indices...\n');    
    % Note: This will add new years after existing... mostly assuming that
    % we will stay in order, though that's not enforced.
    
    sidx = length(tcomplete)+1;
    
    tunit = 'days since 1900-01-01';
    tnc = cftime(ttarget, tunit, [], 'reverse');
    
    ncwrite(outname, 'time', tnc, sidx);
    ncwrite(outname, 'average_bottom_temp', permute(Idx.btemp, [2 3 1]), [1 1 sidx]);
    ncwrite(outname, 'cold_pool_index',     permute(Idx.cpool, [3 2 4 1]), [1 1 1 sidx]);
    
    if length(ttarget) > 1
        ncaddhis(outname, sprintf('Model July 1 data added for years %d-%d', year([min(ttarget) max(ttarget)])));
    else
        ncaddhis(outname, sprintf('Model July 1 data added for year %d', year(ttarget)));
    end

end


