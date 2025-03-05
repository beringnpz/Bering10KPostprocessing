function buildcoldpoolncfile(sim)
%BUILDCOLDPOOLNCFILE Create cold pool index file
%
% buildcoldpoolncfile(sim)
% buildcoldpoolncfile(sim, 'years', yrs)
%
% This function builds the cold pool index file associated with a ROMS
% hindcast simulation.  This is part of the Level 3 postprocessing.  The
% file includes two primary variables: average_bottom_temperature and
% cold_pool_index, with several variants for different polygon masks and
% for different methods of calculation.  See file metadata for details.
%
% Only metadata is added to the file on creation.  To add data, see
% addsurveyreplicatedbtemp.m and addjul1btemp.m in this folder, and
% analysis/cold_pool_for_regional_models.Rmd in the
% afsc-gap-products/coldpool package.
%
% Input variables:
%
%   sim:    Name of simulation in the post-processed dataset.  This is
%           expected to correspond to a hindcast simulation covering the
%           Bering Sea shelf.

% Copyright 2024 Kelly Kearney

% p = inputParser;
% p.addParameter('years', 1970:year(datetime('today')), @(x) validateattributes({'numeric'}, {'vector', 'integer', '>=', 1970, '<=', year(datetime('today'))}));
% 
% p.parse(varargin{:});
% Opt = p.Results;

%--------------------
% Setup
%--------------------

% File name

outname  = fullfile(moxdir, 'roms_for_public', sim, 'Level3', sprintf('%s_coldpool.nc', sim));

if exist(outname, 'file')
    warning('Cold pool file for %s already exists; no action taken', sim);
    return
end

% Thresholds

thresh = [-1 0 1 1.5 2]; % deg C
nthresh = length(thresh);

% Masks

M = cpindexmasks;

%--------------------
% Build file schema
%--------------------

A = struct;

tunit = 'days since 1900-01-01';
% A.time = cftime(ttarget, tunit, [], 'reverse');

A.thresh = thresh;

A.method = char({'Model July', 'Survey', 'Model survey-replicated'});
A.region = char(M.name);

hisstr = sprintf('%s: %s', ...
    datetime('now', 'format', 'eee MMM dd HH:mm:ss yyyy'), ...
    'Cold pool indices file created via buildcoldpoolncfile.m');

% NetCDF file schema setup

F = ncschema_init('classic');

F = ncschema_addatts(F, ...
    'Simulation', sim, ...
    'history', hisstr);

F = ncschema_adddims(F, ...
    'region',    nmask,            false, ...
    'method',    size(A.method,1), false, ...
    'threshold', nthresh,          false, ...
    'len1',      size(A.region,2), false, ...
    'len2',      size(A.method,2), false, ...
    'time',      length(A.time),   true);

% Coordinate variables

F = ncschema_addvars(F, ...
    'time', ...
    {'time'}, ...
    {'long_name', 'time', ...
     'units', tunit, ...
     'calendar', 'standard'}, ...
    'double');
F = ncschema_addvars(F, ...
    'threshold', ...
    {'threshold'}, ...
    {'long_name', 'cold pool index threshold', ...
     'units', 'Celsius'}, ...
    'double');

% Auxiliary coordinate variables

F = ncschema_addvars(F, ...
    'region_label', ...
    {'len1', 'region', }, ...
    {'long_name', 'geographic region'}, ...
    'char');
F = ncschema_addvars(F, ...
    'method_label', ...
    {'len2', 'method'}, ...
    {'long_name', 'index calculation source/method'}, ...
    'char');

% Main data variables

F = ncschema_addvars(F, ...
    'average_bottom_temp', ...
    {'region', 'method', 'time'}, ...
    {'long_name', 'average bottom temperature in region', ...
     'units', 'Celsius', ...
     'coordinates', 'region_label, method_label'}, ...
    'double');
F = ncschema_addvars(F, ...
    'cold_pool_index', ...
    {'threshold', 'region', 'method', 'time'}, ...
    {'long_name', 'fraction of region less than threshold', ...
     'coordinates', 'region_label, method_label'}, ...
    'double');

%--------------------
% Create file and 
% add data
%--------------------

fprintf('Creating file...\n');

ncwriteschema(outname, F);

fprintf('Adding dimension data...\n');

% ncwrite(outname, 'time',                A.time);
ncwrite(outname, 'threshold',           A.thresh);
ncwrite(outname, 'region_label',        A.region');
ncwrite(outname, 'method_label',        A.method');    
% ncwrite(outname, 'average_bottom_temp', permute(Idx.btemp, [2 3 1])); % Note: don't need to specify start b/c all = 1
% ncwrite(outname, 'cold_pool_index',     permute(Idx.cpool, [3 2 4 1]));

% if length(ttarget) > 1
%     ncaddhis(outname, sprintf('Model July 1 data added for years %d-%d', year(minmax(ttarget))));
% else
%     ncaddhis(outname, sprintf('Model July 1 data added for year %d', year(ttarget)));
% end

% Add Fill Value attribute to primary variables, which will include
% missing values (in years where survey and model do not overlap)

ncid = netcdf.open(outname);
vid = netcdf.inqVarID(ncid, 'average_bottom_temp');
[~,fillValue] = netcdf.inqVarFill(ncid,vid);
netcdf.close(ncid);

ncwriteatt(outname, 'average_bottom_temp', '_FillValue', fillValue);
ncwriteatt(outname, 'cold_pool_index', '_FillValue', fillValue);