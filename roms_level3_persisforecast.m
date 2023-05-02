function roms_level3_persisforecast(varargin)
%ROMS_LEVEL3_PERSISFORECAST Calculate bottom temp. persistence forecast
%
% roms_level3_persisforecast
% roms_level3_persisforecast(param, val, ...
%
% This function calculates a persistence forecast for bottom temperature
% extending from the last simulation date of the indicated year through the
% end of summer of that year.  It is intended to be applied to a Bering Sea
% hindcast simulation after each spring (late March through early May)
% update to forecast potential conditions in the coming summer survey.
% Output includes a climatological weekly timeseries of bottom temperature
% (based on all years prior to the indicated year), an anomaly dataset for
% the last time step in the simulation relative to that climatology, and
% a forecast timeseries extending from the last hindcast time step through
% Sep 30.
%
% Output is saved to the Level 3 folder of the specified simulation in a
% file named <simname>_<year>forecast_temp_bottom5m.nc.  If a file with
% this name already exists, no action will be taken.  
%
% Optional input variables (passed as parameter/value pairs)
%
%   simname:    string, name of simulation (in roms_for_public folder).
%               This function assumes that Level 2 bottom temperature files
%               have already been created for the indicated simulation.
%               ['B10K-K20_CORECFS']
%
%   grdfile:    name of grid file for simulation.  This file should include
%               the extended strata mask variables. The function will
%               first check to see if the file name resolves as a full or
%               relative path name, then will look for a file matching the
%               name in the primary roms_for_public folder
%               ['Bering10K_extended_grid.nc']
%
%   year:       scalar indicating the last year of the hindcast.  It is
%               expected that the specified simulation ends midway through
%               the specified year.  Alternatively, you can enter a
%               specific datetime; the forecast will be based on data up to
%               that date.
%               [current year]

% Copyright 2021 Kelly Kearney

%--------------------
% Parse input
%--------------------

p = inputParser;
p.addParameter('sim', 'B10K-K20_CORECFS',  @(x) validateattributes(x, {'char'}, {'scalartext'}));
p.addParameter('grdfile', 'Bering10K_extended_grid.nc',  @(x) validateattributes(x, {'char'}, {'scalartext'}));
p.addParameter('year', year(datetime('today')));

p.parse(varargin{:});
Opt = p.Results;

if exist(Opt.grdfile, 'file')
    grdfile = Opt.grdfile;
else
    grdfile = fullfile(moxdir, 'roms_for_public', Opt.grdfile);
end
if ~exist(grdfile, 'file')
    error('Could not find indicated grid file: %s\n', Opt.grdfile);
end

yrisint = isnumeric(Opt.year) && floor(Opt.year)==Opt.year && Opt.year >= 1970;
yrisdt  = isa(Opt.year, 'datetime') && Opt.year >= datetime(1970,1,15);

if ~(yrisint || yrisdt)
    error('Year input must be an interger >= 1970 or a datetime between 1970/01/15 and today');
end

if yrisdt
    cutoffdate = Opt.year;
    Opt.year = year(Opt.year);
end
    
% Setup

Grd = ncstruct(grdfile);
[nxi, neta] = size(Grd.h);

hcfol = fullfile(moxdir, 'roms_for_public', Opt.sim);

fcfile = fullfile(hcfol, 'Level3', sprintf('%s_%dforecast_temp_bottom5m.nc', Opt.sim, Opt.year));
if exist(fcfile, 'file')
    error('Output file %s already exists; exiting', fcfile);
end

%--------------------
% Calculate forecast
%--------------------

F = dir(fullfile(hcfol, 'Level2', '*temp_bottom5m.nc'));
fname = fullfile({F.folder},{F.name});

% Read bottom temp data

Btemp = ncstruct(fname, 'temp');
Btemp.t = ncdateread(fname, 'ocean_time');
nt = length(Btemp.t);

% Calculate weekly climatology

doybin = linspace(0,1,53);
idx = discretize(doy(Btemp.t, 'remdec'), doybin);

isclim = year(Btemp.t) < Opt.year;
if yrisdt
    iscurrent = Btemp.t <= cutoffdate;
else
    iscurrent = year(Btemp.t) == Opt.year;
end

Data.btemp_clim = splitapply(@(x) nanmean(x,1), ...
    permute(Btemp.temp(:,:,isclim),[3 1 2]), idx(isclim));
Data.btemp_clim = permute(Data.btemp_clim, [2 3 1]);

% Calculate anomaly of last time step

lastidx = find(iscurrent, 1, 'last');
Data.btemp_anom = Btemp.temp(:,:,lastidx) - Data.btemp_clim(:,:,idx(lastidx));
Data.tanom = Btemp.t(lastidx) + days([-3.5 3.5]);
Data.tanom.Format = 'uuuu/MM/dd HH:mm';

isleap = @(x) (mod(x,4)==0 & mod(x,100)~=0) | mod(x,400) == 0;
if isleap(Opt.year)
    Data.clim_time_bnds = datetime(Opt.year,1,1) + days(doybin)*366;
else
    Data.clim_time_bnds = datetime(Opt.year,1,1) + days(doybin)*365;
end
Data.clim_time = Data.clim_time_bnds(1:end-1)+diff(Data.clim_time_bnds)./2;

Data.clim_years = minmax(year(Btemp.t(isclim)));

% Persistence forecast

Data.forecast = Data.btemp_clim + Data.btemp_anom;

Data.forecast(:,:,Data.clim_time < Data.tanom(1)) = NaN;
Data.forecast(:,:,month(Data.clim_time) > 9) = NaN;

%--------------------
% Create file
%--------------------

tunit = 'days since 1900-01-01';

F = ncschema_init('classic');

hisstr = sprintf('%s: Persistence forecast %d data calculated via roms_level3_persisforecast.m', ...
    datestr(now, 'ddd mmm dd HH:MM:SS yyyy'), ...
    Opt.year);

F = ncschema_addatts(F, ...
    'Simulation', Opt.sim, ...
    'history', hisstr);

F = ncschema_adddims(F, ...
    'xi_rho',    nxi,     false, ...
    'eta_rho',   neta,    false, ...
    'time',      52,      true, ...
    'bnd',       2,       false);

% Coordinate variables

F = ncschema_addvars(F, ...
    'time', ...
    {'time'}, ...
    {'climatology', 'climatology_bnds', ...
    'units', tunit, ...
    'calendar', 'standard'}, ...
    'double');

% Auxiliary coordinate variables

F = ncschema_addvars(F, ...
    'climatology_bnds', ...
    {'bnd', 'time'}, ...
    {'units', tunit, ...
    'calendar', 'standard'}, ...
    'double');
F = ncschema_addvars(F, ...
    'lat_rho', ...
    {'xi_rho', 'eta_rho'}, ...
    {'long_name', 'latitude of RHO-points', ...
    'units', 'degree_north'}, ...
    'double');
F = ncschema_addvars(F, ...
    'lon_rho', ...
    {'xi_rho', 'eta_rho'}, ...
    {'long_name', 'longitude of RHO-points', ...
    'units', 'degree_east'}, ...
    'double');

% Main data variables

F = ncschema_addvars(F, ...
    'clima_temp_bottom5m', ...
    {'xi_rho', 'eta_rho', 'time'}, ...
    {'long_name', 'climatological weekly bottom temperature', ...
    'units', 'Celsius', ...
    'cell_methods', sprintf('time: mean within weeks, time: mean over years %d-%d', Data.clim_years), ...
    'coordinates', 'lat_rho lon_rho time'}, ...
    'double');

F = ncschema_addvars(F, ...
    'anom_temp_bottom5m', ...
    {'xi_rho', 'eta_rho'}, ...
    {'long_name', sprintf('bottom temperature weekly anomaly for %s -- %s', Data.tanom), ...
    'units', 'Celsius', ...
    'coordinates', 'lat_rho lon_rho'}, ...
    'double');

F = ncschema_addvars(F, ...
    'forecast_temp_bottom5m', ...
    {'xi_rho', 'eta_rho', 'time'}, ...
    {'long_name', 'bottom temperature persistence forecast (climatology + anomaly))', ...
    'units', 'Celsius', ...
    'coordinates', 'lat_rho lon_rho'}, ...
    'double');


% Create file and add data

ncwriteschema(fcfile, F);

ncwrite(fcfile, 'time',               cftime(Data.clim_time, tunit, [], 'reverse'));
ncwrite(fcfile, 'climatology_bnds',   cftime([Data.clim_time_bnds(1:end-1); Data.clim_time_bnds(2:end)], tunit, [], 'reverse'));
ncwrite(fcfile, 'lat_rho',            Grd.lat_rho);
ncwrite(fcfile, 'lon_rho',            Grd.lon_rho);
ncwrite(fcfile, 'clima_temp_bottom5m', Data.btemp_clim);
ncwrite(fcfile, 'anom_temp_bottom5m',  Data.btemp_anom);
ncwrite(fcfile, 'forecast_temp_bottom5m', Data.forecast);
