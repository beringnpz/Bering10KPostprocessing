function roms_level3_aclimindices(varargin)
%ROMS_LEVEL3_ACLIMINDICES Calculate ACLIM-relevant 1D index timeseries
%
% roms_level3_aclimindices(simname)
% roms_level3_aclimindices(simname, param, val, ...)
%
% This function calculates a variety of one-dimensional index timeseries
% for many of the output variables for a given simulation.  The specific
% indices were developed as part of the Alaska Climate Integrated Modeling
% (ACLIM) project; this function is designed to be applied to the long-term
% forecasts and hindcasts associated with that project.  
%
% The current set of indices includes spatially-averaged versions of
% either two-dimensional model output variables or depth-integrated,
% surface-averaged, or bottom-averaged three-dimensional output variables.
% The regional averages are calculated in two ways:
%   1) Strata-based regions: Model data is spatially averaged within
%      predefined survey strata polygons.  Output is on the same time axis
%      as the simulation average file output.
%   2) Survey-replicated: Model data is resampled annually in time and
%      space based on the mean location and day-of-sampling from the
%      groundfish survey dataset (see Kearney 2021 for details).  The time
%      axis for this set of indices runs from 1970-2100, with placeholder
%      fill values used for periods outside a simulation's extent. 
%
% Output is saved to the Level 3 folder of the specified simulation in a
% files named ACLIMregion_<simname>.nc and ACLIMsurveyrep_<simname>.nc.  If
% a file with either name already exists, no action will be taken. 
%
% Input variables:
%
%   simname:    string, name of simulation (in roms_for_public folder).
%               This function assumes that Level 1 and 2 bottom
%               post-processed files have already been created for the
%               indicated simulation.
%
% Optional input variables (passed as parameter/value pairs, default in
% []):
%
%   svyfile:    name of Excel file holding survey data. The function will
%               first check to see if the file name resolves as a full or
%               relative path name, then will look for a file matching the
%               name in the primary roms_for_public folder.
%               ['AFSC_groundfish_survey_temperature_1982-2020.xlsx']
%
%   grdfile:    name of grid file for simulation.  This file should include
%               the extended strata mask variables. The function will
%               first check to see if the file name resolves as a full or
%               relative path name, then will look for a file matching the
%               name in the primary roms_for_public folder
%               ['Bering10K_extended_grid.nc']

% Copyright 2021 Kelly Kearney

%--------------------
% Parse input
%--------------------

p = inputParser;
p.addRequired('sim', @(x) validateattributes(x, {'char'}, {'scalartext'}));
p.addOptional('svyfile', 'AFSC_groundfish_survey_temperature_1982-2020.xlsx',  @(x) validateattributes(x, {'char'}, {'scalartext'}));
p.addOptional('grdfile', 'Bering10K_extended_grid.nc',  @(x) validateattributes(x, {'char'}, {'scalartext'}));

p.parse(varargin{:});
Opt = p.Results;

% Look for survey and grid files

if exist(Opt.svyfile, 'file')
    svyfile = Opt.svyfile;
else
    svyfile = fullfile(moxdir, 'roms_for_public', Opt.svyfile);
end
if ~exist(svyfile, 'file')
    error('Could not find indicated survey data file: %s\n', Opt.svyfile);
end

if exist(Opt.grdfile, 'file')
    grdfile = Opt.grdfile;
else
    grdfile = fullfile(moxdir, 'roms_for_public', Opt.grdfile);
end
if ~exist(grdfile, 'file')
    error('Could not find indicated grid file: %s\n', Opt.grdfile);
end

% Check for existing output

filereg = fullfile(moxdir, 'roms_for_public', Opt.sim, 'Level3', ...
    sprintf('ACLIMregion_%s.nc', Opt.sim));

filesrep = fullfile(moxdir, 'roms_for_public', Opt.sim, 'Level3', ...
    sprintf('ACLIMsurveyrep_%s.nc', Opt.sim));

if exist(filereg,'file') || exist(filesrep, 'file')
    error('Output file(s) %s and/or %s already exist; exiting', filereg, filesrep);
end

%--------------------
% Setup
%--------------------

fprintf('Setup...\n');

% ACLIM indices

lev1 = {'Ben', 'DetBen', ...
        'Hsbl', ...
        'IceNH4', 'IceNO3', 'IcePhL', ...
        'aice', 'hice', ...
        'shflux', 'ssflux'};
lev2 = {'Cop_integrated'      , 'Cop_surface5m', ...       
        'EupO_integrated'     , 'EupO_surface5m', ...      
        'EupS_integrated'     , 'EupS_surface5m', ...      
        'Iron_bottom5m'       , 'Iron_integrated',     'Iron_surface5m', ...      
        'Jel_integrated'      , 'Jel_surface5m', ...       
        'MZL_integrated'      , 'MZL_surface5m', ...       
        'NCaO_integrated'     , 'NCaO_surface5m', ...      
        'NCaS_integrated'     , 'NCaS_surface5m', ...      
        'NH4_bottom5m'        , 'NH4_integrated',      'NH4_surface5m', ...       
        'NO3_bottom5m'        , 'NO3_integrated',      'NO3_surface5m', ...       
        'PhL_integrated'      , 'PhL_surface5m' , ...      
        'PhS_integrated'      , 'PhS_surface5m' , ...      
        'prod_Cop_integrated' , ... 
        'prod_EupO_integrated', ... 
        'prod_EupS_integrated', ... 
        'prod_Eup_integrated' , ... 
        'prod_Jel_integrated' , ... 
        'prod_MZL_integrated' , ... 
        'prod_NCaO_integrated', ... 
        'prod_NCaS_integrated', ... 
        'prod_NCa_integrated' , ... 
        'prod_PhL_integrated' , ... 
        'prod_PhS_integrated' , ... 
        'salt_surface5m'      , ... 
        'temp_bottom5m'       , 'temp_integrated',     'temp_surface5m', ...      
        'uEast_bottom5m'      , 'uEast_surface5m'  , ...   
        'vNorth_bottom5m'     , 'vNorth_surface5m'};
    
% Idealized Survey-replicated setup: include all stations on the primary
% 20-n.mi grid (including north but not including corner crab-resample
% stations)

Svy = readtable(svyfile, 'sheet', 'StationSummary');
isin = strcmp(Svy.TYPE, 'main');
Svy = Svy(isin,:);

yr = (1970:2100)';

tsvy = datetime(repmat(yr,1,height(Svy)),1,1) + repmat(days(Svy.DOY'-1),length(yr),1);

xisrep = repmat(Svy.B10K_XI', length(yr), 1);
etasrep = repmat(Svy.B10K_ETA', length(yr), 1);
    
% Grid

Grd = ncstruct(grdfile);
[nxi, neta] = size(Grd.h);

% Survey strata for regional masks

stratafile1 = fullfile(moxdir, 'kearney/simAnalysis/gis_updated', 'EBS_NBS_2019.shp');
stratafile2 = fullfile(moxdir, 'kearney/simAnalysis/ebsshelf_all2', 'ebsshelf_all.shp');
Strata1 = shapeprjread(stratafile1);
Strata2 = shapeprjread(stratafile2);
for ii = 1:length(Strata1)
    Strata1(ii).Lon = wrapTo360(Strata1(ii).Lon);
end
for ii = 1:length(Strata2)
    Strata2(ii).Lon = wrapTo360(Strata2(ii).Lon);
end

% Regional masks: Survey strata.  For the outer shelf bits, we include
% versions both with and without the bit beond the model 200-m isobath.
% (Also cut strata 160 and 170, which is the western shelf... only a small
% sliver of that one is within the model-shelf-break). 

isshelf = Grd.h <= 200;

[sidx, stratanum] = findgroups(Grd.surveystrata_comboeast(:));

gmask = arrayfun(@(x) Grd.surveystrata_comboeast == x, stratanum, 'uni', 0);
gmask_shelf = arrayfun(@(x) Grd.surveystrata_comboeast == x & isshelf, stratanum, 'uni', 0);

issame = cellfun(@(x,y) isequal(x,y), gmask, gmask_shelf);

gmask = [gmask(stratanum>0); gmask_shelf(~issame & stratanum>0)];

regname = [compose('stratum%03d', stratanum(stratanum>0)); ...
           compose('stratum%03db', stratanum(~issame & stratanum>0))]; % text label, b = shelf-cut version
regnum = [stratanum(stratanum>0); -stratanum(~issame & stratanum>0)]; % numeric label, negative = shelf-cut version

gmask = cat(3, gmask{:});
nmask = size(gmask,3);

% Set up horizontal subsetting (to minimize data read)

gidx = find(any(gmask,3));
[ixi,ieta] = ind2sub([nxi neta], gidx);

xilim = minmax([ixi; xisrep(:)]);
etalim = minmax([ieta; etasrep(:)]);

gmask = gmask(xilim(1):xilim(2), etalim(1):etalim(2), :);

gweight = Grd.area_feast(xilim(1):xilim(2), etalim(1):etalim(2));

garea = zeros(nmask,1);
for ii = 1:size(gmask,3)
    garea(ii) = sum(gweight(gmask(:,:,ii)));
end
    
Scs = struct('xi_rho',  [xilim(1)  diff(xilim)+1  1], ...
             'eta_rho', [etalim(1) diff(etalim)+1 1]);
         
xisrep = xisrep - xilim(1) + 1;
etasrep = etasrep - etalim(1) + 1;

% Get list of variable attributes from files

cfile = dir(fullfile(moxdir, 'roms_for_public', Opt.sim, 'Level1', ...
    sprintf('%s_*_average_constants.nc', Opt.sim)));
cfile = fullfile(cfile(1).folder, cfile(1).name);
Ic = ncinfo(cfile, 's_rho');
nlayer = Ic.Size;

for ii = length(lev1):-1:1
    V(ii).short = lev1{ii};
     V(ii).internalshort = V(ii).short;
    
    F = dir(fullfile(moxdir, 'roms_for_public', Opt.sim, 'Level1', sprintf('%s_*_average_%s.nc', Opt.sim, V(ii).short)));
    fname = fullfile(F(1).folder, F(1).name);
    I = ncinfo(fname, V(ii).short);
    [tf,loc] = ismember({'long_name', 'units'}, {I.Attributes.Name});
    if tf(1)
        V(ii).long = I.Attributes(loc(1)).Value;
    else
        V(ii).long = '';
    end
    if tf(2)
        V(ii).unit = I.Attributes(loc(2)).Value;
    else
        V(ii).unit = '';
    end 
    V(ii).level = 1;
end

Vtbl = struct2table(V);
clear V;

missing = false(size(lev2));
for ii = length(lev2):-1:1
    vfname = strrep(lev2{ii}, '_integrated', '');
    vfname = strrep(vfname, '_surface5m', '');
    vfname = strrep(vfname, '_bottom5m', ''); % variable name in file
    
    V(ii).short = lev2{ii};
    V(ii).internalshort = vfname;
    
    F = dir(fullfile(moxdir, 'roms_for_public', Opt.sim, 'Level2', sprintf('%s_*_average_%s.nc', Opt.sim, V(ii).short)));
    if isempty(F)
        missing(ii) = true;
    else
        fname = fullfile(F(1).folder, F(1).name);
        I = ncinfo(fname, V(ii).internalshort);
        [tf,loc] = ismember({'long_name', 'units'}, {I.Attributes.Name});
        if tf(1)
            V(ii).long = I.Attributes(loc(1)).Value;
        else
            V(ii).long = '';
        end
        if tf(2)
            V(ii).unit = I.Attributes(loc(2)).Value;
        else
            V(ii).unit = '';
        end 
        V(ii).level = 2;
    end
end

Vtbl = [Vtbl; struct2table(V(~missing))];
nvar = height(Vtbl); % extras will be treated differently

% Extra variables (not just a spatial average of a Level 1/2)

newdata = {...
    'fracbelow0' '' 'fraction of region with bottom temperature below 0 deg C', '', 0
    'fracbelow1' '' 'fraction of region with bottom temperature below 1 deg C', '', 0
    'fracbelow2' '' 'fraction of region with bottom temperature below 2 deg C', '', 0
    };

Vtbl = [...
    Vtbl
    cell2table(newdata, 'variablenames', Vtbl.Properties.VariableNames)];

% Time analysis

fprintf('Reading times...\n');

avgfiles = dir(fullfile(moxdir, 'roms_for_public', Opt.sim, 'Level1', sprintf('%s_*_average_%s.nc', Opt.sim, lev1{1})));
avgfiles = fullfile({avgfiles.folder}, {avgfiles.name});

isextra = contains(avgfiles, '2020-2020'); % TODO
avgfiles = avgfiles(~isextra);

navg = length(avgfiles);
t = cellfun(@(x) ncdateread(x, 'ocean_time'), avgfiles, 'uni', 0);
fileidx = cellfun(@(a,b) ones(size(a))*b, t, num2cell(1:navg), 'uni', 0);
ftimeidx = cellfun(@(a) (1:length(a))', t, 'uni', 0);
t = cat(1, t{:});
fileidx = cat(1, fileidx{:});
ftimeidx = cat(1, ftimeidx{:});

[t,tidx] = unique(t, 'last');

% Match times to survey-rep grid

tsvyidx = interp1(t, tidx, tsvy, 'nearest');
fsvyidx  = nan(size(tsvyidx)); % which file is each time index in?
ftsvyidx = nan(size(tsvyidx)); % which time step in file corresponds to the time?
isn = isnan(tsvyidx);

fsvyidx(~isn) = fileidx(tsvyidx(~isn));
ftsvyidx(~isn) = ftimeidx(tsvyidx(~isn));

%--------------------
% Create output files
%--------------------

% REGIONAL AVERAGES

F = ncschema_init('classic');

hisstr = sprintf('%s: %s', ...
            datestr(now, 'ddd mmm dd HH:MM:SS yyyy'), ...
            'Regional ACLIM indices file created via roms_level3_aclimindices.m');
        
F = ncschema_addatts(F, ...
    'Simulation', Opt.sim, ...
    'history', hisstr, ...
    'Layers', nlayer);

F = ncschema_adddims(F, ...
    'region',        nmask,    false, ...
    'ocean_time', length(t),    true);


Tinfo = ncinfo(avgfiles{1}, 'ocean_time');
tatts = {'long_name', 'units', 'calendar'};
[~,attloc] = ismember(tatts, {Tinfo.Attributes.Name});
tatts = [tatts; {Tinfo.Attributes(attloc).Value}];
tatts = tatts(:);

tunit = tatts{4};
tfile = cftime(t, tunit, [], 'reverse');

F = ncschema_addvars(F, ...
    'ocean_time', ...
    {'ocean_time'}, ...
    tatts, ...
    'double');
F = ncschema_addvars(F, ...
    'region', ...
    {'region'}, ...
    {'long_name', 'region number', ...
     'description', 'Regions based on AFSC groundfish survey strata.  Negative values indice a strata polygon was trimmed to eliminate regions beyond the modeled shelf break'}, ...
    'double');
F = ncschema_addvars(F, ...
    'region_area', ...
    {'region'}, ...
    {'long_name', 'region area', ...
     'units', 'km^2'}, ...
    'double');

for iv = 1:height(Vtbl)
    if isempty(Vtbl.unit{iv})
        atts = {'long_name', Vtbl.long{iv}};
    else
        atts = {'long_name', Vtbl.long{iv}, 'units', Vtbl.unit{iv}};
    end
    
    F = ncschema_addvars(F, ...
        Vtbl.short{iv}, ...
        {'region', 'ocean_time'}, ...
        atts, ...
        'double');
end

% Create file and add coordinate variable data

if ~exist(filereg, 'file')
    ncwriteschema(filereg, F);

    ncwrite(filereg, 'region', regnum);
    ncwrite(filereg, 'region_area', garea);
    ncwrite(filereg, 'ocean_time', tfile);
end

% SURVEY-REP

F = ncschema_init('classic');

hisstr = sprintf('%s: %s', ...
            datestr(now, 'ddd mmm dd HH:MM:SS yyyy'), ...
            'Survey-replicated ACLIM indices file created via roms_level3_aclimindices.m');
        
F = ncschema_addatts(F, ...
    'Simulation', Opt.sim, ...
    'history', hisstr, ...
    'Layers', nlayer);

F = ncschema_adddims(F, ...
    'station',  height(Svy),    false, ...
    'len',      size(char(Svy.STATIONID),2), false, ... 
    'year',     length(yr),     true);

% Coordinate variables: year, station details
        
F = ncschema_addvars(F, ...
     'year', ...
     {'year'}, ...
     {}, ...
     'double');
F = ncschema_addvars(F, ...              
    'station_id', ...
    {'len', 'station'}, ...
    {'long_name', 'station ID'}, ...
    'char');
F = ncschema_addvars(F, ...
    'latitude', ...
    {'station'}, ...
    {'long_name', 'station latitude'}, ...
    'double');
F = ncschema_addvars(F, ...
    'longitude', ...
    {'station'}, ...
    {'long_name', 'station longitude'}, ...
    'double');
F = ncschema_addvars(F, ...
    'stratum', ...
    {'station'}, ...
    {'long_name', 'station survey stratum'}, ...
    'double');
F = ncschema_addvars(F, ...
    'doy', ...
    {'station'}, ...
    {'long_name', 'sampling day-of-year (mean over 1982-2019 surveys)'}, ...
    'double');

for iv = 1:height(Vtbl)
    if isempty(Vtbl.unit{iv})
        atts = {'long_name', Vtbl.long{iv}};
    else
        atts = {'long_name', Vtbl.long{iv}, 'units', Vtbl.unit{iv}};
    end
    
    F = ncschema_addvars(F, ...
        Vtbl.short{iv}, ...
        {'station', 'year'}, ...
        atts, ...
        'double');
end

% Create file and add coordinate variable data

ncwriteschema(filesrep, F);

ncwrite(filesrep, 'year', yr);
ncwrite(filesrep, 'station_id', char(Svy.STATIONID)');
ncwrite(filesrep, 'latitude', Svy.LATITUDE);
ncwrite(filesrep, 'longitude', Svy.LONGITUDE);
ncwrite(filesrep, 'stratum', Svy.STRATUM);
ncwrite(filesrep, 'doy', Svy.DOY);

%--------------------
% Calculate indices
%--------------------

c = ConsoleProgressBar;
c.setMaximum(idxlim(2));
fprintf('Calculating indices...\n');
c.start();

for ii = 1:nvar

    c.setValue(ii);
    c.setText(sprintf('Index %2d/%2d', ii, nvar));

    % Files
    
    avgfiles = dir(fullfile(moxdir, 'roms_for_public', Opt.sim, sprintf('Level%d',Vtbl.level(ii)), sprintf('%s_*_average_%s.nc', Opt.sim, Vtbl.short{ii})));
    avgfiles = fullfile({avgfiles.folder}, {avgfiles.name});
    
    % Set up output arrays

    varreg  = nan(length(t), nmask);
    varsrep = nan(size(tsvy));

    cpflag = strcmp(Vtbl.short{ii}, 'temp_bottom5m');
    if cpflag
        [fracbelow0, fracbelow1, fracbelow2] = deal(nan(length(t),nmask));
    end

    for ia = 1:length(avgfiles)

        Tmp = ncstruct(avgfiles{ia}, Vtbl.internalshort{ii}, Scs);
        Tmp.t = ncdateread(avgfiles{ia}, 'ocean_time');

        % Regional averages

        [tf,loc] = ismember(Tmp.t,t);
        if ~all(tf)
            warning('Out of range time found in %s', avgfiles{ia});
        end
        for ir = 1:nmask
            varreg(loc(tf),ir) = local(Tmp.(Vtbl.internalshort{ii})(:,:,tf), gmask(:,:,ir), 'weight', gweight, 'omitnan');
        end

        % Calculate cold pool index fractions

        if cpflag
            for ir = 1:nmask
                fracbelow0(loc(tf),ir) = local(gweight.*double(Tmp.(Vtbl.internalshort{ii})(:,:,tf)<0), gmask(:,:,ir), @nansum)/garea(ir);
                fracbelow1(loc(tf),ir) = local(gweight.*double(Tmp.(Vtbl.internalshort{ii})(:,:,tf)<1), gmask(:,:,ir), @nansum)/garea(ir);
                fracbelow2(loc(tf),ir) = local(gweight.*double(Tmp.(Vtbl.internalshort{ii})(:,:,tf)<2), gmask(:,:,ir), @nansum)/garea(ir);
            end
        end

        % Survey rep

        usethisfile = fsvyidx == ia;

        xitmp = xisrep(usethisfile);
        etatmp = etasrep(usethisfile);
        ttmp = ftsvyidx(usethisfile);
        idx = sub2ind(size(Tmp.(Vtbl.internalshort{ii})), xitmp, etatmp, ttmp);
        varsrep(usethisfile) = Tmp.(Vtbl.internalshort{ii})(idx);

    end
    
    % Write to file
        
    ncwrite(filereg,  Vtbl.short{ii}, varreg');
    ncwrite(filesrep, Vtbl.short{ii}, varsrep');

    if cpflag
        ncwrite(filereg, 'fracbelow0', fracbelow0');
        ncwrite(filereg, 'fracbelow1', fracbelow1');
        ncwrite(filereg, 'fracbelow2', fracbelow2');
    end
end
c.stop();
fprintf('\n');






