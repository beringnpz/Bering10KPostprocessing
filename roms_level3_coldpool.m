function roms_level3_coldpool(varargin)
%ROMS_LEVEL3_COLDPOOL Annual hindcast cold pool indices
%
% roms_level3_coldpool
% roms_level3_coldpool(param, val, ...)
%
% This function calculates a few different variations on the cold pool
% index.  Variations include the July 1 and survey-replicated model
% variants; the survey-based index is also included for comparison.
%
% This function expects to find survey data (and sampling coordinates) in
% the main roms_for_public folder.
%
% Output from this function is saved to the Level 3 folder of the specified
% simulation in a file named <simname>_coldpool.nc.  If a file with this
% name already exists, no action will be taken.
%
% Optional input variables (passed as parameter/value pairs)
%
%   simname:    string, name of simulation (in roms_for_public folder).
%               This function assumes that Level 2 bottom temperature files
%               have already been created for the indicated simulation.
%               ['B10K-K20_CORECFS']
%
%   svyfile:    name of Excel file holding survey data. The function will
%               first check to see if the file name resolves as a full or
%               relative path name, then will look for a file matching the
%               name in the primary roms_for_public folder.
%               ['AFSC_groundfish_survey_temperature_1982-2020.xlsx']

% Copyright 2021 Kelly Kearney

%--------------------
% Parse input
%--------------------

p = inputParser;
p.addOptional('sim', 'B10K-K20_CORECFS',  @(x) validateattributes(x, {'char'}, {'scalartext'}));
p.addOptional('svyfile', 'AFSC_groundfish_survey_temperature_1982-2020.xlsx',  @(x) validateattributes(x, {'char'}, {'scalartext'}));

p.parse(varargin{:});
Opt = p.Results;

fname = fullfile(moxdir, 'roms_for_public', Opt.sim, 'Level3', sprintf('%s_coldpool.nc', Opt.sim));
if exist(fname, 'file')
    error('Output file %s already exists; exiting', fname);
end

%--------------------
% Setup
%--------------------

grdfile = fullfile(moxdir, '/bering10k/input/grd/Bering_grid_withFeast.nc');

Grd = ncstruct(grdfile);
[nxi, neta] = size(Grd.h);

% [moxdir, Grd, nxi, neta, ~, ~, ~, Mask] = b10kdata;

if exist(Opt.svyfile, 'file')
    svyfile = Opt.svyfile;
else
    svyfile = fullfile(moxdir, 'roms_for_public', Opt.svyfile);
end
if ~exist(svyfile, 'file')
    error('Could not find indicated survey data file: %s\n', Opt.svyfile);
end

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

grdstrata1 = interpshapefile(Strata1, Grd.lat_rho, Grd.lon_rho, 'STRATUM');
grdstrata2 = interpshapefile(Strata2, Grd.lat_rho, Grd.lon_rho, 'EBS_STRATU');
isextra = isnan(grdstrata1) & ~isnan(grdstrata2) & grdstrata2 < 160; % 170 = west shelf, not using that one for this
grdstrata1(isextra) = grdstrata2(isextra);
grdstrata1(isnan(grdstrata1)) = 0;

isshelf = Grd.h <= 200;

gmask = cell(4,1);
gmask{1} = grdstrata1 >= 10 & grdstrata1 <= 62; % SEBS, no NW extension
gmask{2} = gmask{1} | ismember(grdstrata1, [82 90]); % SEBS+NW
gmask{3} = gmask{1} & isshelf; % SEBS, no shelf mismatch
gmask{4} = gmask{2} & isshelf; % SEBS+NW, no shelf mismatch
gmask = cat(3, gmask{:});

isebs = grdstrata1 > 0;

maskname = {...
    'SEBS'
    'SEBS+northwest'
    'SEBS, model shelf'
    'SEBS+northwest, model shelf'};

% Survey-replicated setup: include all stations on the primary
% 20-n.mi grid (for now, SEBS and SEBS+NW only, no crab resample)

SvySummary = readtable(svyfile, 'sheet', 'StationSummary');
isin = strcmp(SvySummary.TYPE, 'main') & ...
       ((SvySummary.STRATUM >= 10 & SvySummary.STRATUM <=62) | ismember(SvySummary.STRATUM, [82 90]));
SvySummary = SvySummary(isin,:);
   
Svy = readtable(svyfile, 'sheet', 'SurveyData');
keep = ismember(Svy.STATIONID, unique(SvySummary.STATIONID)) & Svy.BESTREP;
Svy = Svy(keep,:);

% Group stations into the defined SEBS regions

idx = sub2ind(size(Grd.h), SvySummary.B10K_XI, SvySummary.B10K_ETA);

srepmask(:,1) = SvySummary.STRATUM <= 62;
srepmask(:,2) = true(height(SvySummary),1);
srepmask(:,3) = srepmask(:,1) & Grd.h(idx) <= 200;
srepmask(:,4) = srepmask(:,2) & Grd.h(idx) <= 200;

% Regrid survey data to station x year grid

[SvyG.year, ~, iyr] = unique(year(Svy.DATETIME));
SvyG.station = SvySummary.STATIONID;
[tf,istation] = ismember(Svy.STATIONID, SvyG.station);
if ~all(tf)
    error('Station in main data not matching summary table; check this...');
end

SvyG.stratum = nan(max(iyr),max(istation));
idx = sub2ind([max(iyr) max(istation)], iyr, istation);

fld = {'GEAR_TEMPERATURE', 'STRATUM', 'LATITUDE', 'LONGITUDE', 'DATETIME'};
for ii = 1:length(fld)
    if strcmp(fld{ii}, 'DATETIME')
        SvyG.(fld{ii}) = NaT(max(iyr),max(istation));
    else
        SvyG.(fld{ii}) = nan(max(iyr),max(istation));
    end
    SvyG.(fld{ii})(idx) = Svy.(fld{ii});
end

% Thresholds

thresh = [0 1 1.5 2]; % deg C
nthresh = length(thresh);

% Survey-replicated grid indices

gmasksr = any(gmask,3) & Grd.mask_rho == 1;
gnum = find(gmasksr);

cvt = @(lt,ln) geodetic2enu(lt,ln, 0, ...
    mean(Grd.lat_rho(isebs)), mean(Grd.lon_rho(isebs)), 0, ...
    referenceEllipsoid('earth', 'km'));

[xgrd, ygrd] = cvt(Grd.lat_rho(gmasksr), Grd.lon_rho(gmasksr));
[xs, ys] = cvt(SvyG.LATITUDE, SvyG.LONGITUDE);

[d, gclose] = pdist2([xgrd(:) ygrd(:)], [xs(:) ys(:)], 'euclidean', 'smallest', 1);

[SvyG.xi, SvyG.eta] = ind2sub([nxi neta], reshape(gnum(gclose), size(xs)));
SvyG.xi(isnan(SvyG.LATITUDE)) = NaN;
SvyG.eta(isnan(SvyG.LATITUDE)) = NaN;

% Set up horizontal subsetting (to minimize data read)

gidx = find(any(gmask,3));
[ixi,ieta] = ind2sub([nxi neta], gidx);

xilim = minmax([ixi; SvyG.xi(:)]);
etalim = minmax([ieta; SvyG.eta(:)]);

gmask = gmask(xilim(1):xilim(2), etalim(1):etalim(2), :);

gweight = Grd.area_feast(xilim(1):xilim(2), etalim(1):etalim(2));

nmask = size(gmask,3);
garea = zeros(nmask,1);
for ii = 1:size(gmask,3)
    garea(ii) = sum(gweight(gmask(:,:,ii)));
end
    
Scs = struct('xi_rho',  [xilim(1)  diff(xilim)+1  1], ...
             'eta_rho', [etalim(1) diff(etalim)+1 1]);
         
xisrep = SvyG.xi - xilim(1) + 1;
etasrep = SvyG.eta - etalim(1) + 1;
gidxsrep = sub2ind([Scs.xi_rho(2) Scs.eta_rho(2)], xisrep, etasrep);

%--------------------
% Survey-based
%--------------------

nsvyyear = length(SvyG.year);
cpool1 = nan(nsvyyear, nmask, nthresh);
btemp1 = nan(nsvyyear, nmask);

for imask = 1:nmask
    tmp = SvyG.GEAR_TEMPERATURE(:,srepmask(:,imask));
    btemp1(:,imask) = nanmean(tmp,2);

    for ith = 1:length(thresh)
        cpool1(:,imask,ith) = sum(tmp<thresh(ith),2)./sum(~isnan(tmp),2);
    end
end

%--------------------
% Survey-replicated
%--------------------

% Extract times from files

F = dir(fullfile(moxdir, 'roms_for_public', Opt.sim, 'Level2', '*temp_bottom5m.nc'));
fname = fullfile({F.folder}', {F.name}');

t = cellfun(@(x) ncdateread(x, 'ocean_time'), fname, 'uni', 0);
tidx = cellfun(@(x) (1:length(x))', t, 'uni', 0);
fidx = cellfun(@(x,y) ones(size(x))*y, t, num2cell(1:length(t))', 'uni', 0);

[tunq, iunq] = unique(cat(1, t{:}), 'last');
tidx = cat(1, tidx{:}); % time index in file
fidx = cat(1, fidx{:}); % file index

% Match up survey dataset times to simulation times
% idx is index of t, tidx, fidx

idx = interp1(tunq, iunq, SvyG.DATETIME, 'nearest');

% Build survey-replicated bottom temperature array

btempsrep = nan(size(SvyG.GEAR_TEMPERATURE));

idxlim = minmax(idx);

c = ConsoleProgressBar;
c.setMaximum(idxlim(2));
fprintf('Extracting survey-replicated data...\n');
c.start();

for ii = idxlim(1):idxlim(2)
    c.setValue(ii);
    c.setText(sprintf('%4d/%4d', ii, idxlim(2)));
    
    ismatch = idx == ii & ~isnan(gidxsrep);
    
    if any(ismatch(:))
        
        Scs.ocean_time = [tidx(ii) 1 1];
        
        Data = ncstruct(fname{fidx(ii)}, 'temp', Scs);
        btempsrep(ismatch) = Data.temp(gidxsrep(ismatch));
       
    end
end
c.stop();
fprintf('\n');

% Calculate indices

cpool2 = nan(nsvyyear, nmask, nthresh);
btemp2 = nan(nsvyyear, nmask);

for imask = 1:nmask
    tmp = btempsrep(:,srepmask(:,imask));
    btemp2(:,imask) = nanmean(tmp,2);

    for ith = 1:length(thresh)
        cpool2(:,imask,ith) = sum(tmp<thresh(ith),2)./sum(~isnan(tmp),2);
    end
end

%--------------------
% July 1
%--------------------

yr = unique(year(tunq));

ttarget = datetime(yr,7,1);
idx = interp1(tunq, iunq, ttarget, 'nearest');
nt = length(ttarget);

% Build July temp dataset

btempjuly = nan(Scs.xi_rho(2), Scs.eta_rho(2), nt);

c.setMaximum(nt);
fprintf('Extracting July 1 data...\n');
c.setText(sprintf('%2d/%2d', 0, nt));
c.start();

for ii = 1:nt
    
    c.setValue(ii);
    c.setText(sprintf('%2d/%2d', ii, nt));
    
    Scs.ocean_time = [tidx(idx(ii)) 1 1];
    Data = ncstruct(fname{fidx(idx(ii))}, 'temp', Scs);
    btempjuly(:,:,ii) = Data.temp;
    
end
c.stop();
fprintf('\n');

% Calculate indices

cpool3 = nan(nt, nmask, nthresh);
btemp3 = nan(nt, nmask);

for imask = 1:nmask
    
    btemp3(:,imask) = local(btempjuly, gmask(:,:,imask), 'weight', gweight, 'omitnan'); 

    for ith = 1:length(thresh)
        cpool3(:,imask,ith) = local(gweight.*double(btempjuly<thresh(ith)), gmask(:,:,imask), @nansum)/garea(imask);
    end
end

%--------------------
% Write to file
%--------------------

A = struct;

tunit = 'days since 1900-01-01';
A.time = cftime(ttarget, tunit, [], 'reverse');

A.thresh = thresh;

A.method = strvcat('Model July', 'Survey', 'Model survey-replicated');
A.region = strvcat(maskname{:});

% Combine indices (note: remember, reordered relative to above)

[~,loc] = ismember(SvyG.year, year(ttarget));

A.btemp = nan(nt, nmask, 3);
A.cpool = nan(nt, nmask, nthresh, 3);

A.btemp(  :,:,1) = btemp3; % July
A.btemp(loc,:,2) = btemp1; % survey
A.btemp(loc,:,3) = btemp2; % survey-rep

A.cpool(  :,:,:,1) = cpool3; % July
A.cpool(loc,:,:,2) = cpool1; % survey
A.cpool(loc,:,:,3) = cpool2; % survey-rep

A.btemp = permute(A.btemp, [2 3 1]);   % region, method, time
A.cpool = permute(A.cpool, [3 2 4 1]); % thresh, region, method, time

hisstr = sprintf('%s: %s', ...
    datestr(now, 'ddd mmm dd HH:MM:SS yyyy'), ...
    'Cold pool indices file created via roms_level3_coldpool.m');

% NetCDF file schema setup

F = ncschema_init('classic');

F = ncschema_addatts(F, ...
    'Simulation', Opt.sim, ...
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


% Create file and add data

ncwriteschema(fname, F);

ncwrite(fname, 'time',                A.time);
ncwrite(fname, 'threshold',           A.thresh);
ncwrite(fname, 'region_label',        A.region');
ncwrite(fname, 'method_label',        A.method');
ncwrite(fname, 'average_bottom_temp', A.btemp);
ncwrite(fname, 'cold_pool_index',     A.cpool);

end

















