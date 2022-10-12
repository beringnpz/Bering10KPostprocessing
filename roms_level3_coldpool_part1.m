function roms_level3_coldpool_part1(sim, yrs)
%ROMS_LEVEL3_COLDPOOL_PART1 Build/update cold pool metrics file
%
% roms_level3_coldpool_part1(sim)
%

if nargin < 2
    yrs = 1970:year(datetime('today'));
end

% outname = fullfile(moxdir, 'roms_for_public', sim, 'Level3', sprintf('%s_coldpool.nc', sim));
outname = 'test.nc';

outfileexists = exist(outname, 'file');

%---------------------
% Jul 1 metric
%---------------------

grdfile = fullfile(moxdir, '/roms_for_public/Bering10K_extended_grid.nc');

Grd = ncstruct(grdfile);
[nxi, neta] = size(Grd.h);

% Masks

gmask = cell(4,1);
gmask{1} = Grd.surveystrata_comboeast >= 10 & Grd.surveystrata_comboeast <= 62; % SEBS, no NW extension
gmask{2} = gmask{1} | ismember(Grd.surveystrata_comboeast, [82 90]); % SEBS+NW
gmask{3} = gmask{1} & Grd.h <= 200; % SEBS, no shelf mismatch
gmask{4} = gmask{2} & Grd.h <= 200; % SEBS+NW, no shelf mismatch
gmask = cat(3, gmask{:});

maskname = {...
    'SEBS-northwest'
    'SEBS+northwest'
    'SEBS-northwest, model shelf'
    'SEBS+northwest, model shelf'};

nmask = size(gmask,3);

% Grid weighting

gmaskany = any(gmask,3);

extractblock = @(x) x(any(gmaskany,2),any(gmaskany,1));

w = extractblock(Grd.area_feast);

% Thresholds

thresh = [-1 0 1 1.5 2]; % deg C
nthresh = length(thresh);

% Read data

ttarget = datetime(yrs, 7, 1);

if outfileexists
    tcomplete = ncdateread(outname, 'time');
    [tf, loc] = ismember(ttarget, tcomplete);
    ttarget = ttarget(~tf);
end

addjul1data = ~isempty(ttarget);

if addjul1data
    fprintf('Calculating July1 indices\n');
    [btemp, tmod, mask] = readbtemp('sim', sim, 'dates', ttarget, 'mask', gmask);

    % Calculate indices

    Idx = btemp2indices(btemp, mask, 'w', w, 'thresh', thresh);
end

%---------------------
% Create/update file
% and add Jul1 index
% data
%---------------------

if ~outfileexists
    
    A = struct;

    tunit = 'days since 1900-01-01';
    A.time = cftime(ttarget, tunit, [], 'reverse');

    A.thresh = thresh;

    A.method = char({'Model July', 'Survey', 'Model survey-replicated'});
    A.region = char(maskname);
    

    hisstr = sprintf('%s: %s', ...
        datestr(now, 'ddd mmm dd HH:MM:SS yyyy'), ...
        'Cold pool indices file created via roms_level3_coldpool_part1.m');

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

    % Create file and add data

    fprintf('Creating file...\n');
    
    ncwriteschema(outname, F);
    
    ncwrite(outname, 'time',                A.time);
    ncwrite(outname, 'threshold',           A.thresh);
    ncwrite(outname, 'region_label',        A.region');
    ncwrite(outname, 'method_label',        A.method');    
    ncwrite(outname, 'average_bottom_temp', permute(Idx.btemp, [2 3 1])); % Note: don't need to specify start b/c all = 1
    ncwrite(outname, 'cold_pool_index',     permute(Idx.cpool, [3 2 4 1]));
    
    if length(ttarget) > 1
        ncaddhis(outname, sprintf('Data added for years %d-%d', year(minmax(ttarget))));
    else
        ncaddhis(outname, sprintf('Data added for year %d', year(ttarget)));
    end
    
elseif addjul1data
    
    % Note: This will add new years after existing... mostly assuming that
    % we will stay in order, though that's not enforced.
    
    sidx = length(tcomplete)+1;
    
    tunit = 'days since 1900-01-01';
    tnc = cftime(ttarget, tunit, [], 'reverse');
    
    ncwrite(outname, 'time', tnc, sidx);
    ncwrite(outname, 'average_bottom_temp', permute(Idx.btemp, [2 3 1]), [1 1 sidx]);
    ncwrite(outname, 'cold_pool_index',     permute(Idx.cpool, [3 2 4 1]), [1 1 1 sidx]);
    
    if length(ttarget) > 1
        ncaddhis(outname, sprintf('Data added for years %d-%d', year(minmax(ttarget))));
    else
        ncaddhis(outname, sprintf('Data added for year %d', year(ttarget)));
    end
else
    fprintf('Jul 1 indices already up to date\n');
end

%---------------------
% Survey-rep
%---------------------

rproj = '~/Documents/Research/Working/BeringSea/cold_pool_index';
srepfile = fullfile(rproj, 'data', 'survey_replicates_B10K-K20_CORECFS.csv');

cpoolrepo = '~/Documents/Research/Working/ReposCode/coldpool/';
svyfile = fullfile(cpoolrepo, 'data', 'index_hauls_temperature_data.csv');

srepexists = exist(srepfile, 'file');

if srepexists
    Old = readtable(srepfile);
    New = readtable(svyfile);
    
    isnew = ~ismember(removevars(Old, 'model_bottom_temp'), New);
    Svy = New(isnew,:);
else
    Svy = readtable(svyfile);
end

% Match points to closest model data

if ~isempty(Svy)
    fprintf('Building survey-replicated data file\n');

    F = dir(fullfile(moxdir, 'roms_for_public', sim, 'Level2', '*average_temp_bottom5m.nc'));
    fname = fullfile({F.folder}, {F.name});

    t = ncdateread(fname, 'ocean_time');

    cvt = @(lt,ln) geodetic2enu(lt,ln, 0, ...
        mean(Grd.lat_rho(gmaskany)), mean(Grd.lon_rho(gmaskany)), 0, ...
        referenceEllipsoid('earth', 'km'));

    [xgrd, ygrd] = cvt(extractblock(Grd.lat_rho), extractblock(Grd.lon_rho));
    [xs, ys] = cvt(Svy.latitude, Svy.longitude);

    [~, gclose] = pdist2([xgrd(:) ygrd(:)], [xs(:) ys(:)], 'euclidean', 'smallest', 1);

    tclose = interp1(t, 1:length(t), Svy.start_time, 'nearest');
    checknan = isnan(tclose);
    if any(checknan)
        warning('Survey data found outside of model time bounds; dropping these points from survey-rep dataset');
        Svy = Svy(~checknan,:);
        gclose = gclose(~checknan);
        tclose = tclose(~checknan);
    end
    tunq = unique(tclose);

    % Read data points

    Svy.model_bottom_temp = nan(height(Svy),1);

    for ii = 1:length(tunq)

        isthistime = tclose == tunq(ii);

        btmp = readbtemp('sim', sim, 'dates', t(tunq(ii)), 'mask', gmask);

        Svy.model_bottom_temp(isthistime) = btmp(gclose(isthistime));

    end

    % Append new data to existing and write to file

    if srepexists
        Svy = [Old; Svy];
    end

    writetable(Svy, srepfile);
else
    fprintf('No new survey data found\n');
end










