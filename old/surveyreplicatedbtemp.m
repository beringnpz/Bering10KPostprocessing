function surveyreplicatedbtemp(sim, varargin)
%SURVEYREPLICATEDBTEMP Calculates survey-replicated bottom temperatures
%
% surveyreplicatedbtemp(sim)
% surveyreplicatedbtemp(sim, p1, v1, ...)
%
% This function extracts bottom temperature data from a ROMS simulation
% using the same spatial and temporal sampling pattern as in the AFSC
% Bering Sea shelf groundfish survey.  
% 
% The data is written to a text file in the Level3 folder of the indicated
% simulation; data is a table replicating the .csv data file provided by
% the afsc-gap-products/coldpool index_hauls_temperature_data.csv file with
% an additional column added for model output. 
%
% Input variables:
%
%   sim:        simulation name, corresponding to the name used for
%               post-processing of the data
%
% Optional input variables: 
%
%   cpoolrepo:  path to local clone of the afsc-gap-products/coldpool
%               GitHub repository
%
%   grdfile:    path to ROMS grid file 

% Copyright 2023 Kelly Kearney

p = inputParser;
p.addParameter('cpoolrepo', '~/Documents/Research/Working/ReposCode/coldpool/');
p.addParameter('grdfile', fullfile(moxdir, '/roms_for_public/Bering10K_extended_grid.nc'));

p.parse(varargin{:});
Opt = p.Results;

% The survey-rep file and the actual survey data file

srepfile = fullfile(moxdir, 'roms_for_public', sim, 'Level3', sprintf('survey_replicates_%s.csv', sim));
svyfile = fullfile(Opt.cpoolrepo, 'data', 'index_hauls_temperature_data.csv');

if ~exist(svyfile, 'file')
    warning('Survey file not found; exiting without prepping survey replicates');
    return
end

srepexists = exist(srepfile, 'file');

if srepexists
    Old = readtable(srepfile);
    New = readtable(svyfile);

    % Check for updates

    Otmp = removevars(Old, 'model_bottom_temp');
    [~,ia,ic] = intersect(Otmp(:,{'year','stationid'}), New(:,{'year','stationid'}));
    O1 = sortrows(Otmp(ia,:), {'year','stationid'}); % should be same
    N1 = sortrows(New(ic,:),  {'year','stationid'});
    
    if ~isequaln(O1,N1)
        error('Differences found between stored data and new file');
    end
    
    [~,ia,ic] = setxor(Otmp(:,{'year','stationid'}), New(:,{'year','stationid'}));
    O2 = Otmp(ia,:);
    
    if ~isempty(O2)
        error('There are points in stored data not found in new file');
    end
    
    % TODO: need more robust way to check for updates to survey values
    % isnew = ~ismember(removevars(Old, 'model_bottom_temp'), New);
    Svy = New(ic,:);
else
    Svy = readtable(svyfile);
end

% Grid info

Grd = ncstruct(Opt.grdfile);
[nxi, neta] = size(Grd.h);

% Masks

M = cpindexmasks(Grd);

% Grid weighting

gmask = any(M.mask,3);

extractblock = @(x) x(any(gmask,2),any(gmask,1));

% Match points to closest model data

if ~isempty(Svy)

    F = dir(fullfile(moxdir, 'roms_for_public', sim, 'Level2', '*average_temp_bottom5m.nc'));
    fname = fullfile({F.folder}, {F.name});

    t = ncdateread(fname, 'ocean_time');

    cvt = @(lt,ln) geodetic2enu(lt,ln, 0, ...
        mean(Grd.lat_rho(gmask)), mean(Grd.lon_rho(gmask)), 0, ...
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

        btmp = readbtemp('sim', sim, 'dates', t(tunq(ii)), 'mask', M.mask);

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