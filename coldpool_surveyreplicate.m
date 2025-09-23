function coldpool_surveyreplicate(varargin)
%COLDPOOL_SURVEYREPLICATE Bering Sea survey-replication
%
% This function builds a .csv file holding data from one or more ocean
% models resampled to the same spatial and temporal patterns as seen in the
% annual Bering Sea bottom trawl survey (conducted by NOAA/AFSC/RACE's
% Groundfish Assessment Program). This survey-replication process is used
% to update model-based metrics of the cold pool and Bering Sea bottom
% temperature and to assess the skill of various models.
%
% Optional input variables (passed as parameter/value pairs, defaults in
% []): 
%
%   surveyfile: file name of survey data .csv file.  If empty, this points
%               to the data file in the GitHub-hosted version of the
%               afsc-gap-products/coldpool R package.
%               ['']
%
%   outputfile: file name where survey-replicated data should be written.
%               This file will be identical in format to the survey data
%               file but with additional columns repsenting
%               survey-replicated model data.  If file exists, data will be
%               appended or filled in as possible; otherwise, a new file
%               will be created.
%               ['coldpool_surveyreplicates_data.csv']
%
%   tempfiles:  simulation output file(s) (can be either a single file or a
%               set of files with a shared dimension along which they can
%               be concatenated) 
%               []
%
%   coordfile:  simulation coordinate file, file where lat/lon coordinates
%               for the temperature variable are held.  If empty, the first
%               tempfile will be used.
%               ['']
%
%   label:      Variable name to be used to label the survey-replicated
%               bottom temperature in the output table
%               ['model_bottom_temp']
%
%   vname:      variable name corresponding to temperature.  This must be
%               either a 3D or 4D variable (including the time dimension)
%               ['temp']
%
%   ltname:     variable name corresponding to latitude values
%               ['lat_rho']
%
%   lnname:     variable name corresponding to longitude values
%               ['lon_rho']
%
%   tname:      variable name corresponding to time values.  Time values
%               must conform to CF standards.
%               ['ocean_time']
%
%   bottomlayer:index of z-dimension corresponding to the bottom layer.  If
%               NaN, we assume data is 4D and on a z-coordinate grid, with
%               data below the bottom depth masked by NaNs.
%               [1]
%
%   verbose:    logical scalar, true to print progress to the command
%               window
%               [true]
%
%   mask:       logical array (or numeric that can be recasted) matching
%               the size of the lat/lon coordinate arrays in the file,
%               where 1 indicates water and 0 land.  This mask will be used
%               to remove any land-masked points from the nearest-neighbor
%               search.  If not included, points located close to land may
%               get matched to a grid cell with a placeholder value (e.g.
%               NaN). 
%               []

% Copyright 2024 Kelly Kearney

%--------------------
% Setup
%--------------------

% Input parsing

p = inputParser;
p.addParameter('surveyfile', '',   @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
p.addParameter('outputfile', 'coldpool_surveyreplicates_data.csv', @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
p.addParameter('coordfile', '',    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
p.addParameter('label', 'model_bottom_temp', @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
p.addParameter('tempfiles', [],    @(x) ~isempty(x) && (isstring(x) || iscellstr(x)));
p.addParameter('vname', 'temp',    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
p.addParameter('ltname', 'lat_rho', @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
p.addParameter('lnname', 'lon_rho', @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
p.addParameter('tname', 'ocean_time', @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
p.addParameter('bottomlayer', 1,   @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer'}));
p.addParameter('verbose', true, @(x) validateattributes(x, {'logical'}, {'scalar'}));
p.addParameter('mask', [], @(x) validateattributes(x, {'numeric', 'logical'}, {'nonnan'}));
p.parse(varargin{:});

Opt = p.Results;

% Check label: not conforming can lead to inconsistencies with table
% variable names, so let's not allow it 

if ~isvarname(Opt.label)
    error('Label must conform to Matlab variable name restrictions');
end

% Check mask

if ~islogical(Opt.mask)
    Opt.mask = logical(Opt.mask);
end

% Helper function

minmax = @(x) [min(x(:)) max(x(:))];

% Print start 
 
if Opt.verbose
    fprintf('Survey-replication: %s\n', Opt.label);
end

%--------------------
% Read survey data
%--------------------

if Opt.verbose
    fprintf('  Analyzing survey data...\n');
end

if isempty(Opt.surveyfile)

    % GitHub-hosted version of coldpool R package
    
    rawpath = 'https://raw.githubusercontent.com/afsc-gap-products/coldpool/main/assets/index_hauls_temperature_data.csv';
    Svy = readtable(rawpath);
else
    
    % Local file
    
    if ~exist(Opt.surveyfile, 'file')
        error('Could not locate specified input file');
    else
        Svy = readtable(Opt.surveyfile);
    end
end

%--------------------
% Read output file
%--------------------

if exist(Opt.outputfile, 'file')
    Old = readtable(Opt.outputfile);

    %--------------------
    % Verify survey data
    %--------------------

    % Check for updates between current survey data and the survey data 
    % portion of the output file.  There shouldn't be any changes other
    % than the addition of more recent data!

    simnames = setdiff(Old.Properties.VariableNames, Svy.Properties.VariableNames);
    Otmp = removevars(Old, simnames);

    [~,ia,ic] = intersect(Otmp(:,{'year','stationid'}), Svy(:,{'year','stationid'}));
    O1 = sortrows(Otmp(ia,:), {'year','stationid'}); % should be same
    N1 = sortrows(Svy(ic,:),  {'year','stationid'});
    
    if ~isequaln(O1,N1)
        error('Differences found between survey data stored in output file and survey file');
    end
    
    [~,ia,ic] = setxor(Otmp(:,{'year','stationid'}), Svy(:,{'year','stationid'}));
    O2 = Otmp(ia,:);
    
    if ~isempty(O2)
        error('There are survey data points (year/station ID) in output file not found in survey file');
    end

    %--------------------
    % New table
    %--------------------

    % Our updated table merges the old data, a new column for the new
    % simulation (if not already present), and any new rows of survey data
    % not present in the original.
    
    newval = num2cell(nan(length(ic), length(simnames)), 1);

    Svy = addvars(Svy(ic,:), newval{:}, 'NewVariableNames', simnames);
    Svy = Svy(:, Old.Properties.VariableNames); % make sure columns are in same order
    New = [Old; Svy];

    if ~ismember(Opt.label, simnames)
        New.(Opt.label) = nan(height(New),1);
    end
    
else
    New = Svy;
    New.(Opt.label) = nan(height(New),1);
end

needsdata = isnan(New.(Opt.label)); 
if ~any(needsdata)
    warning('Output file already includes complete replicates for %s', Opt.label);
    return
end

%---------------------
% Match points in time
% and space
%---------------------

if Opt.verbose
    fprintf('  Reading coordinate data...\n');
end

% Read coordinate data from input simulation

if isempty(Opt.coordfile)
    Opt.coordfile = Opt.tempfiles{1};
end

C.x = ncread(Opt.coordfile, Opt.lnname);
C.y = ncread(Opt.coordfile, Opt.ltname);
C.t = ncdateread(Opt.tempfiles, Opt.tname);

latlim = minmax(New.latitude(needsdata));
lonlim = minmax(wrapTo360(New.longitude(needsdata)));

needsdata = needsdata & New.start_time >= min(C.t) & New.start_time <= max(C.t); % skip any points outside sim time limits

if ~any(needsdata)
    warning('Output file already includes complete replicates (within simulation limits) for %s', Opt.label);
    return
end

% Temporal matching: nearest neighbor

if Opt.verbose
    fprintf('  Matching locations/times...\n');
end

tclose = interp1(C.t, 1:length(C.t), New.start_time(needsdata), 'nearest');
if any(isnan(tclose))
    error('Survey data outside model bounds?  Not caught by mask?'); % just in case...
end
tunq = unique(tclose)';
tblidx = find(needsdata);

% Spatial matching: convert to ENU space b/c it's a lot faster than doing
% the spheroid-based calcs, and the region is small enough 

cvt = @(lt,ln) geodetic2enu(lt,ln, 0, ...
    mean(latlim), mean(lonlim), 0, ...
    referenceEllipsoid('earth', 'km'));

if isvector(C.x) && isvector(C.y)
    [C.x, C.y] = ndgrid(C.x, C.y);
end

if isempty(Opt.mask)
    Opt.mask = true(size(C.x));
else
    if ~isequal(size(Opt.mask), size(C.x))
        error('mask size much match that of lat/lon coordinates');
    end
end

[xgrd, ygrd] = cvt(C.y(Opt.mask), C.x(Opt.mask));
[xs, ys] = cvt(New.latitude(needsdata), New.longitude(needsdata));

[~, gclose] = pdist2([xgrd(:) ygrd(:)], [xs(:) ys(:)], 'euclidean', 'smallest', 1);
wateridx = find(Opt.mask);
gclose = wateridx(gclose);

% Hyperslab: Only read the necessary spatial hyperslab
% Note: we assume that the variable in question includes two spatial
% dimensions, a time one, and possibly a depth one

[gr,gc] = ind2sub(size(C.x), gclose);
Ic = ncinfo(Opt.coordfile, Opt.lnname);
Iv = ncinfo(Opt.tempfiles{1}, Opt.vname);

dname = {Ic.Dimensions.Name Opt.tname}; % xyt
tf = ismember({Iv.Dimensions.Name}, dname);
if ~all(tf)
    if sum(~tf)>1
        error('Variable includes unexpected dimensions; can only handle 3D or 4D (including time) variables');
    else
        dname = [dname Iv.Dimensions(~tf).Name]; % xytz
    end
end
[~,prm] = ismember({Iv.Dimensions.Name}, dname);

rlim = minmax(gr);
clim = minmax(gc);

if length(dname) == 3 || (length(dname)==4 && isnan(Opt.bottomlayer))
    Scs = struct(dname{1}, [rlim(1) diff(rlim)+1 1], ...
                 dname{2}, [clim(1) diff(clim)+1 1]);
elseif length(dname) == 4
    Scs = struct(dname{1}, [rlim(1) diff(rlim)+1 1], ...
                 dname{2}, [clim(1) diff(clim)+1 1], ...
                 dname{4}, [Opt.bottomlayer 1 1]);
else
    error('Unexpected dimension (TODO: need to check what leads to this)')
end

gr2 = gr - rlim(1) + 1; % row index in spatial hyperslab
gc2 = gc - clim(1) + 1; % column index in spatial hyperslab
gclose2 = sub2ind([diff(rlim)+1 diff(clim)+1], gr2, gc2); % index in spatial hyperslab

%---------------------
% Read corresponding
% data points
%---------------------

if Opt.verbose
    fprintf('  Reading simulation data...\n');
end

% Read data one time slice at a time

nt = length(tunq);

if Opt.verbose
    c = ConsoleProgressBar;
    c.setMaximum(nt);
    c.start();
    count = 0;
end

for it = tunq
    if Opt.verbose
        count = count+1;
        c.setValue(count);
        c.setText(sprintf('%d/%d', count, nt));
    end

    Scs.(Opt.tname) = [it 1 1];

    % Extract bottom temp  

    A = ncstruct(Opt.tempfiles, Scs, Opt.vname);
    A.(Opt.vname) = permute(A.(Opt.vname), prm); % should be xyz now, with singleton t (and singleton z unless NaN bottomlayer)
    if isnan(Opt.bottomlayer)
        btemp = bottom(A.(Opt.vname));
    else
        btemp = A.(Opt.vname);
    end

    % Add data to points matching this time slice

    New.(Opt.label)(tblidx(tclose == it)) = btemp(gclose2(tclose == it));
end

if Opt.verbose
    c.stop();
    fprintf('\n');
end

%---------------------
% Write to file
%---------------------

writetable(New, Opt.outputfile);
