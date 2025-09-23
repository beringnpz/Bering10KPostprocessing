function process5yrblocks(simbase, varargin)
%PROCESS5YRBLOCKS Produce post-processed level 1-2 output for ROMS sims
%
% This function runs the primary level 1 and level 2 processing steps for
% any ROMS simulation.  It is primarily intended to be used with the
% Bering10K hindcasts and ACLIM long term forecasts.  For each simulation,
% it collects output into 5-year blocks ranging from 1970-2100 (as
% applicable).  The following types of processing are available:
%
%   - Extract any variable from the raw output of average, history, station
%     files to the Level 1 collection
%   - Add pre-calculated depth values to the Level 2 collection
%   - Calculate surface-5m and/or bottom-5m average values for any 3D rho,
%     u, or v variable in the Level 1 or Level 2 collection, added to Level
%     2 collection
%   - Calculate depth-integrated values for any 3D rho, u, or v variable in
%     the Level 1 or Level 2 collection, added to Level 2 collection
%   - Calculate depth-integrated production values for all BESTNPZ
%     phytoplankton and zooplankton groups (note that this is handled
%     separately from other depth-integration due to the change in output
%     units between the H16 and K20 versions of BESTNPZ), added to Level 2
%     collection 
%   - Calculate geographically-rotated uEast and vNorth fields, added to
%     Level 2 collection 
%
% Input variables:
%
%   simbase:            base name for simulation, character array.  This
%                       will be used as the name of the folder under  
%                       the posprocessing folder location (see below) as
%                       well as for the base within all file names.  
%
% Optional input variables (passed as parameter/value pairs).  Defaults
% below each in [].
%
%   [a/h/s]files:       cell array of strings, full filename of raw output
%                       averages/history/station files to be processed
%                       [{}]
%
%   [a/h/s]lev1:        cell array of strings, names of variables to be
%                       extracted from the raw output files to the Level 1
%                       collection.  The special string 'allfluxes' can be
%                       used to extract all BESTNPZ output variables
%                       matching the 'fluxtype_source_target' format. 
%                       [{}]
%
%   [a/h/s]rotateuv:    logical scalar, true to add rotated uEast and
%                       vNorth files to Level 2 collection
%                       [false]
%
%   [a/h/s]addz:        logical scalar, true to add calculated depth values
%                       to the Level 2 collection
%                       [false]
%
%   [a/h/s]srf:         cell array of strings, names of 3D rho, u, or v
%                       variables from the exisiting Level 1 and 2
%                       collection for which surface-(dz)m averages should be
%                       calculated and added to the Level 2 collection 
%                       [{}]
%
%   [a/h/s]bot:         cell array of strings, names of 3D rho, u, or v
%                       variables from the exisiting Level 1 and 2
%                       collection for which bottom-(dz)m averages should be
%                       calculated and added to the Level 2 collection 
%                       [{}]
%
%   [a/h/s]int:         cell array of strings, names of 3D rho, u, or v
%                       variables from the exisiting Level 1 and 2
%                       collection for which depth-integrated values should
%                       be calculated and added to the Level 2 collection 
%                       [{}]
%
%   [a/h/s]intprod:     logical scalar, true to calculate depth-integrated
%                       versions of the production variables.  These need
%                       to be integrated as a special case, rather than
%                       being included in the [a/h/s]int list.
%
%   dryrun:             logical scalar, true to perform a dry run where no
%                       actions are actually taken, but instead the planned
%                       post-processing steps are printed.
%                       [true]
%
%   logid:              file ID (see fopen) where log messages should be
%                       written.  The default is to print to standard
%                       output (i.e. the command window).
%                       [1]
%
%   verbose:            logical scalar, true to use verbose mode when
%                       calling the main postprocessing steps.  This prints
%                       more detailed progress indicators as the functions
%                       run.
%                       [false]
%
%   dz:                 thickness of surface/bottom layer used for
%                       [a/h/s][srf/bot] extraction, in meters
%                       [5]
%
%   ppbase:             postprocessing folder base location, i.e. where the
%                       new simbase folder will be added.  Note that the
%                       default value here points to a deprecated machine;
%                       going forward, all calls should include an updated
%                       value for this.
%                       ['/gscratch/bumblereem/roms_for_public/']
%
%   grdfile:            path to grid file. Note that the default value here
%                       points to a deprecated machine; going forward, all
%                       calls should include an updated value for this.
%                       ['/gscratch/bumblereem/roms_for_public/Bering10K_extended_grid.nc']


% Copyright 2021 Kelly Kearney

%--------------
% Input parsing
%--------------

isstrings = @(x) isempty(x) || iscellstr(x) || isstring(x);

p = inputParser;
p.addParameter('afiles', {},  isstrings); 
p.addParameter('hfiles', {},  isstrings);
p.addParameter('sfiles', {},  isstrings);
p.addParameter('dfiles', {},  isstrings);

p.addParameter('alev1', {},  isstrings); % 'allfluxes' special
p.addParameter('hlev1', {},  isstrings);
p.addParameter('slev1', {},  isstrings);
p.addParameter('dlev1', {},  isstrings);

p.addParameter('arotateuv', false, @(x) validateattributes(x, {'logical'}, {'scalar'}));
p.addParameter('hrotateuv', false, @(x) validateattributes(x, {'logical'}, {'scalar'}));
p.addParameter('srotateuv', false, @(x) validateattributes(x, {'logical'}, {'scalar'}));
p.addParameter('drotateuv', false, @(x) validateattributes(x, {'logical'}, {'scalar'}));

p.addParameter('aaddz', false, @(x) validateattributes(x, {'logical'}, {'scalar'}));
p.addParameter('haddz', false, @(x) validateattributes(x, {'logical'}, {'scalar'}));
p.addParameter('saddz', false, @(x) validateattributes(x, {'logical'}, {'scalar'}));
p.addParameter('daddz', false, @(x) validateattributes(x, {'logical'}, {'scalar'}));

p.addParameter('asrf', {},  isstrings);
p.addParameter('hsrf', {},  isstrings);
p.addParameter('ssrf', {},  isstrings);
p.addParameter('dsrf', {},  isstrings);

p.addParameter('abot', {},  isstrings);
p.addParameter('hbot', {},  isstrings);
p.addParameter('sbot', {},  isstrings);
p.addParameter('dbot', {},  isstrings);

p.addParameter('aint', {},  isstrings);
p.addParameter('hint', {},  isstrings);
p.addParameter('sint', {},  isstrings);
p.addParameter('dint', {},  isstrings);

p.addParameter('aintprod', false, @(x) validateattributes(x, {'logical'}, {'scalar'}));
p.addParameter('hintprod', false, @(x) validateattributes(x, {'logical'}, {'scalar'}));
p.addParameter('sintprod', false, @(x) validateattributes(x, {'logical'}, {'scalar'}));
p.addParameter('dintprod', false, @(x) validateattributes(x, {'logical'}, {'scalar'}));

p.addParameter('dryrun', true, @(x) validateattributes(x, {'logical'}, {'scalar'}));
p.addParameter('logid', 1, @(x) validateattributes(x, {'double'}, {'scalar','integer'}));
p.addParameter('verbose', false, @(x) validateattributes(x, {'logical'}, {'scalar'}));
p.addParameter('dz', 5, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'nonnegative'}));

p.addParameter('ppbase', '/gscratch/bumblereem/roms_for_public/', @(x) validateattributes(x, {'char','string'}, {'scalartext'}))
p.addParameter('grdfile', fullfile(mounteddir('klone'), 'GR011846_reem/postprocessed', 'Bering10K_extended_grid.nc'))

p.parse(varargin{:});
Opt = p.Results;

% Collect parameters by file type

ntype = 4;

filetype = {'average', 'history', 'station', 'diagnos'};
prefix = {'a', 'h', 's', 'd'};

param = {'files', 'lev1', 'srf', 'bot', 'int', 'addz', 'rotateuv', 'intprod'};
for itype = ntype:-1:1
    for ip = 1:length(param)
        P(itype).(param{ip}) = p.Results.([prefix{itype} param{ip}]);
    end
end

% Post-processed file destination

if ~exist(Opt.ppbase, 'dir')
    error('Postprocessing folder base location (%s) not found.  Check ppbase input option.', Opt.ppbase);
end
ppfol = fullfile(Opt.ppbase, simbase);

% Log setup

if ~exist(Opt.grdfile, 'file')
    error('Could not locate grid file (%s)', Opt.grdfile);
end

for itype = 1:ntype
    
    fprintf(Opt.logid, 'Beginning post-processing for %s %s files...\n', simbase, filetype{itype});
    
    if isempty(P(itype).files)
        fprintf('  No raw output files found\n');
        continue;
    end
    
    %------------
    % Setup
    %------------
    
    % Categorize variables that exist in the raw data files
    % Note: production variables are treated specially because of the
    % way they need to be integrated/summed based on model version.

    I = ncinfo(P(itype).files{1});

    if ismember('allfluxes', P(itype).lev1) % special flag for all flux variables

        isflux = regexpfound({I.Variables.Name}, '(\w*)_(\w*)_(\w*)');

        P(itype).lev1 = union(...
            setdiff(P(itype).lev1, 'allfluxes', 'stable'), ...
            {I.Variables(isflux).Name}, 'stable');
    end

    haslev1 = ismember(P(itype).lev1, {I.Variables.Name});

    prodvars = P(itype).lev1(haslev1 & ...
               startsWith(P(itype).lev1, 'prod') & ...
               ~contains(P(itype).lev1, 'Ben') & ...
               ~contains(P(itype).lev1, 'Ice'));

    P(itype).int = setdiff(P(itype).int, prodvars, 'stable');
    
    % Set up 5-year blocks

    tedge = datetime(1970:5:2100,1,1);

    ngrp = length(tedge)-1;
    
    % Query time bounds for all files

    fprintf(Opt.logid, '  Querying raw file time bounds...\n');

    filetlim = ncdatelim(P(itype).files, 'ocean_time');
    
    %------------
    % Processing
    %------------

    for ii = 1:ngrp

        try
            basegrp = sprintf('%s_%d-%d', simbase, year(tedge(ii)), year(tedge(ii+1)-days(1)));

            lev1folder = fullfile(ppfol, 'Level1');
            lev2folder = fullfile(ppfol, 'Level2');

            % Check for existing files

            [l1v,l2v] = existinglev1and2(ppfol, basegrp, filetype{itype});

            % Figure out which station, history, and average files fall in the
            % specified range 

%             frac = fracoverlap(filetlim, tedge([ii ii+1]));
%             if ~any(frac>0)
%                 continue; % Skip remaining calcs if no files in this time reange
%             end
            frac = fileininterval(filetlim, tedge([ii ii+1]));
            if ~any(frac)
                continue;
            end

            grpavg = P(itype).files(frac); 

            % Call level 1 postprocessing for indicated variables
            % (Note: it would be a little more efficient to do all file
            % types at once, but this keeps the script more flexible).

            l1missing = setdiff(P(itype).lev1(haslev1), l1v, 'stable');

            if ~isempty(l1missing)
                fprintf(Opt.logid, '  Level 1 processing, time block %d/%d: %s\n', ii, ngrp, basegrp);
                if Opt.dryrun
                    fprintf(Opt.logid, '   Will add: ');
                    fprintf(Opt.logid, '%s, ', l1missing{:});
                    fprintf(Opt.logid, '\n');
                else
                    if itype == 1
                        roms_level1(lev1folder, basegrp, Opt.grdfile, ...
                            'avg', grpavg, ...
                            'his', [], ...
                            'sta', [], ...
                            'variables', l1missing, ...
                            'tbound', [tedge(ii) tedge(ii+1)], ...
                            'verbose', Opt.verbose);
                    elseif itype == 2
                        roms_level1(lev1folder, basegrp, Opt.grdfile, ...
                            'avg', [], ...
                            'his', grpavg, ...
                            'sta', [], ...
                            'variables', l1missing, ...
                            'tbound', [tedge(ii) tedge(ii+1)], ...
                            'verbose', Opt.verbose);
                    elseif itype == 3
                        roms_level1(lev1folder, basegrp, Opt.grdfile, ...
                            'avg', [], ...
                            'his', [], ...
                            'sta', grpavg, ...
                            'variables', l1missing, ...
                            'tbound', [tedge(ii) tedge(ii+1)], ...
                            'verbose', Opt.verbose);
                    elseif itype == 4
                        roms_level1(lev1folder, basegrp, Opt.grdfile, ...
                            'avg', [], ...
                            'his', [], ...
                            'sta', [], ...
                            'dia', grpavg, ...
                            'variables', l1missing, ...
                            'tbound', [tedge(ii) tedge(ii+1)], ...
                            'verbose', Opt.verbose);
                        
                    end
                end
            end

            % Level 2: Add pre-calculated depths variables

            if P(itype).addz
                hasz = ismember({'z_psi', 'z_psi_w', 'z_rho', 'z_u', 'z_u_w', ...
                                 'z_v', 'z_v_w', 'z_w'}, l2v);
                hasz = (itype<3 && all(hasz)) || (itype==3 && all(hasz([3 8])));
                if ~hasz
                    fprintf(Opt.logid, '  Level 2 depths, time block %d/%d: %s\n', ii, ngrp, basegrp);
                    if Opt.dryrun
                        fprintf(Opt.logid, '   Will add depth variables\n');
                    else
                        roms_level2_depths(Opt.grdfile, lev1folder, lev2folder, basegrp, ...
                            'verbose', Opt.verbose, 'ftype', filetype{itype});
                    end
                end
            end

            % Level 2: Rotate u and v and move to common grid

            [l1v,l2v] = existinglev1and2(ppfol, basegrp, filetype{itype});
            
            if P(itype).rotateuv

                uvcheck = ismember({'u','v'}, l1v) & ~ismember({'uEast', 'vNorth'}, l2v);
                if all(uvcheck)
                    fprintf(Opt.logid, '  Level 2 u/v rotation, time block %d/%d: %s\n', ii, ngrp, basegrp);
                    
                    if Opt.dryrun
                        fprintf(Opt.logid, '   Will add uEast and vNorth\n');
                    else
                        roms_level2_rotateuv(Opt.grdfile, lev1folder, lev2folder, basegrp, ...
                            'verbose', Opt.verbose, 'ftype', 'average', 'vtype', 'uv');
                    end
                end
                
                uvcheck = ismember({'sustr','svstr'}, l1v) & ~ismember({'sustrEast', 'svstrNorth'}, l2v);
                if all(uvcheck)
                    fprintf(Opt.logid, '  Level 2 sustr/svstr rotation, time block %d/%d: %s\n', ii, ngrp, basegrp);
                    
                    if Opt.dryrun
                        fprintf(Opt.logid, '   Will add sustrEast and svstrNorth\n');
                    else
                        roms_level2_rotateuv(Opt.grdfile, lev1folder, lev2folder, basegrp, ...
                            'verbose', Opt.verbose, 'ftype', 'average', 'vtype', 'suvstr');
                    end
                end
                
            end

            % Level 2: extract bottom/surface values

            [l1v,l2v] = existinglev1and2(ppfol, basegrp, filetype{itype});
            
            srfbot = reshape(union(P(itype).srf, P(itype).bot), 1, []);
            dzstr = sprintf('%dm', Opt.dz);
            isdone =  ismember([cellfun(@(x) [x '_surface' dzstr], srfbot', 'uni', 0), ...
                                cellfun(@(x) [x '_bottom' dzstr ], srfbot', 'uni', 0)], ...
                                l2v);
            if Opt.dryrun
                l1possible = [l1v l1missing];
            else
                l1possible = l1v;
            end
            canbedone = ismember(srfbot', l1possible) | ismember(srfbot', l2v); % parent file
            isneeded = [ismember(srfbot', P(itype).srf) ismember(srfbot', P(itype).bot)] & ...
                        ~isdone & canbedone;         

            if any(any(isneeded,2))  
                fprintf(Opt.logid, '  Level 2 surface/bottom averaging, time block %d/%d: %s\n', ii, ngrp, basegrp);
                
                if Opt.dryrun
                    srfstr = sprintf('%s,', srfbot{isneeded(:,1)});
                    botstr = sprintf('%s,', srfbot{isneeded(:,2)});
                    fprintf(Opt.logid, '   Will add surface for: %s\n', srfstr(1:end-1));
                    fprintf(Opt.logid, '   Will add bottom for: %s\n', botstr(1:end-1));
                else
                    for iv = 1:length(srfbot)
                        if any(isneeded(iv,:))

                            if all(isneeded(iv,:))
                                loc = 'both';
                            elseif isneeded(iv,1)
                                loc = 'top';
                            elseif isneeded(iv,2)
                                loc = 'bottom';
                            end

                            par1 = fullfile(lev1folder, sprintf('%s_%s_%s.nc', basegrp, filetype{itype}, srfbot{iv}));
                            par2 = fullfile(lev2folder, sprintf('%s_%s_%s.nc', basegrp, filetype{itype}, srfbot{iv}));
                            if exist(par1, 'file')
                                parentfile = par1;
                                hasparent = true;
                            elseif exist(par2, 'file')
                                parentfile = par2;
                                hasparent = true;
                            else
                                hasparent = false;
                            end

                            if hasparent
                                roms_level2_bottomsurf(...
                                    parentfile, ...
                                    lev2folder, ...
                                    [basegrp '_' filetype{itype}], ...
                                    fullfile(lev2folder, sprintf('%s_%s', basegrp, filetype{itype})), ...
                                    srfbot{iv}, ...
                                    'verbose', Opt.verbose, ...
                                    'location', loc, ...
                                    'dz', Opt.dz);
                            end
                        end
                    end 
                end
            end

            % Level 2: depth-integrated values

            canbedone = ismember(P(itype).int, l1possible) | ismember(P(itype).int, l2v); % parent file
            isdone =  ismember(cellfun(@(x) [x '_integrated'], P(itype).int, 'uni', 0), l2v) | ~canbedone;

            if ~all(isdone)
                fprintf(Opt.logid, '  Level 2 depth integration, time block %d/%d: %s\n', ii, ngrp, basegrp);
                if Opt.dryrun
                    totstr = sprintf('%s,', P(itype).int{~isdone});
                    fprintf(Opt.logid, '   Will add integrated for: %s\n', totstr(1:end-1));
                else
                    for iv = 1:length(P(itype).int)
                        if ~isdone(iv)

                            par1 = fullfile(lev1folder, sprintf('%s_%s_%s.nc', basegrp, filetype{itype}, P(itype).int{iv}));
                            par2 = fullfile(lev2folder, sprintf('%s_%s_%s.nc', basegrp, filetype{itype}, P(itype).int{iv}));
                            if exist(par1, 'file')
                                parentfile = par1;
                                hasparent = true;
                            elseif exist(par2, 'file')
                                parentfile = par2;
                                hasparent = true;
                            else
                                hasparent = false;
                            end

                            if hasparent
                                roms_level2_depthintegrate(...
                                    parentfile, ...
                                    lev2folder, ...
                                    sprintf('%s_%s', basegrp, filetype{itype}), ...
                                    fullfile(lev2folder, sprintf('%s_%s', basegrp, filetype{itype})), ...
                                    P(itype).int{iv}, ...
                                    'verbose', Opt.verbose);
                            end
                        end
                    end
                end
            end

            % Level 2: depth-integrated production.
            % Note: Code for prod. calculations changed significantly
            % between H16 and K20.  Most notably, K20 uses integrated units
            % per layer, while H16 used volumetric.  Also, there was some
            % weirdness in the NCa/Eup calcs, where a single value summed
            % over both sub-groups... not necessarily incorrect but perhaps
            % a little confusing?

            if P(1).intprod || P(4).intprod

                ish16 = contains(simbase, 'H16');
                isdone =  ismember(cellfun(@(x) [x '_integrated'], prodvars, 'uni', 0), l2v);

                if ~all(isdone)
                    fprintf(Opt.logid, '  Level 2 prod. integration, time block %d/%d: %s\n', ii, ngrp, basegrp);

                    if Opt.dryrun
                        totstr = sprintf('%s,', prodvars{~isdone});
                        fprintf(Opt.logid, '   Will add integrated for: %s\n', totstr(1:end-1));
                    else
                        for iv = 1:length(prodvars)
                            if ~isdone(iv)
                                roms_level2_depthintegrate(...
                                    fullfile(lev1folder, sprintf('%s_%s_%s.nc', basegrp, filetype{itype}, prodvars{iv})), ...
                                    lev2folder, ...
                                    sprintf('%s_%s', basegrp, filetype{itype}), ...
                                    fullfile(lev2folder, sprintf('%s_%s', basegrp, filetype{itype})), ...
                                    prodvars{iv}, ...
                                    'verbose', Opt.verbose, ...
                                    'simplesum', ~ish16);
                            end
                        end
                    end
                end
            end

        catch ME
            fprintf(Opt.logid, '%s: An error occurred processing %s files, time block %d (%s)\n', datestr(now), filetype{itype}, ii, basegrp);
            rethrow(ME);
        end
    end
    
end

end

%-----------------------
% Subfunctions
%-----------------------

function isin = fileininterval(tlim, tedge)

isin = ~(tlim(:,2) < tedge(1) | tlim(:,1) > tedge(2)) & ~any(isnat(tlim),2);

end

% Subfunction: fraction of file in time interval

function f = fracoverlap(tlim, tedge)

daysin = max(days(min(tedge(2),tlim(:,2)) - max(tedge(1),tlim(:,1))),0);
f = daysin./days(tlim(:,2)-tlim(:,1));
    
end

% Subfunction: which variables are in the Level 1-2 folders?

function [l1v,l2v] = existinglev1and2(folder, basegrp, ftype)

L1 = dir(fullfile(folder, 'Level1', [basegrp '*.nc']));
L2 = dir(fullfile(folder, 'Level2', [basegrp '*.nc']));

filter = sprintf('_%s_', ftype);
is1 = contains({L1.name}, filter);
is2 = contains({L2.name}, filter);

l1v = strrep(strrep({L1(is1).name}, [basegrp filter], ''), '.nc', '');
l2v = strrep(strrep({L2(is2).name}, [basegrp filter], ''), '.nc', '');

end

