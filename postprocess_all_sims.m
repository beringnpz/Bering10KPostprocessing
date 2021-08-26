% Level 1-2 Postprocessing for all simulations
%
% This is the primary post-processing script for all simulations destined
% for the roms_for_public folder.
%
% This script was originally developed to run postprocessing on all the
% simulations under the ACLIM umbrella, but has now been extended so it
% should be able to apply to any simulation.  

% Copyright 2021 Kelly Kearney

% Note: To run, user should modify/check "User flags" and "What to process"
% sections.

%------------
% User flags
%------------

% true to print progress bar (turn off in parallel)

vflag = true;  

% true for dry run (check for files and print intended action, but make no
% changes) 

dryrunflag = true;  

% true to run in parallel
 
parflag = false;   

% true to print to log file instead of screen (dry run mode only,
% recommended if processing more than a handful of variables or
% simulations, since the log gets pretty long)

uselog = false;      
      
%------------
% Setup
%------------

% Parallel setup

if parflag
    numCores = feature('numcores');
    parpool(numCores);
    maxworker = numCores;
else
    maxworker = 0;
end

if dryrunflag && uselog
    logfile = sprintf('postprocess_%s.log', datestr(now, 'yyyymmddHHMM'));
    fid = fopen(logfile, 'wt');
else
    fid = 1;
end

% Some setup info

% moxdir = b10kdata;
grdfile = fullfile(moxdir, '/bering10k/input/grd/Bering_grid_withFeast.nc');

% Add processing toolbox to path

if isempty(which('roms_level1'))
    addpath(fullfile(moxdir, 'bering10k', 'output', 'processing_scripts'));
end

%-----------------
% What to process
%-----------------

% Processing details

% Set the things to be extracted/calculated
%
%   lev1:       cell array of strings, short names of Level 1 variables to
%               be extracted  
%   rotateuv:   scalar logical, true to rotate u/v fields to geographic
%   addz:       scalar logical, true to add precalculated depth variables
%   srf:        cell array of strings, short names of Level 1 or 2
%               variables for which you want to extract surface-5m values
%   bot:        cell array of strings, short names of Level 1 or 2
%               variables for which you want to extract bottom-5m values
%   tot:        cell array of strings, short names of Level 1 or 2
%               variables for which you want to extract depth-integrated
%               values (note: do not include production variables here;
%               they are handled separately)
%   integratedprod: scalar logical, true to calculate integrated values for
%               all production variables

% TODO: need to change this to not use structure (problem for parfor)

% Average file processing

Vars.lev1{1} = {'temp', 'salt', ...
             'u', 'v', 'ubar', 'vbar', ...
             'aice', 'hice', ...
             'Hsbl', ...
             'shflux', 'ssflux', ...
             'IceNH4', 'IceNO3', 'IcePhL', ...
             'PhS', 'PhL', ...
             'MZL', 'Cop', 'NCaS', 'NCaO', 'EupS', 'EupO', ...
             'NO3', 'NH4', 'Iron', ...
             'Det', 'DetF', 'DetBen', ...
             'Jel', 'Ben', ...
             'prod_PhS', 'prod_PhL', 'prod_MZL', 'prod_Cop', 'prod_NCaS', ...
             'prod_EupS', 'prod_NCaO', 'prod_EupO', 'prod_Jel', 'prod_Ben', ...
             'prod_IcePhL', 'prod_NCa', 'prod_Eup', ...
             'TIC', 'alkalinity', 'oxygen', ...
             'sustr', 'svstr'};
Vars.rotateuv(1) = true;
Vars.addz(1) = true;
Vars.srf{1} = {'temp', 'salt', 'PhS', 'PhL', 'MZL', 'Cop', 'NCaS', 'NCaO', 'EupS', 'EupO', 'Jel', 'NO3', 'NH4', 'Iron', 'uEast', 'vNorth'};
Vars.bot{1} = {'temp', 'NO3', 'NH4', 'Iron', 'uEast', 'vNorth'};
Vars.tot{1} = {'temp', 'PhS', 'PhL', 'MZL', 'Cop', 'NCaS', 'NCaO', 'EupS', 'EupO', 'Jel', 'NO3', 'NH4', 'Iron'}; 
Vars.integrateprod(1) = false;

% Vars.lev1{1} = {};
% Vars.rotateuv(1) = false;
% Vars.addz(1) = false;
% Vars.srf{1} = {};
% Vars.bot{1} = {};
% Vars.tot{1} = {};
% Vars.integrateprod(1) = false;

% History file processing

Vars.lev1{2} = {'temp', 'u', 'v'};
Vars.rotateuv(2) = false;
Vars.addz(2) = true;
Vars.srf{2} = {};
Vars.bot{2} = {};
Vars.tot{2} = {};
Vars.integrateprod(2) = false;

% Station file processing

Vars.lev1{3} = {'temp', 'salt', ...
                'PhS', 'PhL', 'MZL', 'Cop', 'NCaS', 'NCaO', ...
                'EupS', 'EupO', 'Jel', 'Ben', ...
                'NO3', 'NH4', 'Iron'};
Vars.rotateuv(3) = false;
Vars.addz(3) = true;
Vars.srf{3} = {};
Vars.bot{3} = {};
Vars.tot{3} = {};
Vars.integrateprod(3) = false;

% Vars.lev1{3} = {};
% Vars.rotateuv(3) = false;
% Vars.addz(3) = false;
% Vars.srf{3} = {};
% Vars.bot{3} = {};
% Vars.tot{3} = {};
% Vars.integrateprod(3) = false;

% Simulations: base name, path to output folder, filters for files

b10koutfol = fullfile(moxdir, 'bering10k', 'output');

simdetails = {
% 'B10K-K20P19_CMIP6_cesm_historical'  , fullfile(b10koutfol, 'forecasts/cmip6/cesm_historical/Out')        ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
% 'B10K-K20P19_CMIP6_cesm_ssp126'      , fullfile(b10koutfol, 'forecasts/cmip6/cesm_ssp126/Out')            ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
% 'B10K-K20P19_CMIP6_cesm_ssp585'      , fullfile(b10koutfol, 'forecasts/cmip6/cesm_ssp585/Out')            ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
% 'B10K-K20P19_CMIP6_gfdl_historical'  , fullfile(b10koutfol, 'forecasts/cmip6/gfdl_historical/Out')        ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
% 'B10K-K20P19_CMIP6_gfdl_ssp126'      , fullfile(b10koutfol, 'forecasts/cmip6/gfdl_ssp126/Out')            ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
% 'B10K-K20P19_CMIP6_gfdl_ssp585'      , fullfile(b10koutfol, 'forecasts/cmip6/gfdl_ssp585/Out')            ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
% 'B10K-K20P19_CMIP6_miroc_historical' , fullfile(b10koutfol, 'forecasts/cmip6/miroc_historical/Out')       ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
% 'B10K-K20P19_CMIP6_miroc_ssp126'     , fullfile(b10koutfol, 'forecasts/cmip6/miroc_ssp126/Out')           ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
% 'B10K-K20P19_CMIP6_miroc_ssp585'     , fullfile(b10koutfol, 'forecasts/cmip6/miroc_ssp585/Out')           ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
'B10K-K20_CORECFS'                   , fullfile(b10koutfol, 'hindcasts/npz_201904_delta30/Out')           ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
% 'B10K-H16_CORECFS'                   , fullfile(b10koutfol, 'hindcasts/npz_201904_aclim/Out/')            ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
% 'B10K-K20P19_CMIP5_CESM_rcp45'       , fullfile(b10koutfol, 'cmip5-avg-carbon/CESM_rcp45/')               ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
% 'B10K-K20P19_CMIP5_CESM_rcp85'       , fullfile(b10koutfol, 'cmip5-avg-carbon/CESM_rcp85/')               ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
% 'B10K-K20P19_CMIP5_GFDL_rcp45'       , fullfile(b10koutfol, 'cmip5-avg-carbon/GFDL_rcp45/')               ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
% 'B10K-K20P19_CMIP5_GFDL_rcp85'       , fullfile(b10koutfol, 'cmip5-avg-carbon/GFDL_rcp85/')               ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
% 'B10K-K20P19_CMIP5_MIROC_rcp45'      , fullfile(b10koutfol, 'cmip5-avg-carbon/MIROC_rcp45/')              ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
% 'B10K-K20P19_CMIP5_MIROC_rcp85'      , fullfile(b10koutfol, 'cmip5-avg-carbon/MIROC_rcp85/')              ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
% 'B10K-H16_CMIP5_CESM_BIO_rcp85'      , fullfile(b10koutfol, 'cmip5-avg/CESM_BIO_rcp85/')                  ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
% 'B10K-H16_CMIP5_CESM_rcp45'          , fullfile(b10koutfol, 'cmip5-avg/CESM_rcp45/')                      ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
% 'B10K-H16_CMIP5_CESM_rcp85'          , fullfile(b10koutfol, 'cmip5-avg/CESM_rcp85/')                      ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
% 'B10K-H16_CMIP5_GFDL_BIO_rcp85'      , fullfile(b10koutfol, 'cmip5-avg/GFDL_BIO_rcp85/')                  ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
% 'B10K-H16_CMIP5_GFDL_rcp45'          , fullfile(b10koutfol, 'cmip5-avg/GFDL_rcp45/')                      ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
% 'B10K-H16_CMIP5_GFDL_rcp85'          , fullfile(b10koutfol, 'cmip5-avg/GFDL_rcp85/')                      ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
% 'B10K-H16_CMIP5_MIROC_rcp45'         , fullfile(b10koutfol, 'cmip5-avg/MIROC_rcp45/')                     ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
% 'B10K-H16_CMIP5_MIROC_rcp85'         , fullfile(b10koutfol, 'cmip5-avg/MIROC_rcp85/')                     ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
% 'B10K-K20nobio_CFS'                  , fullfile(b10koutfol, 'hindcasts/phys_201910_30_singlesource/Out/') ,'physhccfs*_avg_*.nc'   ,'physhccfs*_his_*.nc'  ,'physhccfs*_sta_*.nc' 
% 'B10K-K20nobio_CORE'                 , fullfile(b10koutfol, 'hindcasts/phys_201910_30_singlesource/Out/') ,'physhccore*_avg_*.nc'  ,'physhccore*_his_*.nc' ,'physhccore*_sta_*.nc'
};

simdetails = cell2table(simdetails, 'VariableNames', {'base','path', 'afilter', 'hfilter', 'sfilter'});
nsim = height(simdetails);

for isim = 1:nsim
    simdetails.afiles{isim} = dirfull(fullfile(simdetails.path{isim}, simdetails.afilter{isim}));
    simdetails.hfiles{isim} = dirfull(fullfile(simdetails.path{isim}, simdetails.hfilter{isim}));
    simdetails.sfiles{isim} = dirfull(fullfile(simdetails.path{isim}, simdetails.sfilter{isim}));
end

% New postprocessing parent folders

simdetails.ppfol = fullfile(moxdir, 'roms_for_public', simdetails.base);

% Make sure all necessary folders exist, creating Level1/2 folders if
% needed

ppfol = unique(simdetails.ppfol);
for ii = 1:length(ppfol)
    if ~exist(ppfol{ii}, 'dir')
        error('Check postprocessing parent folders; %s not found', ppfol{ii});
    end
    if ~exist(fullfile(ppfol{ii}, 'Level1'), 'dir')
        mkdir(fullfile(ppfol{ii}, 'Level1'));
    end
    if ~exist(fullfile(ppfol{ii}, 'Level2'), 'dir')
        mkdir(fullfile(ppfol{ii}, 'Level2'));
    end
end

%--------------------
% Processing
%--------------------

filetype = {'average', 'history', 'station'};

% parfor (isim = 1:nsim, maxworker)
for isim = 1:nsim % for debugging

    for itype = 1:3
        if itype == 1
            outfiles = simdetails.afiles{isim};
        elseif itype == 2
            outfiles = simdetails.hfiles{isim};
        elseif itype == 3
            outfiles = simdetails.sfiles{isim};
        end
    
        % Categorize variables that exist in the raw data files
        % Note: production variables are treated specially because of the
        % way they need to be integrated/summed based on model version.

        I = ncinfo(outfiles(1).name);
        
        if ismember('allfluxes', Vars.lev1{itype}) % special flag for all flux variables
            
            isflux = regexpfound({I.Variables.Name}, '(\w*)_(\w*)_(\w*)');
 
            Vars.lev1{itype} = union(...
                setdiff(Vars.lev1{itype}, 'allfluxes', 'stable'), ...
                {I.Variables(isflux).Name}, 'stable');
        end
            
        haslev1 = ismember(Vars.lev1{itype}, {I.Variables.Name});

        prodvars = Vars.lev1{itype}(haslev1 & ...
                            startsWith(Vars.lev1{itype}, 'prod') & ...
                            ~contains(Vars.lev1{itype}, 'Ben') & ...
                            ~contains(Vars.lev1{itype}, 'Ice'));
                        
        Vars.tot{itype} = setdiff(Vars.tot{itype}, prodvars, 'stable');          
                        
        % Group into 5-year intervals

        tedge = datetime(1970:5:2100,1,1);

        ngrp = length(tedge)-1;

        % Query time bounds for all files

        fprintf(fid, 'Gathering file info (%s, %s, %d/%d)...\n', filetype{itype}, simdetails.base{isim}, isim, nsim);

        filetlim = ncdatelim({outfiles.name}', 'ocean_time');

        % Now process!

        for ii = 1:ngrp

            try
                basegrp = sprintf('%s_%d-%d', simdetails.base{isim}, year(tedge(ii)), year(tedge(ii+1)-days(1)));

                lev1folder = fullfile(simdetails.ppfol{isim}, 'Level1');
                lev2folder = fullfile(simdetails.ppfol{isim}, 'Level2');

                % Check for existing files

                [l1v,l2v] = existinglev1and2(simdetails.ppfol{isim}, ...
                                             basegrp, filetype{itype});

                % Figure out which station, history, and average files fall in the
                % specified range 

                frac = fracoverlap(filetlim, tedge([ii ii+1]));
                if ~any(frac>0)
                    continue; % Skip remaining calcs if no files in this time reange
                end

                grpavg = {outfiles(frac>0).name}';

                % Call level 1 postprocessing for indicated variables
                % (Note: it would be a little more efficient to do all file
                % types at once, but this keeps the script more flexible).

                l1missing = setdiff(Vars.lev1{itype}(haslev1), l1v, 'stable');

                if ~isempty(l1missing)
                    fprintf(fid, '  Level 1 processing, %s (%d/%d), group %d/%d: %s\n', filetype{itype}, isim, nsim, ii, ngrp, basegrp);
                    if dryrunflag
                        fprintf(fid, '   Will add: ');
                        fprintf(fid, '%s, ', l1missing{:});
                        fprintf(fid, '\n');
                    else
                        if itype == 1
                            roms_level1(lev1folder, basegrp, grdfile, ...
                                'avg', grpavg, ...
                                'his', [], ...
                                'sta', [], ...
                                'variables', l1missing, ...
                                'tbound', [tedge(ii) tedge(ii+1)], ...
                                'verbose', vflag);
                        elseif itype == 2
                            roms_level1(lev1folder, basegrp, grdfile, ...
                                'avg', [], ...
                                'his', grpavg, ...
                                'sta', [], ...
                                'variables', l1missing, ...
                                'tbound', [tedge(ii) tedge(ii+1)], ...
                                'verbose', vflag);
                        elseif itype == 3
                            roms_level1(lev1folder, basegrp, grdfile, ...
                                'avg', [], ...
                                'his', [], ...
                                'sta', grpavg, ...
                                'variables', l1missing, ...
                                'tbound', [tedge(ii) tedge(ii+1)], ...
                                'verbose', vflag);
                        end
                    end
                end

                % Level 2: Add pre-calculated depths variables

                if Vars.addz(itype)
                    hasz = ismember({'z_psi', 'z_psi_w', 'z_rho', 'z_u', 'z_u_w', ...
                                     'z_v', 'z_v_w', 'z_w'}, l2v);
                    hasz = (itype<3 && all(hasz)) || (itype==3 && all(hasz([3 8])));
                    if ~hasz
                        fprintf(fid, '  Level 2 depths, %s (%d/%d), group %d/%d: %s\n', filetype{itype}, isim, nsim, ii, ngrp, basegrp);
                        if dryrunflag
                            fprintf(fid, '   Will add depth variables\n');
                        else
                            roms_level2_depths(grdfile, lev1folder, lev2folder, basegrp, ...
                                'verbose', vflag, 'ftype', filetype{itype});
                        end
                    end
                end

                % Level 2: Rotate u and v and move to common grid
                
                if Vars.rotateuv(itype)

                    hasuv = ismember({'uEast', 'vNorth'}, l2v);
                    if ~all(hasuv)
                        fprintf(fid, '  Level 2 u/v rotation, %s (%d/%d), group %d/%d: %s\n', filetype{itype}, isim, nsim, ii, ngrp, basegrp);

                        if dryrunflag
                            fprintf(fid, '   Will add uEast and vNorth\n');
                        else
                            roms_level2_rotateuv(grdfile, lev1folder, lev2folder, basegrp, ...
                                'verbose', vflag, 'ftype', 'average');
                        end
                    end
                end

                % Level 2: extract bottom/surface values

                srfbot = union(Vars.srf{itype}, Vars.bot{itype});
                isdone =  ismember([cellfun(@(x) [x '_surface5m'], srfbot', 'uni', 0), ...
                                    cellfun(@(x) [x '_bottom5m' ], srfbot', 'uni', 0)], ...
                                    l2v);
                canbedone = ismember(srfbot', l1v) | ismember(srfbot', l2v); % parent file
                isneeded = [ismember(srfbot', Vars.srf{itype}) ismember(srfbot', Vars.bot{itype})] & ...
                            ~isdone & canbedone;         
                
                if any(isneeded)  
                    fprintf(fid, '  Level 2 surface/bottom averaging, %s (%d/%d), group %d/%d: %s\n', filetype{itype}, isim, nsim, ii, ngrp, basegrp);

                    if dryrunflag
                        srfstr = sprintf('%s,', srfbot{isneeded(:,1)});
                        botstr = sprintf('%s,', srfbot{isneeded(:,2)});
                        fprintf(fid, '   Will add surface for: %s\n', srfstr(1:end-1));
                        fprintf(fid, '   Will add bottom for: %s\n', botstr(1:end-1));
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
                                        'verbose', vflag, ...
                                        'location', loc);
                                end
                            end
                        end 
                    end
                end

                % Level 2: depth-integrated values

                canbedone = ismember(Vars.tot{itype}, l1v) | ismember(Vars.tot{itype}, l2v); % parent file
                isdone =  ismember(cellfun(@(x) [x '_integrated'], Vars.tot{itype}, 'uni', 0), l2v) | ~canbedone;
                
                if ~all(isdone)
                    fprintf(fid, '  Level 2 depth integration (%d/%d), group %d/%d: %s\n', isim, nsim, ii, ngrp, basegrp);
                    if dryrunflag
                        totstr = sprintf('%s,', Vars.tot{itype}{~isdone});
                        fprintf(fid, '   Will add integrated for: %s\n', totstr(1:end-1));
                    else
                        for iv = 1:length(Vars.tot{itype})
                            if ~isdone(iv)
                                
                                par1 = fullfile(lev1folder, sprintf('%s_%s_%s.nc', basegrp, filetype{itype}, Vars.tot{itype}{iv}));
                                par2 = fullfile(lev2folder, sprintf('%s_%s_%s.nc', basegrp, filetype{itype}, Vars.tot{itype}{iv}));
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
                                        Vars.tot{itype}{iv}, ...
                                        'verbose', vflag);
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
                
                if Vars.integrateprod(itype)

                    ish16 = contains(simdetails.base{isim}, 'H16');
                    isdone =  ismember(cellfun(@(x) [x '_integrated'], prodvars, 'uni', 0), l2v);
                    
                    if ~all(isdone)
                        fprintf(fid, '  Level 2 prod. integration, %s (%d/%d), group %d/%d: %s\n', filetype{itype}, isim, nsim, ii, ngrp, basegrp);
                        if dryrunflag
                            totstr = sprintf('%s,', prodvars{~isdone});
                            fprintf(fid, '   Will add integrated for: %s\n', totstr(1:end-1));
                        else
                            for iv = 1:length(prodvars)
                                if ~isdone(iv)
                                    roms_level2_depthintegrate(...
                                        fullfile(lev1folder, sprintf('%s_%s_%s.nc', basegrp, filetype{itype}, prodvars{iv})), ...
                                        lev2folder, ...
                                        sprintf('%s_%s', basegrp, filetype{itype}), ...
                                        fullfile(lev2folder, sprintf('%s_%s', basegrp, filetype{itype})), ...
                                        prodvars{iv}, ...
                                        'verbose', vflag, ...
                                        'simplesum', ~ish16);
                                end
                            end
                        end
                    end
                end

            catch ME
                fprintf(fid, '%s: An error occurred processing %s, sim %d, group %d (%s)\n', datestr(now), filetype{itype}, isim, ii, basegrp);
                rethrow(ME);
            end
        end
    end

    
end

if dryrunflag && uselog
    fclose(fid);
end

%-----------------------
% Subfunctions
%-----------------------

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
