%% What variables are available?
%
% This script does most of the work to build the variable list spreadsheet
% Bering10K_simulation_variables_YYYYMMDD.xlsx, saved to
% ~/Documents/Research/ManuscriptsReviews/2020_ROMS_use_docs.  After
% running this script, the following manual aesthetic tweaks to the
% spreadsheet are usually made:  
%
%   1) Bold row 3 (table headings)
%   2) Rotate simulation names to vertical
%   3) Resize simulation columns to a width of 4
%   4) Merge average, history, and station cells above simulation names,
%      and shade (light blue, light orange, light yellow)
%   5) Color simulation names blue, orange, and gold
%
% The final spreadsheet documents which variables are available within each
% simulation that has been postprocessed.

onmox = contains(moxdir, 'gscratch'); % read on mox, build file on local

%--------------------
% Variable summaries
%--------------------

varfile = sprintf('varsummary_%s.mat', datestr(datetime('today'), 'yyyymmdd'));

if onmox
    V = processed_variable_summary(true);
    R = raw_variable_summary;
 
    % Check that the same simulations are in the two

    test = isequal(sort(V.sname), sort(R.sname));
    if ~test
        error('Sim mismatch... check raw_variables_summary table?');
    end

    %---------------------
    % Build combined table
    %---------------------

    sims = sort(V.sname);
    vars = setdiff(union(V.vname, R.vname), 'constants');

    nv = length(vars);
    ns = length(sims);

    % In raw files?  Fill in variable details from raw where possible

    Vtbl = table;
    Vtbl.name = vars;
    Vtbl.lname = cell(nv,1);
    Vtbl.units = cell(nv,1);
    Vtbl.grid_ah = cell(nv,1);
    Vtbl.grid_s = cell(nv,1);

    Vtbl.status = nan(nv,ns,3);

    [vtf,vloc] = ismember(R.vname, Vtbl.name);
    [~,sloc] = ismember(R.sname, sims);
    isinraw = false(nv,ns,3);
    isinraw(vloc(vtf),sloc,:) = R.varsaved;

    Vtbl.lname(vloc(vtf))   = R.lname(vtf);
    Vtbl.units(vloc(vtf))   = R.units(vtf);
    Vtbl.grid_ah(vloc(vtf)) = R.grid_ah(vtf);
    Vtbl.grid_s(vloc(vtf))  = R.grid_s(vtf);

    Vtbl.status(isinraw) = 0;

    % In processed files?  Fill in variable details where needed

    [vtf,vloc] = ismember(V.vname, Vtbl.name);
    [~,sloc] = ismember(V.sname, sims);

    hasfile = ~cellfun(@isempty, V.files);
    levfile = V.files(hasfile);
    levfile = regexp(levfile, 'Level(\d)', 'tokens', 'once');
    levfile = cellfun(@str2double, cat(1, levfile{:}));
    lev = nan(size(V.files));
    lev(hasfile) = levfile;
    lev = permute(max(lev,[],2), [3 1 4 2]);

    levinprocessed = nan(nv,ns,3);
    levinprocessed(vloc(vtf),sloc,:) = lev(vtf,:,:);

    mask = levinprocessed>0;
    Vtbl.status(mask) = levinprocessed(mask);

    isemp = cellfun(@isempty, Vtbl.lname);
    [~,loc] = ismember(Vtbl.name(isemp), V.vname);
    Vtbl.lname(isemp) = V.lname(loc);
    Vtbl.units(isemp) = V.units(loc);

    gtype = V.coord(loc);
    dims = V.ndim(loc);
    for ii = 1:length(gtype)
        if ~isempty(gtype{ii})
            gtype{ii} = sprintf('%s%dd', gtype{ii}(1), dims(ii)-1);
        end
    end
    Vtbl.grid_ah(isemp) = gtype;

    for iv = 1:nv
        if isempty(Vtbl.grid_ah{iv})
            Vtbl.grid_ah{iv} = '';
        end
        if isempty(Vtbl.grid_s{iv})
            Vtbl.grid_s{iv} = '';
        end
    end

    % ... constants

    isinlev1c = false(nv,ns,3); % in constants file

    for is = 1:ns
        isc = strcmp(V.vname, 'constants');
        for it = 1:3
            cfiles = V.files(is,:,isc,it);
            mask = cellfun(@isempty, cfiles);
            if ~all(mask)
                cfiles = cfiles(~mask);
                I = ncinfo(cfiles{1});
                cvars = {I.Variables.Name};
                isinlev1c(:,is,it) = ismember(Vtbl.name, cvars);
            end
        end
    end

    mask = isinlev1c & Vtbl.status == 0;
    Vtbl.status(mask) = 1;

    % Compare grid variables and classify variables as well as possible
    % (Note: anything 2d and in station files only won't be able to be
    % classified as rho, u, v, or psi).

    vgrid1 = Vtbl.grid_ah;
    gridtranslate = {'2d'   '*2d'
                     'r'    'r3d'
                     'w'    'w3d'};
    [tf,loc] = ismember(Vtbl.grid_s, gridtranslate(:,1));
    vgrid2 = Vtbl.grid_s;
    vgrid2(tf) = gridtranslate(loc(tf),2);

    isemp = cellfun(@isempty, [vgrid1 vgrid2]);
    use2 = isemp(:,1) & ~isemp(:,2);
    vgrid1(use2) = vgrid2(use2);

    Vtbl.grid = vgrid1;
    Vtbl = removevars(Vtbl, {'grid_ah','grid_s'});

    % Sort

    Vtbl = sortrows(Vtbl, {'grid', 'name'});
    
    % Save

    save(varfile, 'Vtbl', '-append');
    return
else
    Vtbl = load(varfile);
    Vtbl = Vtbl.Vtbl;
end

%% Test plot

if ~onmox
    h = spycolor(reshape(Vtbl.status, nv, []), NaN);
    set(h(1), 'sizedata', 5);
    gridxy(ns.*(1:3)+0.5);
    cmap = cptcmap('Set1_09');
    cmap = cmap(1:4,:);
    set(gca, 'clim', [-0.5 3.5], 'colormap', cmap);
end


%% 

%---------------------
% Write to file
%---------------------

Main = removevars(Vtbl, 'status');

Status = reshape(Vtbl.status, nv, []);
Status = array2table(Status);

% The template sets up the header info, and removes unnecessary extra
% sheets.

outdir = '~/Documents/Research/ManuscriptsReviews/2020_ROMS_use_docs';
outname = fullfile(outdir, sprintf('Bering10K_simulation_variables_%s.xlsx', datestr(datetime('today'), 'yyyymmdd')));
templatename = fullfile(outdir, 'Bering10K_simulation_variables_template.xlsx');

copyfile(templatename, outname);

% Add update date

writematrix(datetime('now'), outname, 'Sheet', 'Variables', 'Range', 'B2');

% Add main table

writetable([Main Status], outname, 'Sheet', 'Variables', 'Range', 'A3');

% Simulation header info

header = [...
    repmat({''}, 1, ns*3)
    repmat(sims', 1, 3)];
header(1,(0:2)*ns+1) = V.tname';
header = cell2table(header);

writetable(header, outname, 'Sheet', 'Variables', 'Range', 'E2', ...
    'WriteVariableNames', false);


return



% Simulations: base name, postprocessing parent folder, list of average files
% Note: col 2 (ppfol) is legacy, now overwritten by roms_for_public
% location.
% Note 2: Always make sure this list includes all runs from
% postprocess_aclim

simdetails = {
'B10K-K20P19_CMIP6_cesm_historical'  './Postprocessed'                       dirfull(fullfile('.', 'cesm_historical',  'Out', '*avg*.nc'))
'B10K-K20P19_CMIP6_cesm_ssp126'      './Postprocessed'                       dirfull(fullfile('.', 'cesm_ssp126',      'Out', '*avg*.nc'))
'B10K-K20P19_CMIP6_cesm_ssp585'      './Postprocessed'                       dirfull(fullfile('.', 'cesm_ssp585',      'Out', '*avg*.nc'))
'B10K-K20P19_CMIP6_gfdl_historical'  './Postprocessed'                       dirfull(fullfile('.', 'gfdl_historical',  'Out', '*avg*.nc'))
'B10K-K20P19_CMIP6_gfdl_ssp126'      './Postprocessed'                       dirfull(fullfile('.', 'gfdl_ssp126',      'Out', '*avg*.nc'))
'B10K-K20P19_CMIP6_gfdl_ssp585'      './Postprocessed'                       dirfull(fullfile('.', 'gfdl_ssp585',      'Out', '*avg*.nc'))
'B10K-K20P19_CMIP6_miroc_historical' './Postprocessed'                       dirfull(fullfile('.', 'miroc_historical', 'Out', '*avg*.nc'))
'B10K-K20P19_CMIP6_miroc_ssp126'     './Postprocessed'                       dirfull(fullfile('.', 'miroc_ssp126',     'Out', '*avg*.nc'))
'B10K-K20P19_CMIP6_miroc_ssp585'     './Postprocessed'                       dirfull(fullfile('.', 'miroc_ssp585',     'Out', '*avg*.nc'))
'B10K-K20_CORECFS'                   '../../hindcasts/npz_201904_delta30/',  dirfull('../../hindcasts/npz_201904_delta30/Out/*avg*.nc')
'B10K-H16_CORECFS'                   '../../hindcasts/npz_201904_aclim/',    dirfull('../../hindcasts/npz_201904_aclim/Out/*avg*.nc')
'B10K-K20P19_CMIP5_CESM_rcp45'       '../../cmip5-avg-carbon/Postprocessed', dirfull('../../cmip5-avg-carbon/CESM_rcp45/*avg*.nc')
'B10K-K20P19_CMIP5_CESM_rcp85'       '../../cmip5-avg-carbon/Postprocessed', dirfull('../../cmip5-avg-carbon/CESM_rcp85/*avg*.nc')
'B10K-K20P19_CMIP5_GFDL_rcp45'       '../../cmip5-avg-carbon/Postprocessed', dirfull('../../cmip5-avg-carbon/GFDL_rcp45/*avg*.nc')
'B10K-K20P19_CMIP5_GFDL_rcp85'       '../../cmip5-avg-carbon/Postprocessed', dirfull('../../cmip5-avg-carbon/GFDL_rcp85/*avg*.nc')
'B10K-K20P19_CMIP5_MIROC_rcp45'      '../../cmip5-avg-carbon/Postprocessed', dirfull('../../cmip5-avg-carbon/MIROC_rcp45/*avg*.nc')
'B10K-K20P19_CMIP5_MIROC_rcp85'      '../../cmip5-avg-carbon/Postprocessed', dirfull('../../cmip5-avg-carbon/MIROC_rcp85/*avg*.nc')
'B10K-H16_CMIP5_CESM_BIO_rcp85'      '../../cmip5-avg/Postprocessed',        dirfull('../../cmip5-avg/CESM_BIO_rcp85/*avg*.nc')
'B10K-H16_CMIP5_CESM_rcp45'          '../../cmip5-avg/Postprocessed',        dirfull('../../cmip5-avg/CESM_rcp45/*avg*.nc')
'B10K-H16_CMIP5_CESM_rcp85'          '../../cmip5-avg/Postprocessed',        dirfull('../../cmip5-avg/CESM_rcp85/*avg*.nc')
'B10K-H16_CMIP5_GFDL_BIO_rcp85'      '../../cmip5-avg/Postprocessed',        dirfull('../../cmip5-avg/GFDL_BIO_rcp85/*avg*.nc')
'B10K-H16_CMIP5_GFDL_rcp45'          '../../cmip5-avg/Postprocessed',        dirfull('../../cmip5-avg/GFDL_rcp45/*avg*.nc')
'B10K-H16_CMIP5_GFDL_rcp85'          '../../cmip5-avg/Postprocessed',        dirfull('../../cmip5-avg/GFDL_rcp85/*avg*.nc')
'B10K-H16_CMIP5_MIROC_rcp45'         '../../cmip5-avg/Postprocessed',        dirfull('../../cmip5-avg/MIROC_rcp45/*avg*.nc')
'B10K-H16_CMIP5_MIROC_rcp85'         '../../cmip5-avg/Postprocessed',        dirfull('../../cmip5-avg/MIROC_rcp85/*avg*.nc')
'B10K-K20nobio_CFS'                  ''                                      dirfull('../../hindcasts/phys_201910_30_singlesource/Out/physhccfs*_avg_*.nc')
'B10K-K20nobio_CORE'                 ''                                      dirfull('../../hindcasts/phys_201910_30_singlesource/Out/physhccore*_avg_*.nc')
};

simdetails = cell2table(simdetails, 'VariableNames', {'base','ppfol','files'});
nsim = height(simdetails);

simdetails.ppfol = fullfile(moxdir, 'roms_for_public', simdetails.base);

% Scan the postprocessed files

V = aclim_variable_summary(false, true);

test = isequal(sort(V.sname), sort(simdetails.base));
if ~test
    error('Sim mismatch... check simdetails table vs aclim_variable_summary');
end

% Scan the raw files

for is = 1:height(simdetails)
    S(is) = simvariableinfo({simdetails.files{is}.name}, 'average');
end
[tf,loc] = ismember(V.sname, simdetails.base);
S = S(loc);

refidx = find(strcmp(V.sname, 'B10K-K20_CORECFS')); % for long name and unit clashes
for is = 1:length(S)
    S(is).priority = repmat(is==refidx, length(S(is).name), 1);
end

% Level 0: what's in the raw files?

lev1 = arrayfun(@struct2table, S, 'uni', 0);
lev1 = cat(1, lev1{:});

lev1 = sortrows(lev1, 'priority', 'descend');

[unq, iu, iunq] = unique(lev1(:,{'name','grid'}), 'rows', 'first');
lev1 = removevars(lev1(iu,:), 'priority');

lev1 = sortrows(lev1, {'grid', 'name'});

nlev1 = height(lev1);
nsim = length(V.sname);

isinraw = false(nlev1,nsim);
for is = 1:nsim
    isinraw(:,is) = ismember(lev1.name, S(is).name);
end

% Level 1: What variables have been processed so far?

isinlev1c = false(nlev1,nsim); % in constants file
isinlev1  = false(nlev1,nsim); % in their own file

mask = ~cellfun(@isempty,V.files);
isc = strcmp(V.vname, 'constants');

vmask = permute(any(mask,2), [3 1 2]);

for is = 1:nsim
    cfiles = V.files(is,:,isc);
    cfiles = cfiles(~cellfun(@isempty,cfiles));
    I = ncinfo(cfiles{1});
    cvar = {I.Variables.Name};
    isinlev1c(:,is) = ismember(lev1.name, cvar);
    
    isinlev1(:,is) = ismember(lev1.name, V.vname(vmask(:,is)));
end

% Level 2: What's been added to that?

lev2 = struct;
[lev2.name, vidx] = setdiff(V.vname, [lev1.name; 'constants']);
lev2.grid = cell(size(lev2.name));
lev2.lname = V.lname(vidx);
lev2.units = V.units(vidx);

vargroup = {...
            'r2d' 'xi_rho,eta_rho,ocean_time,'      
            'r3d' 'xi_rho,eta_rho,s_rho,ocean_time,'
            'w3d' 'xi_rho,eta_rho,s_w,ocean_time,'  
            'u2d' 'xi_u,eta_u,ocean_time,'          
            'u3d' 'xi_u,eta_u,s_rho,ocean_time,'    
            'v2d' 'xi_v,eta_v,ocean_time,'          
            'v3d' 'xi_v,eta_v,s_rho,ocean_time,'
            'p3d' 'xi_psi,eta_psi,s_rho,ocean_time'
            'p3dw' 'xi_psi,eta_psi,s_w,ocean_time'
            'u3dw' 'xi_u,eta_u,s_w,ocean_time'
            'v3dw' 'xi_v,eta_v,s_w,ocean_time'
            };
for iv = 1:size(vargroup,1)
    tmp = regexp(vargroup{iv,2}, ',', 'split');
    tmp = sort(tmp(~cellfun(@isempty,tmp)));
    vargroup{iv,2} = sprintf('%s,', tmp{:});
end

mask = cellfun(@isempty,V.files);
dstr = cell(size(lev2.name));
[~,loc] = ismember(lev2.name, V.vname);
for ii = 1:length(loc)
    tmp = V.files(:,:,loc(ii));
    tmp = tmp(~mask(:,:,loc(ii)));
    I = ncinfo(tmp{1});
    dname = sort({I.Dimensions.Name});
    dstr{ii} = sprintf('%s,', dname{:});
end
[tf,gloc] = ismember(dstr, vargroup(:,2));
lev2.grid = vargroup(gloc,1);

lev2 = struct2table(lev2);

nlev2 = height(lev2);
isinlev2 = false(nlev2,nsim);
for is = 1:nsim
    isinlev2(:,is) = ismember(lev2.name, V.vname(vmask(:,is)));
end

% Combine

status = nan(size(isinraw));
status(isinraw & ~(isinlev1 | isinlev1c)) = 0;
status(isinlev1c) = 1;
status(isinlev1) = 1;

s2 = nan(size(isinlev2));
s2(isinlev2) = 2;

status = [status; s2];

Status = [[lev1; lev2] ...
          array2table(status, 'variablenames', strrep(V.sname, 'B10K-',''))];
      
       
%--------------
% Write to file
%--------------

% The template sets up the header info, and removes unnecessary extra
% sheets.

outdir = '~/Documents/Research/ManuscriptsReviews/2020_ROMS_use_docs';
outname = fullfile(outdir, sprintf('Bering10K_simulation_variables_%s.xlsx', datestr(datetime('today'), 'yyyymmdd')));
templatename = fullfile(outdir, 'Bering10K_simulation_variables_template.xlsx');

copyfile(templatename, outname);

% Add update date

writematrix(datetime('now'), outname, 'Sheet', 'Variables', 'Range', 'B2');

% Add main table

writetable(Status, outname, 'Sheet', 'Variables', 'Range', 'A3');
