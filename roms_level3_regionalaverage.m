function roms_level3_regionalaverage(sim, Opt)

arguments
    sim {mustBeTextScalar}
    Opt.grdfile {mustBeTextScalar} =fullfile(mounteddir('klone'), 'GR011846_reem/postprocessed/Bering10K_extended_grid.nc')
    Opt.verbose (1,1) {mustBeNumericOrLogical} =true
    Opt.mode {mustBeTextScalar, mustBeMember(Opt.mode, ["skip", "append", "dryrun"])} ="dryrun" 
    Opt.ppbase {mustBeTextScalar} =fullfile(mounteddir('klone'), 'GR011836_aclim/postprocessed')
    Opt.lev1 {mustBeText} =''
    Opt.lev2 {mustBeText} =''
    Opt.depthavg {mustBeText} =''
    Opt.addcoldpool (1,1) {mustBeNumericOrLogical} =true

end

% Extra input check

if ~exist(Opt.grdfile, 'file')
    error('Could not find indicated grid file: %s\n', Opt.grdfile);
end

filereg = fullfile(Opt.ppbase, sim, 'Level3', sprintf('ACLIMregion_%s.nc', sim));

appendflag = false;
if exist(filereg, 'file')
    switch Opt.mode
        case 'skip'
            fprintf('Skipping and exiting: (%s already exists)\n', filereg);
            return
        % case 'replace'
        %     fprintf('Deleting and replacing: (%s already exists)\n', filereg);
        %     delete(filereg);
        otherwise
            fprintf('Appending as necessary: (%s already exists)\n', filereg);
            appendflag = true;
    end
end

dryrun = strcmp(Opt.mode, 'dryrun');

% Default variables

% if isempty(Opt.lev1)
%     Opt.lev1 = {'Ben', 'DetBen', ...
%                 'Hsbl', ...
%                 'IceNH4', 'IceNO3', 'IcePhL', ...
%                 'aice', 'hice', ...
%                 'shflux', 'ssflux'};
% end
% if isempty(Opt.lev2)
%     Opt.lev2 = {'Cop_integrated'      , 'Cop_surface5m', ...       
%                 'EupO_integrated'     , 'EupO_surface5m', ...      
%                 'EupS_integrated'     , 'EupS_surface5m', ...      
%                 'Iron_bottom5m'       , 'Iron_integrated',     'Iron_surface5m', ...    
%                 'Fe_bottom5m'         , 'Fe_integrated',       'Fe_surface5m', ... 
%                 'Jel_integrated'      , 'Jel_surface5m', ...       
%                 'MZL_integrated'      , 'MZL_surface5m', ...       
%                 'NCaO_integrated'     , 'NCaO_surface5m', ...      
%                 'NCaS_integrated'     , 'NCaS_surface5m', ...      
%                 'NH4_bottom5m'        , 'NH4_integrated',      'NH4_surface5m', ...       
%                 'NO3_bottom5m'        , 'NO3_integrated',      'NO3_surface5m', ...       
%                 'PhL_integrated'      , 'PhL_surface5m' , ...      
%                 'PhS_integrated'      , 'PhS_surface5m' , ...      
%                 'prod_Cop_integrated' , ... 
%                 'prod_EupO_integrated', ... 
%                 'prod_EupS_integrated', ... 
%                 'prod_Eup_integrated' , ... 
%                 'prod_Jel_integrated' , ... 
%                 'prod_MZL_integrated' , ... 
%                 'prod_NCaO_integrated', ... 
%                 'prod_NCaS_integrated', ... 
%                 'prod_NCa_integrated' , ... 
%                 'prod_PhL_integrated' , ... 
%                 'prod_PhS_integrated' , ... 
%                 'salt_surface5m'      , ... 
%                 'temp_bottom5m'       , 'temp_integrated',     'temp_surface5m', ...      
%                 'uEast_bottom5m'      , 'uEast_surface5m'  , ...   
%                 'vNorth_bottom5m'     , 'vNorth_surface5m', ...
%                 'alkalinity_bottom5m' , 'alkalinity_integrated'  , 'alkalinity_surface5m', ...
%                 'TIC_bottom5m'        , 'TIC_integrated'         , 'TIC_surface5m', ...
%                 'oxygen_bottom5m'     , 'oxygen_integrated'      , 'oxygen_surface5m', ...
%                 'pH_bottom5m'         , 'pH_integrated'          , 'pH_surface5m', ...
%                 'calcite_bottom5m'    , 'calcite_integrated'     , 'calcite_surface5m', ...
%                 'arag_bottom5m'       , 'arag_integrated'        , 'arag_surface5m'};
% end
% if isempty(Opt.depthavg)
%     Opt.depthavg = {'temp', 'Iron', 'NH4', 'NO3', 'TIC', 'alkalinity', 'oxygen', ...
%                     'pH', 'calcite', 'arag'}';
% end

%---------------------------
% Setup
%---------------------------

if Opt.verbose
    fprintf('Setup...\n');
end

% Build masks

G = buildmasks(Opt.grdfile);


% Grd = ncstruct(Opt.grdfile, 'h', 'surveystrata_comboeast', 'area_feast');
% [nxi, neta] = size(Grd.h);
% 
% isshelf = Grd.h <= 200;
% 
% [sidx, stratanum] = findgroups(Grd.surveystrata_comboeast(:));
% 
% gmask = arrayfun(@(x) Grd.surveystrata_comboeast == x, stratanum, 'uni', 0);
% gmask_shelf = arrayfun(@(x) Grd.surveystrata_comboeast == x & isshelf, stratanum, 'uni', 0);
% 
% issame = cellfun(@(x,y) isequal(x,y), gmask, gmask_shelf);
% 
% gmask = [gmask(stratanum>0); gmask_shelf(~issame & stratanum>0)];
% 
% regname = [compose('stratum%03d', stratanum(stratanum>0)); ...
%            compose('stratum%03db', stratanum(~issame & stratanum>0))]; % text label, b = shelf-cut version
% regnum = [stratanum(stratanum>0); -stratanum(~issame & stratanum>0)]; % numeric label, negative = shelf-cut version
% nreg = length(regnum);
% 
% gmask = cat(3, gmask{:});
% nmask = size(gmask,3);
% 
% % Set up horizontal subsetting (to minimize data read)
% 
% [G.Scs, xilim, etalim] = romsmask2scs(any(gmask,3));
% 
% G.mask   = gmask(         xilim(1):xilim(2), etalim(1):etalim(2),:);
% G.weight = Grd.area_feast(xilim(1):xilim(2), etalim(1):etalim(2));
% G.depth  = Grd.h(         xilim(1):xilim(2), etalim(1):etalim(2));
% 
% G.area = zeros(nmask,1);
% for ii = 1:size(gmask,3)
%     G.area(ii) = sum(gweight(gmask(:,:,ii)));
% end

%---------------------------
% Time and variable analysis
%---------------------------

if Opt.verbose
    fprintf('Checking time and variables in Level 1/2 data...\n');
end

% What variables have been processed?

[Vpp, fsample] = parsefiles(Opt.ppbase, sim);
Vpp2d = Vpp(strcmp(Vpp.dimtype, '2r'),:);
Vpp3d = Vpp(strcmp(Vpp.dimtype, '3r'),:);

nv = height(Vpp2d);

% What is the current simulation time extent?

[T.t, T.num] = ncdateread(fsample, 'ocean_time');
T.long_name = ncreadatt(fsample{1}, 'ocean_time', 'long_name');
T.units     = ncreadatt(fsample{1}, 'ocean_time', 'units');
T.calendar  = ncreadatt(fsample{1}, 'ocean_time', 'calendar');

[T.t,T.idx] = unique(T.t, 'last'); 
T.num = T.num(T.idx);

Iz = ncinfo(fsample{1}, 's_rho');

% Compare to existing file, if present

Vpp2d.action = zeros(nv, 1); % 1 = add, 2 = extend, 0 = nothing

if exist(filereg, 'file')

    Iold = ncinfo(filereg);
    told = ncdateread(filereg, 'ocean_time');

    [vnew, inew] = setdiff(Vpp2d.var, {Iold.Variables.Name});

    if isequal(told, T.t)
        Vpp2d.action(inew) = 1;
    elseif isequal(told, T.t(1:length(told))) % new times appended
 
        tnew = T.t(T.t>told(end));

        Vpp2d.action(inew) = 1; % Add missing variables ...
        Vpp2d.action(Vpp2d.action ~= 1) = 2; % ... and extend the rest

        fprintf('  File will be extended from %s to %s', tnew(1), tnew(end));
        if ~dryrun
            ncwrite(filereg, 'ocean_time', T.num);
        end
    else
        error('Times in existing file do not match Level 1/2 data; cannot modify existing file');
    end

else
    Vpp2d.action = ones(nv,1);    
end

%--------------------
% Calculate indices
%--------------------

% Create file if necessary

if ~exist(filereg, 'file')
    if dryrun
        fprintf('Dry run: new file will be created\n');
    else
        createregionalaveragefile(filereg, T, Iz.size);
        ncwrite(filereg, 'region', G.regnum);
        ncwrite(filereg, 'region_area', G.area);
        ncwrite(filereg, 'ocean_time', T.num);
    end
end

% Add new variables if necessary

for iv = 1:height(Vpp2d)
    if Vpp2d.action(iv) == 1 % add new variable

        Vatts = queryvariableinfo(Vpp2d(iv,:), Opt.ppbase, sim);

        if dryrun
            fprintf('Dry run: %s will be added to file\n', Vatts.var);
        else
            nccreate(filereg, Vatts.var, ...
                'Dimensions', {'region', 'ocean_time'});
            if ~isempty(Vatts.long)
                ncwriteatt(filereg, Vatts.var, 'long_name', Vatts.long);
            end
            if ~isempty(Vatts.unit)
                ncwriteatt(filereg, Vatts.var, 'units', Vatts.unit);
            end
        end

        % TODO: add cold pool, depth average variables

    end
end

% Add/extend data

if ~dryrun
    w = warning('off', 'CDT:ncstructExtraScsField');
    
    fprintf('Adding data: \n')
    for iv = 1:height(Vpp2d)
    
        fprintf('  %s\n', Vpp2d.var{iv});
        adddata(Vpp2d(iv,:), Opt, sim, filereg, T, G)
    
    end
    
    warning(w);
end












end

%----------------------------
% Subfunction: variable table
%----------------------------

function G = buildmasks(grdfile)

    Grd = ncstruct(grdfile, 'h', 'surveystrata_comboeast', 'area_feast');
    [nxi, neta] = size(Grd.h);
    
    isshelf = Grd.h <= 200;
    
    [sidx, stratanum] = findgroups(Grd.surveystrata_comboeast(:));
    
    gmask = arrayfun(@(x) Grd.surveystrata_comboeast == x, stratanum, 'uni', 0);
    gmask_shelf = arrayfun(@(x) Grd.surveystrata_comboeast == x & isshelf, stratanum, 'uni', 0);
    
    issame = cellfun(@(x,y) isequal(x,y), gmask, gmask_shelf);
    
    gmask = [gmask(stratanum>0); gmask_shelf(~issame & stratanum>0)];
    
    G.regname = [compose('stratum%03d', stratanum(stratanum>0)); ...
               compose('stratum%03db', stratanum(~issame & stratanum>0))]; % text label, b = shelf-cut version
    G.regnum = [stratanum(stratanum>0); -stratanum(~issame & stratanum>0)]; % numeric label, negative = shelf-cut version
    
    gmask = cat(3, gmask{:});
    nmask = size(gmask,3);
    
    % Set up horizontal subsetting (to minimize data read)
    
    [G.Scs, xilim, etalim] = romsmask2scs(any(gmask,3));
    
    G.mask   = gmask(         xilim(1):xilim(2), etalim(1):etalim(2),:);
    G.weight = Grd.area_feast(xilim(1):xilim(2), etalim(1):etalim(2));
    G.depth  = Grd.h(         xilim(1):xilim(2), etalim(1):etalim(2));

    G.area = zeros(nmask,1);
    for ii = 1:size(G.mask,3)
        G.area(ii) = sum(G.weight(G.mask(:,:,ii)));
    end
end



function V = queryvariableinfo(V, ppbase, sim)
%QUERYVARIABLEINFO Query variable name in file, long name, and units

    V.internalshort = strrep(strrep(strrep(V.var, '_integrated', ''), '_surface5m', ''), '_bottom5m', '');

    F = dir(fullfile(ppbase, sim, sprintf('Level%d', V.level), sprintf('%s_*_average_%s.nc', sim, V.var)));
    fname = fullfile(F(1).folder, F(1).name);

    I = ncinfo(fname, V.internalshort);
    [tf,loc] = ismember({'long_name', 'units'}, {I.Attributes.Name});
    if tf(1)
        V.long = string(I.Attributes(loc(1)).Value);
    else
        V.long = "";
    end
    if tf(2)
        V.unit = string(I.Attributes(loc(2)).Value);
    else
        V.unit = "";
    end 
end


function Vtbl = variabletable(sim, ppbase, lev1, lev2, depthavg)

    % Get list of variable attributes from files
    
    ppfol = fullfile(ppbase, sim);

    cfile = dir(fullfile(ppbase, sim, 'Level1', sprintf('%s_*_average_constants.nc', sim)));
    cfile = fullfile(cfile(1).folder, cfile(1).name);
    Ic = ncinfo(cfile, 's_rho');
    nlayer = Ic.Size;
    
    for ii = length(lev1):-1:1
        V(ii).short = lev1{ii};
        V(ii).internalshort = V(ii).short;
        
        F = dir(fullfile(ppbase, sim, 'Level1', sprintf('%s_*_average_%s.nc', sim, V(ii).short)));
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
        
        F = dir(fullfile(ppbase, sim, 'Level2', sprintf('%s_*_average_%s.nc', sim, V(ii).short)));
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
    
    depthavg1 = cellfun(@(x) [x '_depthavg'], depthavg, 'uni', 0);
    depthavg2 = cellfun(@(x) [x '_integrated'], depthavg, 'uni', 0);
    
    [tf, loc] = ismember(depthavg2, Vtbl.short);
    unittmp = Vtbl.unit(loc(tf));
    unittmp = strrep(unittmp, ')*m', '');
    unittmp = regexprep(unittmp, '^\(', '');
    nnew = sum(tf);
    newdata = [depthavg1(tf), ...
               repmat({''}, nnew, 1), ...
               strrep(Vtbl.long(loc(tf)), 'integrated over depth', 'averaged over depth'), ...
               unittmp, ...
               num2cell(zeros(nnew,1))];
    
    Vtbl = [...
        Vtbl
        cell2table(newdata, 'variablenames', Vtbl.Properties.VariableNames)];
        
end

function [V, fsample] = parsefiles(ppbase, sim)
%PARSEFILES What variables exist in the Level 1/2 collection?

    fname = dir(fullfile(ppbase, sim, 'Level*', [sim '_*-*_average*.nc']));
    pattern = [sim '_\d\d\d\d-\d\d\d\d_average_(\w*)\.nc'];
    isavgpp = matches({fname.name}, regexpPattern(pattern));
    fname = fname(isavgpp); % just in case...

    vars = regexp({fname.name}, pattern, 'tokens', 'once');
    vars = cellfun(@(x) x{1}, vars, 'uni', 0);
    [vunq, ia, ic] = unique(vars, 'stable');

    [~,levpth] = fileparts({fname.folder});
    pplev = str2double(strrep(levpth(ia), 'Level', ''));

    fname = fullfile({fname.folder}, {fname.name});

    dimtype = cell(max(ic),1);
    for ii = 1:max(ic)
        I = ncinfo(fname{find(ic==ii,1)});
        if     isequal(sort({I.Dimensions.Name}), sort({'ocean_time', 'eta_rho', 'xi_rho'}))
            dimtype{ii} = '2r';
        elseif isequal(sort({I.Dimensions.Name}), sort({'ocean_time', 's_rho', 'eta_rho', 'xi_rho'}))
            dimtype{ii} = '3r';
        else
            dimtype{ii} = 'other';
        end
        

        % if all(ismember({'s_rho', 'ocean_time'}, {I.Dimensions.Name}))
        %     fsample = fname(ic == ii);
        %     return
        % end
    end

    V = table(string(vunq'), pplev', string(dimtype), 'variablenames', {'var','level','dimtype'});

    idx = find(strcmp(V.dimtype, '3r'), 1);
    fsample = fname(ic == idx);

end


function createregionalaveragefile(fname, T, nlayer) %, ppbase, sim, Vtbl, nlayer, nt)
%BUILDREGIONALAVERAGEFILE Builds file and add basic dimensions, attributes
 
    F = ncschema_init('classic');
    
    hisstr = sprintf('%s: %s', ...
                datetime('now', 'format', 'eee MMM dd HH:mm:SS yyyy'), ...
                'Regional ACLIM indices file created via roms_level3_regionalaverage.m');
            
    F = ncschema_addatts(F, ...
        'Simulation', Opt.sim, ...
        'history', hisstr, ...
        'Layers', nlayer);
    
    F = ncschema_adddims(F, ...
        'region',      nmask,     false, ...
        'ocean_time',  0, true);
    
    F = ncschema_addvars(F, ...
        'ocean_time', ...
        {'ocean_time'}, ...
        {'long_name', T.long_name, 'units', T.units, 'calendar', T.calendar}, ...
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
    
    % Create file and add coordinate variable data
    
    if ~exist(fname, 'file')
        ncwriteschema(fname, F);
    end
end


function adddata(V, Opt, sim, filereg, T, G)

    % Files
    
    avgfiles = dir(fullfile(Opt.ppbase, sim, sprintf('Level%d',V.level), sprintf('%s_*_average_%s.nc', sim, V.var)));
    avgfiles = fullfile({avgfiles.folder}, {avgfiles.name});
    
    % Set up output arrays

    varreg = ncread(filereg, V.var);

    cpflag = Opt.addcoldpool && strcmp(V.var, 'temp_bottom5m');
    if cpflag
        fracbelow0 = ncread(filereg, 'fracbelow0');
        fracbelow1 = ncread(filereg, 'fracbelow1');
        fracbelow2 = ncread(filereg, 'fracbelow2');
    end

    davgflag = endsWith(V.var, 'integrated') && ...
               ismember(strrep(V.var, '_integrated', ''), Opt.depthavg);
    if davgflag
        davgvarreg  = ncread(filereg, strrep(V.var, '_integrated', '_depthavg'));
    end

    % Where does data need to be added?

    ncid = netcdf.open(filereg,'NOWRITE');
    varid = netcdf.inqVarID(ncid,V.var);
    [~,fillValue] = netcdf.inqVarFill(ncid,varid);
    netcdf.close(ncid);
    ismissing = all(isnan(varreg) | varreg==fillValue, 1);
    ismissing = ismissing(:);

    tlim = ncdatelim(avgfiles, 'ocean_time');
    isin = any(T.t(ismissing)' >= tlim(:,1) & T.t(ismissing)' <= tlim(:,2),2);
    avgfiles = avgfiles(isin);

    % Read data from file

    Vatts = queryvariableinfo(V, Opt.ppbase, sim);
    nmask = size(G.mask,3);

    for ia = 1:length(avgfiles)

        if davgflag
            Tmp = ncstruct(avgfiles{ia}, Vatts.internalshort{1}, 'zeta', G.Scs);
        else
            Tmp = ncstruct(avgfiles{ia}, Vatts.internalshort{1}, G.Scs);
        end
        Tmp.t = ncdateread(avgfiles{ia}, 'ocean_time');

        % Regional averages

        [tf,loc] = ismember(Tmp.t,T.t);

        if ~all(tf)
            error('Out of range time found in %s', avgfiles{ia});
        end
        tf = tf & ismissing(loc(tf)); % only add missing values
        for ir = 1:nmask
            varreg(ir,loc(tf)) = local(Tmp.(Vatts.internalshort)(:,:,tf), G.mask(:,:,ir), 'weight', G.weight, 'omitnan');
        end

        % Calculate cold pool index fractions

        if cpflag
            for ir = 1:nmask
                fracbelow0(ir,loc(tf)) = local(G.weight.*double(Tmp.(Vatts.internalshort)(:,:,tf)<0), G.mask(:,:,ir), @nansum)/G.area(ir);
                fracbelow1(ir,loc(tf)) = local(G.weight.*double(Tmp.(Vatts.internalshort)(:,:,tf)<1), G.mask(:,:,ir), @nansum)/G.area(ir);
                fracbelow2(ir,loc(tf)) = local(G.weight.*double(Tmp.(Vatts.internalshort)(:,:,tf)<2), G.mask(:,:,ir), @nansum)/G.area(ir);
            end
        end
        
        % Depth-averaged
        
        if davgflag
            davgval = Tmp.(Vatts.internalshort{1})./(Tmp.zeta + G.depth);
            for ir = 1:nmask
                davgvarreg(ir,loc(tf)) = local(davgval(:,:,tf), G.mask(:,:,ir), 'weight', G.weight, 'omitnan');
            end
        end

    end
    
    % Write to file
        
    % vtmp = strrep(Vtbl.short{ii}, 'Iron_', 'Fe_'); % fix Fe/Iron mix
    
    ncwrite(filereg, Vatts.var, varreg);

    if cpflag
        ncwrite(filereg, 'fracbelow0', fracbelow0);
        ncwrite(filereg, 'fracbelow1', fracbelow1);
        ncwrite(filereg, 'fracbelow2', fracbelow2);
    end
    
    if davgflag
        ncwrite(filereg,  strrep(Vatts.var, 'integrated', 'depthavg'), davgvarreg);
    end
end



