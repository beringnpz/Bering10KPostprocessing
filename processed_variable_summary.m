function V = processed_variable_summary(getattrib)
%PROCESSED_VARIABLE_SUMMARY Returns details of postprocessed files
%
% V = processed_variable_summary(getattrib)
%
% This function summarizes the variables found in the Level 1-2
% post-processed files of the primary folder (roms_for_public).
%
% Input variables:
%
%   getattrib:  scalar logical, true to read certain details from the
%               files, false to skip this slower step (fields with an
%               asterisk below are skipped).  Default is false.  
%
% Output variables:
%
%   V:          structure  with the following fields:
%
%               sname:      cell array of strings, simulation names
%
%               yname:      cell array of strings, year-blocks (assumes
%                           similar processing across simulations) 
%
%               vname:      cell array of strings, output variables
%
%               tname:      cell array of strings, output types (average,
%                           history, station) 
%
%               files:      ns x ny x nv x nt cell array, files categorized
%                           by simulation, year block, variable, and type
%
%               size:       ns x ny x nv x nt, size of each file, in MB
%
%               vreffile:   nv x 1 cell array, reference file used to
%                           extract variable name and units corresponding
%                           to each variable
%
%               lname*:     nv x 1, variable long name
%
%               units*:     nv x 1, variable units
%
%               coord*:     nv x 1, variable coordinate type (rho, psi, u,
%                           v, w) 
%
%               ndim*:      nv x 1, number of dimensions variable has
%                           internally
%
%               trefidx*:   ns x 1, index of first year-block period
%                           available for each simulation
%
%               nlayer*:    ns x 1, number of layers in each simulation

% Copyright 2021 Kelly Kearney

if nargin < 1
    getattrib = false;
end

%--------------------
% List and categorize
% files
%--------------------

Ftbl = listprocessedfiles;

% Pull out level 1 and 2 files

isnative = Ftbl.level <= 2;
Ftbl = Ftbl(isnative,:);

% Classify files into simulation, type, year-block, and variable

[V.sname, ~, sidx] = unique(Ftbl.simulation);
[V.yname, ~, yidx] = unique(Ftbl.years);
[V.vname, ~, vidx] = unique(Ftbl.variable);
[V.tname, ~, tidx] = unique(Ftbl.type);

ns = max(sidx); % sims
ny = max(yidx); % year-blocks
nv = max(vidx); % variables
nt = max(tidx); % types

V.files = cell(ns,ny,nv,nt);
idx = sub2ind(size(V.files), sidx, yidx, vidx, tidx);

tmp = strfind(Ftbl.folder, 'roms_for_public');
shortfol = cellfun(@(a,b) a(b:end), Ftbl.folder, tmp, 'uni', 0);

V.files(idx) = fullfile(Ftbl.folder, Ftbl.name);

sref = 'B10K-K20_CORECFS';
[~,srefidx] = ismember(sref, V.sname);

% File size (MB)

V.size = nan(size(V.files));
V.size(idx) = Ftbl.bytes/1024/1024; % mb

%--------------------
% Pull attributes
%--------------------

% For attributes and units, we want to rely on the K20 hindcast, except
% where that sim is missing a variable.  Figure out reference file for each
% variable.  Also, prioritize average file, then station.

[~,tploc] = ismember({'average','station'}, V.tname);

hasvar = ~cellfun(@isempty, V.files);
V.vreffile = cell(nv,1);
for iv = 1:nv
    
    flag1 = hasvar(srefidx,:,iv,tploc(1)); % hindcast average
    flag2 = hasvar(srefidx,:,iv,tploc(2)); % hindcast station
    flag3 = hasvar(:,:,iv,tploc(1)); % any average
    flag4 = hasvar(:,:,iv,tploc(2)); % any station
    
    if any(flag1(:))
        tmp = V.files(srefidx,:,iv,tploc(1));
        idx = find(hasvar(srefidx,:,iv,tploc(1)), 1);
        V.vreffile{iv} = tmp{idx};
    elseif any(flag2(:))
        tmp = V.files(srefidx,:,iv,tploc(2));
        idx = find(hasvar(srefidx,:,iv,tploc(2)), 1);
        V.vreffile{iv} = tmp{idx};
    elseif any(flag3(:))
        tmp = V.files(:,:,iv,tploc(1));
        idx = find(hasvar(:,:,iv,tploc(1)), 1);
        V.vreffile{iv} = tmp{idx};
    elseif any(flag4(:))
        tmp = V.files(:,:,iv,tploc(2));
        idx = find(hasvar(:,:,iv,tploc(2)), 1);
        V.vreffile{iv} = tmp{idx};
    else
        error('No file found for variable %d (%s)?', iv, V.vname{iv});
    end    
end
    
if getattrib
    
    % Variable details
    
    V.lname = cell(nv,1);
    V.units = cell(nv,1);
    V.coord = cell(nv,1);
    V.ndim = zeros(nv,1);
    
    [V.units{:}] = deal('');
    [V.lname{:}] = deal('');
    [V.coord{:}] = deal('');
    
    for iv = 1:nv
        var = strrep(V.vname{iv}, '_integrated', '');
        var = strrep(var, '_surface5m', '');
        var = strrep(var, '_bottom5m', '');

        if ~strcmp(var, 'constants')  
            Vinfo = ncinfo(V.vreffile{iv}, var);
            
            vatt = struct2table(Vinfo.Attributes);
            [~,vloc] = ismember({'long_name', 'units'}, vatt.Name);
            if vloc(1)>0
                V.lname{iv} = vatt.Value{vloc(1)};
            end
            if vloc(2)>0
                V.units{iv} = vatt.Value{vloc(2)};
            end
            
            if contains(V.vreffile{iv}, 'average')
                V.ndim(iv) = length(Vinfo.Size);

                isxi = contains({Vinfo.Dimensions.Name}, 'xi');
                if any(isxi)
                    V.coord{iv} = strrep(Vinfo.Dimensions(isxi).Name, 'xi_', '');
                end
            else
%                 error('Need to write station-parsing coordinates bit still...');
                V.ndim(iv) = length(Vinfo.Size) + 1;
            end
        end
    end
    V.lname = strrep(V.lname, 'time-averaged ', ''); 

    % Simulation number of layers 

    [~,loc] = ismember('temp', V.vname);

    V.trefidx = zeros(ns,1);
    V.nlayer = zeros(ns,1);

    for is = 1:ns

        V.trefidx(is) = find(~cellfun(@isempty, V.files(is,:,loc)), 1);

        Vinfo = ncinfo(V.files{is,V.trefidx(is),loc}, 's_rho');

        V.nlayer(is) = Vinfo.Size;
    end

end






