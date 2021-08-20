function S = simvariableinfo(fname, ftype)
%SIMVARIABLEINFO Extract list of variables in set of output files
%
% S = simvariableinfo(fname, ftype)
%
% This function list the variables in a set of output files, along with the
% grid each variable falls on.  This more or less replicates the first bit
% of the roms_level1.m post-procssing script.
%
% Input variables:
%
%   fname:  cell array of strings, path names to a set of average, history,
%           or station files from a ROMS simulation
%
%   ftype:  type of simulation ('average', 'history', or 'station')
%
% Output variables:
%
%   S:      structure with the following fields:
%
%           name:   cell array of strings, name of output variables found
%                   in the file
%
%           grid:   cell array of strings, grid on which each variable
%                   falls. 'constant' indicates a non-time-varying variable
%                   of any sort.  History and average files can hold
%                   3d variables on 5 grids ([r]ho, w, u, v,
%                   [p]si) and 2d variables on rho, u, or v); station file
%                   variables can fall onto [r]ho, w, or 2d (i.e. no extra
%                   dimensions beyond station and time).
%
%           lname:  cell array of strings, long name of variables
%
%           units:  cell array of strings, variable units
               
% Copyright 2020 Kelly Kearney

fname = fname(:);
    
% First, group variables based on grid they use

cnt = 1;
I = ncinfo(fname{cnt});
while (I.Dimensions(strcmp({I.Dimensions.Name},'ocean_time')).Length == 0) % make sure we have some time steps in the reference file
    cnt = cnt+1;
    I = ncinfo(fname{cnt});
end    

vname = {I.Variables.Name}';
dname = cell(size(vname));
for iv = 1:length(vname)
    if isempty(I.Variables(iv).Dimensions)
        dname{iv} = '';
    else
        dname{iv} = sprintf('%s,', I.Variables(iv).Dimensions.Name);
    end
end

[dgrp, didx, id] = unique(dname);
Drep = arrayfun(@(X) X.Dimensions, I.Variables(didx), 'uni', 0);

% Groups of variables: constants, time, and everything else

isconstant = ~cellfun(@(x) ~isempty(x) && ismember('ocean_time', {x.Name}), Drep);
vars_constant = vname(ismember(id, find(isconstant)));

ist = strcmp(vname, 'ocean_time');

switch ftype
    case {'average', 'history'}
        vargroup = {...
            'r2d' 'xi_rho,eta_rho,ocean_time,'      
            'r3d' 'xi_rho,eta_rho,s_rho,ocean_time,'
            'w3d' 'xi_rho,eta_rho,s_w,ocean_time,'  
            'u2d' 'xi_u,eta_u,ocean_time,'          
            'u3d' 'xi_u,eta_u,s_rho,ocean_time,'    
            'v2d' 'xi_v,eta_v,ocean_time,'          
            'v3d' 'xi_v,eta_v,s_rho,ocean_time,'
            };
    case {'station'}
        vargroup = {...
            '2d'  'station,ocean_time,'      
            'r'   's_rho,station,ocean_time,'
            'w'   's_w,station,ocean_time,'
            };
end

% Determine which grid each variable is on

vgrp = cell(size(vname));
[tf,loc] = ismember(dgrp(id), vargroup(:,2));
vgrp(tf) = vargroup(loc(tf),1);
[vgrp{~tf}] = deal('constants');

istime = strcmp(vname, 'ocean_time');

S.name = vname(~istime);
S.grid = vgrp(~istime);

% Long names and units

[~,loc] = ismember(S.name, {I.Variables.Name});

[S.lname, S.units] = deal(cell(size(S.name)));
[S.lname{:}] = deal('');
[S.units{:}] = deal('');

atts = {I.Variables(loc).Attributes};
for ii = 1:length(atts)
    [tf,loc] = ismember({'long_name', 'units'}, {atts{ii}.Name});
    if tf(1); S.lname{ii} = atts{ii}(loc(1)).Value; end
    if tf(2); S.units{ii} = atts{ii}(loc(2)).Value; end
end
S.lname = strrep(S.lname, 'time-averaged ', '');

