function roms_level1(varargin)
%ROMS_LEVEL1 Level 1 ROMS output postprocessing
% 
% roms_level1(outfolder, outbase, gridfile, 'folder', folder)
% roms_level1(outfolder, outbase, gridfile, 'avg', avgfiles, ...
%                                           'his', hisfiles, ...
%                                           'sta', stafiles, ...)
% roms_level1(..., param, val, ...)
%
% This function reformats the output from a ROMS simulation to make it
% easier for the end user to work with.  It creates files that each contain
% a single primary output variable, along with its associated coordinate
% variables (including lat/lon from the grid file).  This is in contrast to
% native ROMS output, where all variables are in a single file, split over
% time.  The time concatenation also locates and removes any time
% replicates (an artifact of some restarts), favoring the post-restart
% version.
%
% This function currently applies to average, history, and station file
% output.  It should be applicable to any Bering10K ROMS simulation (and
% likely other domains as well, though it hasn't been tested for that yet).
%
% Input variables:
%
%   outfolder:  path to folder where postprocessed output will be saved
%
%   outbase:    base filename for newly-created files.  New files will be
%               named <outbase>_<ftype>_<variable>.nc, where <ftype> is
%               either 'average', 'history', or 'station' (based on ROMS
%               output data type), and <variable> indicates either the
%               indicated ROMS variable or 'constants' (the latter holds
%               all non-time-dependent variables).
%
%   gridfile:   path to ROMS grid file used for the simulation
%
% Optional input variables (passed as parameter/value pairs):
%
%   folder:     path to folder where ROMS raw output is saved. If this
%               option is used, all files matching the *avg*.nc, *his*.nc,
%               and *sta.nc patterns will be included in the averages,
%               history, and station datasets, respectively.  To include a
%               selected subset instead, use the 'avg', 'his', and 'sta'
%               input options and leave this empty. ['']
%
%   avg:        cell array of file names, indicating averages output files
%               to process.  If empty and no 'folder' input provided, no
%               averages data will be processed. [{}] 
%
%   his:        cell array of file names, indicating history output files
%               to process. If empty and no 'folder' input provided, no
%               history data will be processed. [{}]   
%
%   sta:        cell array of file names, indicating station output files
%               to process. If empty and no 'folder' input provided, no
%               stations data will be processed. [{}]  
%
%   variables:  cell array of variable names, indicating subset of
%               time-dependent variables to process.  In all cases,
%               additional time-independent variables and relevant
%               coordinate variables will also be added to the
%               postprocessed output.  If empty, all variables in the raw
%               output files will be processed. {[]} 
%
%   verbose:    scalar logical, true to print progress to screen [true]
%
%   tbound:     1 x 2 datetime array, indicating time interval to include
%               in the output files.  Times will be filtered inclusive of
%               the lower bound and exclusive of the upper bound.  If not
%               provided (or a [NaT NaT] array is passed as input), then
%               all unique times found in the input set of files will be
%               included.  Note that the bounds apply only to data that is
%               copied, not to global file attributes (all file's
%               attributes will be listed in the new attributes regardless
%               of whether their data are included in the final file). [NaT
%               NaT]    
%
%   constants:  logical scalar, true to create a file holding all
%               non-time-dependant variables from the original simulation
%               [true]

% Copyright 2020 Kelly Kearney

%----------------------
% Parse and check input
%----------------------

% Input parser

isfilecell = @(x) iscellstr(x) & all(cellfun(@(a) exist(a,'file'), x));

p = inputParser;
p.addRequired('outfolder',   @(x) validateattributes(x, {'char'}, {'scalartext'}));
p.addRequired('outbase',     @(x) validateattributes(x, {'char'}, {'scalartext'}));
p.addRequired('gridfile',    @(x) validateattributes(x, {'char'}, {'scalartext'}));
p.addParameter('folder', '', @(x) validateattributes(x, {'char'}, {'scalartext'}));
p.addParameter('avg', {}, @(x) isempty(x) || isfilecell(x));
p.addParameter('his', {}, @(x) isempty(x) || isfilecell(x));
p.addParameter('sta', {}, @(x) isempty(x) || isfilecell(x));
p.addParameter('dia', {}, @(x) isempty(x) || isfilecell(x));
p.addParameter('variables', {}, @(x) isempty(x) || iscellstr(x));
p.addParameter('verbose', true, @(x) validateattributes(x, {'logical'}, {'scalar'}));
p.addParameter('tbound', [NaT NaT], @(x) validateattributes(x, {'datetime'}, {'numel',2}));
p.addParameter('constants', true, @(x) validateattributes(x, {'logical'}, {'scalar'}));

p.parse(varargin{:});
Opt = p.Results;

% TODO: checks on tbound (increasing, all or no NaTs)

% Check that either a single folder or 3 separate lists of avg/his/sta
% input were provided, and set flags indicating which categories get
% processed.

hasfolder = ~isempty(Opt.folder) && isempty(Opt.avg) && isempty(Opt.his) && isempty(Opt.sta);

flag(1) = ~isempty(Opt.avg);
flag(2) = ~isempty(Opt.his);
flag(3) = ~isempty(Opt.sta);
flag(4) = ~isempty(Opt.dia);

hasfiles = isempty(Opt.folder) && any(flag);

if (hasfolder && hasfiles) || (~hasfolder && ~hasfiles)
    error('Must provide either folder or avg, his, sta, and dia inputs');
end

% If folder provided as input, build file list

if hasfolder
    F = dir(fullfile(Opt.folder, '*_avg_*.nc'));
    Opt.avg = arrayfun(@(X) fullfile(X.folder,X.name), F, 'uni', 0);
    F = dir(fullfile(Opt.folder, '*_his_*.nc'));
    Opt.his = arrayfun(@(X) fullfile(X.folder,X.name), F, 'uni', 0);
    F = dir(fullfile(Opt.folder, '*_sta_*.nc'));
    Opt.sta = arrayfun(@(X) fullfile(X.folder,X.name), F, 'uni', 0);
    F = dir(fullfile(Opt.folder, '*_dia_*.nc'));
    Opt.dia = arrayfun(@(X) fullfile(X.folder,X.name), F, 'uni', 0);
    
    flag = true(1,4);
end

% Output folder must exist previously (safety precaution...)

if ~exist(Opt.outfolder, 'dir')
    error('Output folder (%s) not found; must exist prior to running this script', Opt.outfolder);
end

% Check that grid file exists

if ~exist(Opt.gridfile, 'file')
    error('Grid file (%s) not found', Opt.gridfile);
end

% Check time bounds are properly formatted

if ~(all(isnat(Opt.tbound)) || (~any(isnat(Opt.tbound)) && Opt.tbound(2)>Opt.tbound(1)))
    error('tbound input must be datetimes, where tbound(2)>tbound(1)');
end

%----------------------
% Transfer data
%----------------------

% Main data transfer

if flag(1)
    if Opt.verbose
        fprintf('AVERAGE\n');
    end
    transferdata(Opt.avg, 'average', Opt.outfolder, Opt.outbase, Opt.gridfile, Opt.variables, Opt.tbound, Opt.verbose, Opt.constants);
end
if flag(2)
    if Opt.verbose
        fprintf('HISTORY\n');
    end
    transferdata(Opt.his, 'history', Opt.outfolder, Opt.outbase, Opt.gridfile, Opt.variables, Opt.tbound, Opt.verbose, Opt.constants);
end
if flag(3)
    if Opt.verbose
        fprintf('STATIONS\n');
    end
    transferdata(Opt.sta, 'station', Opt.outfolder, Opt.outbase, Opt.gridfile, Opt.variables, Opt.tbound, Opt.verbose, Opt.constants);
end
if flag(4)
    if Opt.verbose
        fprintf('DIAGNOSTICS\n');
    end
    transferdata(Opt.dia, 'diagnos', Opt.outfolder, Opt.outbase, Opt.gridfile, Opt.variables, Opt.tbound, Opt.verbose, Opt.constants);
end


%----------------------
% Transfer data: 
% subfunction
%----------------------

function transferdata(fname, ftype, outfolder, outbase, gridfile, variables, tbnd, vflag, cflag)    
    
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
        case {'average', 'history', 'diagnos'}
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
        otherwise
            warning('New file type? roms_level1 may not be set up to handle this yet.');
    end

    % Get attributes
    
    Att = simattributes(fname);
    
    % Figure out which time steps and/or files need to be dropped

    if vflag
        fprintf('Analyzing time values...');
    end

    Itmp = cellfun(@(x) ncinfo(x, 'ocean_time'), fname);
    isemp = [Itmp.Size] == 0;
    fname = fname(~isemp);
    nfile = length(fname);

    t = cellfun(@(x) ncread(x, 'ocean_time'), fname, 'uni', 0); 
    tall = cat(1, t{:}); % time in ROMS units
    tdate = ncdateread(fname{1}, 'ocean_time', tall); % datetimes

    if isnat(tbnd(1))
        [tunq, tidx] = unique(tall, 'last');
    else
        [tunq, tidx] = unique(tall, 'last');
        tdate = tdate(tidx);
        isin = tdate >= tbnd(1) & tdate < tbnd(2);
        tunq = tunq(isin);
        tidx = tidx(isin);
    end
    nt = length(tunq);
    
    % Calculate time indices determining how to read just the unique time
    % values, and where to place them in the new files

    indexfull = mat2cell(1:length(tall), 1, cellfun(@length, t)); % index in tall

    isin = cell(size(t));
    for ifl = 1:length(t)
        isin{ifl} = ismember(indexfull{ifl}, tidx);
    end

    nadd = cellfun(@sum, isin);
    sidx = [0; cumsum(nadd(1:end-1))]+1;

    if vflag
        fprintf('done\n');
    end

    % Determine which grid each variable is on

    vgrp = cell(size(vname));
    [tf,loc] = ismember(dgrp(id), vargroup(:,2));
    vgrp(tf) = vargroup(loc(tf),1);

    % All constant variables go into a single file; they can just be sliced
    % out directly via ncks

    if vflag && cflag
        fprintf('Writing constants.nc...\n');
    end

    if cflag
        newfile = fullfile(outfolder, sprintf('%s_%s_constants.nc', outbase, ftype));    
        if ~exist(newfile, 'file')
            slicevars(vars_constant, fname{cnt}, newfile, false);
            updateattributes(Att, newfile);
        end
    end
    
    % Transfer data for all time-dependent variables
    
    if vflag
        fprintf('Writing variable data...\n');
    end

    skip = ismember(vname, [vars_constant; 'ocean_time']) | ...
           ~ismember(vname, variables);

    for iv = 1:length(vname)
        if ~skip(iv)

            if vflag
                fprintf('%15s: ', vname{iv});
            end
            
            % Start by slicing out the indicated variable from an output
            % file; this is the easiest way to copy all the dimensions,
            % attributes, metadata, etc. The main variable will also be
            % accompanied by the necessary coordinate variables.   

            newfile = fullfile(outfolder, sprintf('%s_%s_%s.nc', outbase, ftype, vname{iv}));

            slicevars(vname(iv), fname{cnt}, newfile);
 
            % 3D variables should be accompanied by zeta.  However, in some
            % of our sims, zeta-archiving was accidentally turned off.  For
            % those, we need to manually add a zeta variable and set it to
            % all-0s.  For now I'm only doing this for average/history
            % files, since the simulations in question didn't save station
            % files.
            
            if ismember(ftype, {'average','history','diagnos'}) && contains(vgrp{iv},'3d')
                In = ncinfo(newfile);
                haszeta = ismember('zeta', {In.Variables.Name});
                if ~haszeta
                    Ig = ncinfo(gridfile);
                    [~,rholoc] = ismember({'xi_rho','eta_rho'}, {Ig.Dimensions.Name});
                    nxirho = Ig.Dimensions(rholoc(1)).Length;
                    netarho = Ig.Dimensions(rholoc(2)).Length;
                    
                    if ismember('xi_rho', {In.Dimensions.Name})
                        % Dimensions exist, just need to add the variable
                        nccreate(newfile, 'zeta', ...
                            'Dimensions', {'xi_rho','eta_rho','ocean_time'}, ...
                            'Datatype', 'single');
                    else
                        % Need to add xi_rho and eta_rho dimensions
                        nccreate(newfile, 'zeta', ...
                            'Dimensions', {'xi_rho', nxirho, ...
                                           'eta_rho', netarho, ...
                                           'ocean_time', Inf}, ...
                            'Datatype', 'single');
                    end
                    
                    % Add attributes
                    
                    switch ftype
                        case 'average'
                            ncwriteatt(newfile, 'zeta', 'long_name', 'time-averaged free-surface');
                        case {'history', 'diagnos'}
                            ncwriteatt(newfile, 'zeta', 'long_name', 'free-surface');
                    end
                    ncwriteatt(newfile, 'zeta', 'units', 'meter');
                    ncwriteatt(newfile, 'zeta', 'comment', 'zeta not archived in this simulation; 0s used as placeholder');
                end
            else
                haszeta = true; % assumed
            end
            
            % If an average or history file, append the relevant lat/lon
            % coordinates to the file
            
            if ismember(ftype, {'average','history','diagnos'})
                switch vgrp{iv}
                    case {'v2d', 'v3d'}
                        cmd = sprintf('ncks -A -v lat_v,lon_v %s %s', gridfile, newfile);
                        system(cmd);
                    case {'u2d', 'u3d'}
                        cmd = sprintf('ncks -A -v lat_u,lon_u %s %s', gridfile, newfile);
                        system(cmd);
                    case {'r2d', 'r3d', 'w3d'}
                        cmd = sprintf('ncks -A -v lat_rho,lon_rho %s %s', gridfile, newfile);
                        system(cmd);
                end
            end
            
            % The new file will need the full record of attributes from all
            % the output files
            
            updateattributes(Att, newfile);
            
            % Check history attribute to see whether data has already been
            % transferred.
            
            hisatt = ncreadatt(newfile, '/', 'history');
            iscomplete = contains(hisatt, 'Data traferred from raw output files via roms_level1.m');
            
            % Copy data from output files to new file, one file at a time
            
            V = ncinfo(newfile, vname{iv});
            ntv = V.Size(strcmp({V.Dimensions.Name}, 'ocean_time'));
            if ntv > nt
                warning('Existing file has more time steps than unique values in dataset');
            end
            
            if ~iscomplete
           
                if vflag
                    c = ConsoleProgressBar;
                    c.setMaximum(nfile);
                    c.start();
                end

                for ifl = 1:nfile
                    if vflag
                        c.setValue(ifl);
                        c.setText(sprintf('%d/%d', ifl, nfile));
                    end

                    if any(isin{ifl})
                        switch ftype
                            case {'average','history','diagnos'}
                                switch vgrp{iv}
                                    case {'r2d', 'u2d', 'v2d'} % variable + time
                                        s = [1 1 sidx(ifl)];
                                        data = ncread(fname{ifl}, vname{iv});
                                        ncwrite(newfile, vname{iv}, data(:,:,isin{ifl}), s);

                                        tdata = ncread(fname{ifl}, 'ocean_time');
                                        ncwrite(newfile, 'ocean_time', tdata(isin{ifl}), sidx(ifl));

                                    case {'r3d', 'w3d', 'u3d', 'v3d'} % variable + time + zeta
                                        s = [1 1 1 sidx(ifl)];
                                        data = ncread(fname{ifl}, vname{iv});
                                        ncwrite(newfile, vname{iv}, data(:,:,:,isin{ifl}), s);

                                        tdata = ncread(fname{ifl}, 'ocean_time');
                                        ncwrite(newfile, 'ocean_time', tdata(isin{ifl}), sidx(ifl));

                                        if haszeta
                                            zdata = ncread(fname{ifl}, 'zeta');
                                            zdata = zdata(:,:,isin{ifl});
                                        else
                                            zdata = zeros(nxirho, netarho, sum(isin{ifl}));
                                        end
                                        ncwrite(newfile, 'zeta', zdata, [1 1 sidx(ifl)]);

                                    otherwise
                                        error('New grid type?');
                                end
                            case {'station'}
                                switch vgrp{iv}
                                    case '2d' % variable + time
                                        s = [1 sidx(ifl)];
                                        data = ncread(fname{ifl}, vname{iv});
                                        ncwrite(newfile, vname{iv}, data(:,isin{ifl}), s);

                                        tdata = ncread(fname{ifl}, 'ocean_time');
                                        ncwrite(newfile, 'ocean_time', tdata(isin{ifl}), sidx(ifl));

                                    case {'r', 'w'} % variable + time + zeta
                                        s = [1 1 sidx(ifl)];
                                        data = ncread(fname{ifl}, vname{iv});
                                        ncwrite(newfile, vname{iv}, data(:,:,isin{ifl}), s);

                                        tdata = ncread(fname{ifl}, 'ocean_time');
                                        ncwrite(newfile, 'ocean_time', tdata(isin{ifl}), sidx(ifl));

                                        zdata = ncread(fname{ifl}, 'zeta');
                                        ncwrite(newfile, 'zeta', zdata(:,isin{ifl}), [1 sidx(ifl)]);

                                end
                        end
                    end
                end
                if vflag
                    c.stop();
                end
                ncaddhis(newfile, 'Data traferred from raw output files via roms_level1.m');
            end
            if vflag
                fprintf('\n');
            end
        end
    end
end

end

%----------------------
% System command 
% wrappers
%----------------------

% Slice out certain variables from a file, via ncks.  Slice only a single
% time slice.  Only runs if potential new file doesn't already exist.

function slicevars(vars, oldfile, newfile, addcoord)

    vstr = sprintf('%s,', vars{:});
    vstr = vstr(1:end-1);
    if nargin<4 || addcoord
        cmd = sprintf('ncks -F -d ocean_time,1,1 -v %s %s %s', vstr, oldfile, newfile);
    else
        cmd = sprintf('ncks -F -C -d ocean_time,1,1 -v %s %s %s', vstr, oldfile, newfile);
    end
    if ~exist(newfile, 'file')
        system(cmd);
    end
end

% Deprecated

function saferun(newfile, cmd)
    
    if exist(newfile, 'file')
        fprintf('(file exists)...');
%         warning('postprocessroms:fileExists', '%s exists, skipping', newfile);
    else
%         fprintf('Will run: %s\n', cmd);
        system(cmd);
    end
end
       
%----------------------
% Depth calc
%----------------------

% Generic s-coordinate transforms

function [zr,zw] = romsdepth(Vt, Cr, Cw, sr, sw, hc, h, zeta)

    switch Vt
        case 1
            z0r = hc .* (sr - Cr) + h .* Cr;
            z0w = hc .* (sw - Cw) + h .* Cw;

            zr = z0r + zeta .* (1 + z0r./h);
            zw = z0w + zeta .* (1 + z0w./h);

        case 2
            z0r = (hc .* sr + Cr .* h)./(h + hc);
            z0w = (hc .* sw + Cw .* h)./(h + hc);

            zr = zetar + (zeta + h) .* z0r;
            zw = zetaw + (zeta + h) .* z0w;
    end

end

%----------------------
% New output file stuff
%----------------------

% % Variable structure
% 
% function V = varstruct(vname, dimnames, atts, type, len)
% 
%     V.Name = vname;
%     if isempty(dimnames)
%         V.Dimensions = [];
%         V.Size = len;
%     else
%         V.Dimensions = struct('Name', dimnames, 'Length', [], 'Unlimited', []);
%         V.Size = [];
%     end
%     V.Attributes = attribstruct(atts{:});
%     V.Datatype = type;
% 
% end
% 
% % Attributes structure
% 
% function A = attribstruct(varargin)
% 
%     atts = reshape(varargin, 2, []);
%     A = struct('Name', atts(1,:), 'Value', atts(2,:));
% 
% end

%-----------------------
% Attribute consolidator
%-----------------------

% Gather global attributes from collection of files

function A = simattributes(oldfl)

    I = cellfun(@ncinfo, oldfl);
    nfl = length(oldfl);

    aname = arrayfun(@(X) {X.Attributes.Name}, I, 'uni', 0);
    aname = unique(cat(2, aname{:}));

    natt = length(aname);

    aval = cell(natt, nfl);
    for ii = 1:nfl

        [~,loc] = ismember({I(ii).Attributes.Name}, aname);
        aval(loc,ii) = {I(ii).Attributes.Value};

    end
    isemp = cellfun(@isempty, aval);

    for ia = 1:natt
        avaltmp = aval(ia,~isemp(ia,:));
        if iscellstr(avaltmp)
            [tmp, itmp] = unique(avaltmp);
            if length(tmp) == 1
                A.(aname{ia}) = tmp{1};
            else
                [~,isrt] = sort(itmp);
                tmp = tmp(isrt);
                A.(aname{ia}) = sprintf('%s\n', tmp{:});
            end
        else
            avaltmp = cat(1, avaltmp{:}); % assuming numeric and can be
            [tmp, itmp] = unique(avaltmp);
            if length(tmp) == 1
                A.(aname{ia}) = tmp(1);
            else
                [~,isrt] = sort(itmp);
                tmp = tmp(isrt);
                A.(aname{ia}) = sprintf('%d, \n', tmp); % if different, convert to character array
            end
        end
            
    end

end

% Update file global attributes to match the old collection

function updateattributes(A, newfl)

    Inew = ncinfo(newfl);

    for ii = 1:length(Inew.Attributes)
        if isfield(A, Inew.Attributes(ii).Name) && ~isequal(Inew.Attributes(ii).Value, A.(Inew.Attributes(ii).Name))
            if strcmp(Inew.Attributes(ii).Name, 'history')
                his1 = Inew.Attributes(ii).Value;
                his2 = A.history;

                his1 = regexp(his1, '\n', 'split');
                his2 = regexp(his2, '\n', 'split');

                isextra = ~ismember(his1, his2);
                his = [his1(isextra) his2(end:-1:1)];
                his = regexprep(sprintf('%s\n', his{:}), '\n+$', '');

                ncwriteatt(newfl, '/', Inew.Attributes(ii).Name, his);
                
%                 fprintf('%s = %s\n\n', Inew.Attributes(ii).Name, his);

            else
                ncwriteatt(newfl, '/', Inew.Attributes(ii).Name, A.(Inew.Attributes(ii).Name));
%                 fprintf('%s = %s\n\n', Inew.Attributes(ii).Name, A.(Inew.Attributes(ii).Name));
            end
        end
    end
end
    



