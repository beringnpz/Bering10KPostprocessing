function roms_level2_rotateuv(varargin)
%ROMS_LEVEL2_ROTATEUV Rotate u/v to east/north
%
% roms_level2_rotateuv(gridfile, lev1folder, outfolder, outbase)
% roms_level2_rotateuv(gridfile, lev1folder, outfolder, outbase, param, value, ...)
%
% This function rotates the u/v from ROMS xi/eta-oriented to east/north
% oriented components, with the new components labeled as uEast and vNorth.
% Also, for average and history files, the values are moved from the u/v
% grids to the rho grid so that the values are co-located rather than
% staggered.
%
% Note: right now only applies to averages and history files, not stations.
%
% Input variables:
%
%   gridfile:   path to ROMS grid file used for the simulation
%
%   lev1folder: path to level 1 postprocessed ROMS output folder (see
%               roms_level1.m).
%
%   outfolder:  path to folder where level 2 postprocessed output from this
%               function will be saved
%
%   outbase:    base filename for newly-created files.  Files will be named
%               <outbase>_<ftype>_<var>.nc, where <var> will be 'uEast' and
%               'vNorth'. The same outbase/ftype formula is used to
%               determine where to find u and v data for the calculations.
%
% Optional input variables (passed as parameter/value pairs):
%
%   verbose:    scalar logical, true to print progress to screen [true]
%
%   ntper:      scalar, integer, number of time steps to process at one
%               time [10]
%
%   ftype:      string or cell array of strings, holding 'average' or
%               'history' to indicate the type of Level 1 output files for
%               which rotated coordinates should be calculated.
%
%   vtype:      string or cell array of strings, indicating rotation of the
%               following variables:
%               'uv':      u/v water column velocity
%               'suvstr':  sustr/svstr surface wind stress
%               'buvstr':  bustr/bvstr bottom stress
%               'uvbar':   ubar/vbar vertically-averaged velocity
%               'uvice':   uice/vice ice velocity


% Copyright 2020 Kelly Kearney

p = inputParser;
p.addRequired('gridfile',    @(x) validateattributes(x, {'char'}, {'scalartext'}));
p.addRequired('lev1folder', @(x) validateattributes(x, {'char'}, {'scalartext'}));
p.addRequired('outfolder',   @(x) validateattributes(x, {'char'}, {'scalartext'}));
p.addRequired('outbase',   @(x) validateattributes(x, {'char'}, {'scalartext'}));
p.addParameter('verbose', true, @(x) validateattributes(x, {'logical'}, {'scalar'}));
p.addParameter('ntper', 10, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer'}));
p.addParameter('ftype', {'average', 'history'}, @(x) ischar(x) || isstring || iscellstr(x));
p.addParameter('vtype', {'velocity'}, @(x) ischar(x) || isstring || iscellstr(x));

p.parse(varargin{:});
Opt = p.Results;

% Verify ftype

if ischar(Opt.ftype)
    Opt.ftype = {Opt.ftype};
end
tf = ismember(Opt.ftype, {'average', 'history'});
if ~all(tf)
    str = sprintf('%s,', Opt.ftype{~tf});
    error('Unrecognized file type: %s', str(1:end-1));
end

% Verify vtype

uvvars = {...
    'uv'        'u'     'v'
    'uvbar'     'ubar'  'vbar'
    'suvstr'    'sustr' 'svstr'
    'buvstr'    'bustr' 'bvstr'
    'uvice'     'uice'  'vice'};

if ischar(Opt.vtype)
    Opt.vtype = {Opt.vtype};
end
[tf,loc] = ismember(Opt.vtype, uvvars(:,1));
if ~all(tf)
    str = sprintf('%s,', Opt.vtype{~tf});
    error('Unrecognized file type: %s', str(1:end-1));
end
uvar = uvvars(loc,2);
vvar = uvvars(loc,3);

% Angles from grid file

angle = ncread(Opt.gridfile, 'angle');

for ii = 1:length(Opt.ftype)
    
    for iv = 1:length(uvar)
        
        is3d = strcmp(uvar{iv}, 'u');
        
        ufile = fullfile(Opt.lev1folder, sprintf('%s_%s_%s.nc', Opt.outbase, Opt.ftype{ii}, uvar{iv}));
        vfile = fullfile(Opt.lev1folder, sprintf('%s_%s_%s.nc', Opt.outbase, Opt.ftype{ii}, vvar{iv}));

        if ~exist(ufile, 'file') && exist(vfile, 'file')
            if Opt.verbose
                fprintf('  %s/%s file(s) missing; skipping rotation', uvar{iv}, vvar{iv});
            end
            continue;
        end
        
        unew = fullfile(Opt.outfolder, sprintf('%s_%s_%sEast.nc',  Opt.outbase, Opt.ftype{ii}, uvar{iv}));
        vnew = fullfile(Opt.outfolder, sprintf('%s_%s_%sNorth.nc', Opt.outbase, Opt.ftype{ii}, vvar{iv}));

        % New files: use ocean_time and s_rho from u/v-file, lat/lon from grid
        % Slice only one time step from ocean_time for now

        if is3d
            vstr = 'ocean_time,s_rho';
        else
            vstr = 'ocean_time';
        end
        
        if exist(ufile,'file') && ~exist(unew,'file')
            cmd = sprintf('ncks -F -v %s -d ocean_time,1,1 %s %s', vstr, ufile, unew);
            system(cmd);

            cmd = sprintf('ncks -A -v lat_rho,lon_rho %s %s', Opt.gridfile, unew);
            system(cmd);
        end
        if exist(vfile,'file') && ~exist(vnew,'file')
            cmd = sprintf('ncks -F -v %s-d ocean_time,1,1 %s %s', vstr, vfile, vnew);
            system(cmd);

            cmd = sprintf('ncks -A -v lat_rho,lon_rho %s %s', Opt.gridfile, vnew);
            system(cmd); 
        end
    
        % Loop over data to shift and rotate

        I = ncinfo(ufile);
        nt = I.Dimensions(strcmp({I.Dimensions.Name}, 'ocean_time')).Length;
        nxi = I.Dimensions(strcmp({I.Dimensions.Name}, 'xi_rho')).Length;
        neta = I.Dimensions(strcmp({I.Dimensions.Name}, 'eta_rho')).Length;
    
        % Create new variables

        uunit  = ncreadatt(ufile, uvar{iv}, 'units');
        ulname = ncreadatt(ufile, uvar{iv}, 'long_name');
        vunit  = ncreadatt(vfile, vvar{iv}, 'units');
        vlname = ncreadatt(vfile, vvar{iv}, 'long_name');
    
        if is3d
            dims = {'xi_rho', 'eta_rho', 's_rho', 'ocean_time'};
        else
            dims = {'xi_rho', 'eta_rho', 'ocean_time'};
        end
        
        uvarnew = [uvar{iv} 'East'];
        vvarnew = [vvar{iv} 'North'];
        
        Iu = ncinfo(unew);
        if ~ismember(uvarnew, {Iu.Variables.Name})
            nccreate(unew, uvarnew, ...
                'Dimensions', dims, ...
                'Datatype', 'single');
            ncwriteatt(unew, uvarnew, 'long_name', sprintf('%s, geo-rotated', ulname));
            ncwriteatt(unew, uvarnew, 'units', uunit);
        end
    
        Iv = ncinfo(vnew);
        if ~ismember(vvarnew, {Iv.Variables.Name})
            nccreate(vnew, vvarnew, ...
                'Dimensions', dims, ...
                'Datatype', 'single');
            ncwriteatt(vnew, vvarnew, 'long_name', sprintf('%s, geo-rotated', vlname));
            ncwriteatt(vnew, vvarnew, 'units', vunit);
        end
    
    
        % Loop over data to shift and rotate

        st = 1:Opt.ntper:nt; % start index
        endi = min(st+Opt.ntper, nt);
        cn = endi - st + 1;
        nchunk = length(st);

        if Opt.verbose
            fprintf('Rotating %s/%s...\n', uvar{iv}, vvar{iv});
            c = ConsoleProgressBar;
            c.setMaximum(nchunk);
            c.start();
        end
    
        for ic = 1:nchunk

            if Opt.verbose
                c.setValue(ic);
                c.setText(sprintf('%d/%d', ic, nchunk));
            end

            if is3d
                scin  = {[1 1 1 st(ic)], [Inf Inf Inf cn(ic)]}; % start, count
                scout = [1 1 1 st(ic)]; % start
            else
                scin  = {[1 1 st(ic)], [Inf Inf Inf cn(ic)]}; % start, count
                scout = [1 1 st(ic)]; % start
            end

            u = ncread(ufile, uvar{iv}, scin{:});
            v = ncread(vfile, vvar{iv}, scin{:});

            % Move to rho grid

            if is3d

                urho = nan(nxi, neta, size(u,3), size(u,4));
                vrho = nan(nxi, neta, size(v,3), size(v,4));

                urho(2:end-1,:,:,:) = 0.5 .* (u(1:end-1,:,:,:) + u(2:end,:,:,:));
                vrho(:,2:end-1,:,:) = 0.5 .* (v(:,1:end-1,:,:) + v(:,2:end,:,:));
            else
                urho = nan(nxi, neta, size(u,3));
                vrho = nan(nxi, neta, size(v,3));

                urho(2:end-1,:,:) = 0.5 .* (u(1:end-1,:,:) + u(2:end,:,:));
                vrho(:,2:end-1,:) = 0.5 .* (v(:,1:end-1,:) + v(:,2:end,:));
            end

            % Rotate

            uvitheta = (urho + 1i.*vrho).*exp(1i.*angle);
            uEast = real(uvitheta);
            vNorth = imag(uvitheta);

            % Write

            ncwrite(unew, uvarnew,  uEast,  scout);
            ncwrite(vnew, vvarnew, vNorth, scout);

        end
    
        if Opt.verbose
            c.stop();
            fprintf('\n');
        end
    
        tdata = ncread(ufile, 'ocean_time');
        ncwrite(unew, 'ocean_time', tdata);
        ncwrite(vnew, 'ocean_time', tdata);
    
%     uveitheta = (u_roms+1i*v_roms).*exp(1i*angle); % 1i = sqrt(-1) in Matlab
%     u_east = real(uveitheta);
%     v_north = imag(uveitheta);
    
%     Urho(i,j) = 0.5 * ( U(i,j) + U(i+1,j) )
%     Vrho(i,j) = 0.5 * ( V(i,j) + V(i,j+1) )
    end
end



