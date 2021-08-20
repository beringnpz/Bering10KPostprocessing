function roms_level2_depths(varargin)
%ROMS_LEVEL2_DEPTHS Calculate s-coordinate depth values
%
% roms_level2_depths(gridfile, lev1folder, outfolder, outbase);
% roms_level2_depths(gridfile, lev1folder, outfolder, outbase, param, value, ...);
%
% This function adds a file holding depth coordinate values.  While these
% values can be calculated from the original coordinate data coupled with
% the ROMS and CF-standard s-coordinate formulas, we pre-calculate the
% values here for ease of use.   
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
%               <outbase>_<ftype>_<var>.nc, where <var> indicates one 8
%               depth variables: z_rho, z_w, z_psi, z_u, z_v, z_psi_w,
%               z_u_w, z_v_w.  The same outbase/ftype formula is used to
%               determine where to find zeta and ocean_time data for the
%               depth calculations.
%
% Optional input variables (passed as parameter/value pairs):
%
%   verbose:    scalar logical, true to print progress to screen [true]
%
%   ntper:      scalar, integer, number of time steps to process at one
%               time [10]
%
%   ftype:      string or cell array of strings, holding 'average',
%               'history', and/or 'station' to indicate the type of Level 1
%               output files for which depths should be calculated
%
%   zvars:      string or cell array of strings, indicating which depth
%               variables to calculate.  If not included, all potential
%               depth variables for the specified file type will be
%               calculated.  

% Copyright 2020 Kelly Kearney


%----------------------
% Parse and check input
%----------------------

% Input parsing

p = inputParser;
p.addRequired('gridfile',    @(x) validateattributes(x, {'char'}, {'scalartext'}));
p.addRequired('lev1folder', @(x) validateattributes(x, {'char'}, {'scalartext'}));
p.addRequired('outfolder',   @(x) validateattributes(x, {'char'}, {'scalartext'}));
p.addRequired('outbase',   @(x) validateattributes(x, {'char'}, {'scalartext'}));
p.addParameter('verbose', true, @(x) validateattributes(x, {'logical'}, {'scalar'}));
p.addParameter('ntper', 10, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer'}));
p.addParameter('ftype', {'average', 'history', 'station'}, @(x) ischar(x) || iscellstr(x));
p.addParameter('zvars', {}, @(x) ischar(x) || iscellstr(x));

p.parse(varargin{:});
Opt = p.Results;

% Output folder must exist previously (safety precaution...)

if ~exist(Opt.outfolder, 'dir')
    error('Output folder (%s) not found; must exist prior to running this script', Opt.outfolder);
end

% Verify ftype

if ischar(Opt.ftype)
    Opt.ftype = {Opt.ftype};
end
tf = ismember(Opt.ftype, {'average', 'history', 'station'});
if ~all(tf)
    str = sprintf('%s,', Opt.ftype{~tf});
    error('Unrecognized file type: %s', str(1:end-1));
end

%----------------------
% File creation
%----------------------

% Info about new depth variables: type, variable, variables from grid file,
% variable from constants file

zvars = {...
    'average'   'z_rho'     {'lat_rho','lon_rho'}   {'s_rho'}
    'average'   'z_w'       {'lat_rho','lon_rho'}   {'s_w'}
    'average'   'z_psi'     {'lat_psi','lon_psi'}   {'s_rho'}   
    'average'   'z_u'       {'lat_u'  ,'lon_u'  }   {'s_rho'}
    'average'   'z_v'       {'lat_v'  ,'lon_v'  }   {'s_rho'}
    'average'   'z_psi_w'   {'lat_psi','lon_psi'}   {'s_w'}
    'average'   'z_u_w'     {'lat_u'  ,'lon_u'  }   {'s_w'}
    'average'   'z_v_w'     {'lat_v'  ,'lon_v'  }   {'s_w'}
    'history'   'z_rho'     {'lat_rho','lon_rho'}   {'s_rho'}
    'history'   'z_w'       {'lat_rho','lon_rho'}   {'s_w'}
    'history'   'z_psi'     {'lat_psi','lon_psi'}   {'s_rho'}
    'history'   'z_u'       {'lat_u'  ,'lon_u'  }   {'s_rho'}
    'history'   'z_v'       {'lat_v'  ,'lon_v'  }   {'s_rho'}
    'history'   'z_psi_w'   {'lat_psi','lon_psi'}   {'s_w'}
    'history'   'z_u_w'     {'lat_u'  ,'lon_u'  }   {'s_w'}
    'history'   'z_v_w'     {'lat_v'  ,'lon_v'  }   {'s_w'}
    'station'   'z_rho'     {}                      {'lat_rho','lon_rho','s_rho'}
    'station'   'z_w'       {}                      {'lat_rho','lon_rho','s_w'}
};

zisin = ismember(zvars(:,1), Opt.ftype);

if ~isempty(Opt.zvars)
    if ischar(Opt.zvars)
        Opt.zvars = {Opt.zvars};
    end
    tf = ismember(Opt.zvars, zvars(zisin,2));
    if ~all(tf)
        str = sprintf('%s,', Opt.zvars{~tf});
        error('Input zvar value(s) do not match file type: %s', str(1:end-1));
    end
    zisin = zisin & ismember(zvars(:,2), Opt.zvars);
end

zvars = zvars(zisin,:);

% Figure out reference files (constants, zeta) for each type

for itp = 1:length(Opt.ftype)

    % Constants file
    
    Ref.(Opt.ftype{itp}).cfile = fullfile(Opt.lev1folder, sprintf('%s_%s_constants.nc', Opt.outbase, Opt.ftype{itp}));
  
    % Find zeta file
    
    Ref.(Opt.ftype{itp}).zetafile = '';
    
    ztest = fullfile(Opt.lev1folder, sprintf('%s_%s_zeta.nc', Opt.outbase, Opt.ftype{itp}));
    if exist(ztest, 'file')
        Ref.(Opt.ftype{itp}).zetafile = ztest;
    else
        pth = fullfile(Opt.lev1folder, sprintf('%s_%s_*.nc', Opt.outbase, Opt.ftype{itp})); 
        F = dir(pth);
        for ii = 1:length(F)
            Itmp = ncinfo(fullfile(F(ii).folder,F(ii).name));
            haszeta = ismember('zeta', {Itmp.Variables.Name});
            if haszeta
                Ref.(Opt.ftype{itp}).zetafile = fullfile(F(ii).folder,F(ii).name);
                break
            end
        end
    end
    
    Ref.(Opt.ftype{itp}).flag = exist(Ref.(Opt.ftype{itp}).cfile, 'file') && ...
                            ~isempty(Ref.(Opt.ftype{itp}).zetafile);
    
end

% Create new input files by gathering dimensions and dimension variables
% from the constants file and grid file, as necessary.  Then add a new
% variable for the to-be-calculated depths.

for ii = 1:size(zvars,1)
    
    newfile = fullfile(Opt.outfolder, sprintf('%s_%s_%s.nc', Opt.outbase, zvars{ii,1}, zvars{ii,2}));
    
    if ~exist(newfile,'file') && Ref.(zvars{ii,1}).flag
        
        % Slice constants variables (also brings along global atts)
        
        cvarstr = sprintf('%s,', zvars{ii,4}{:});
        cvarstr = cvarstr(1:end-1);
        cmd = sprintf('ncks -v %s %s %s', cvarstr, Ref.(zvars{ii,1}).cfile, newfile);
        system(cmd);
        
        % Append grid variables 
        
        if ~isempty(zvars{ii,3})
            gvarstr = sprintf('%s,', zvars{ii,3}{:});
            gvarstr = gvarstr(1:end-1);
            cmd = sprintf('ncks -A -v %s %s %s', gvarstr, Opt.gridfile, newfile);
            system(cmd);
        end
        
        % Append ocean_time (for now, just one time step... adding it all
        % at once causes the new variables to be created in full via
        % nccreate, which is much slower than growing it as we add data).
        
        cmd = sprintf('ncks -A -F -d ocean_time,1,1 -v ocean_time %s %s', Ref.(zvars{ii,1}).zetafile, newfile);
        system(cmd);
        
        % At this point, all necessary dimensions are in the file, and we
        % can create the new variable and calculate values
        
        switch sprintf('%s_%s', zvars{ii,1:2})
            case {'average_z_rho', 'history_z_rho'}
                nccreate(newfile, zvars{ii,2}, ...
                    'Dimensions', {'xi_rho', ...
                                   'eta_rho', ...
                                   's_rho', ...
                                   'ocean_time'}, ...
                    'Datatype', 'single');
                ncwriteatt(newfile, zvars{ii,2}, 'long_name', 'depth at horizontal rho points, rho depth');
                ncwriteatt(newfile, zvars{ii,2}, 'units', 'm');
                ncwriteatt(newfile, zvars{ii,2}, 'coordinates', 'lon_rho lat_rho s_rho ocean_time');
                                    
            case {'average_z_w', 'history_z_w'}
                nccreate(newfile, zvars{ii,2}, ...
                    'Dimensions', {'xi_rho', ...
                                   'eta_rho', ...
                                   's_w', ...
                                   'ocean_time'}, ...
                    'Datatype', 'single'); 
                ncwriteatt(newfile, zvars{ii,2}, 'long_name', 'depth at horizontal rho points, w depth');
                ncwriteatt(newfile, zvars{ii,2}, 'units', 'm');
                ncwriteatt(newfile, zvars{ii,2}, 'coordinates', 'lon_rho lat_rho s_w ocean_time');
                               
            case {'average_z_psi', 'history_z_psi'}
                nccreate(newfile, zvars{ii,2}, ...
                    'Dimensions', {'xi_psi', ...
                                   'eta_psi', ...
                                   's_rho', ...
                                   'ocean_time'}, ...
                    'Datatype', 'single');
                ncwriteatt(newfile, zvars{ii,2}, 'long_name', 'depth at horizontal psi points, rho depth');
                ncwriteatt(newfile, zvars{ii,2}, 'units', 'm');
                ncwriteatt(newfile, zvars{ii,2}, 'coordinates', 'lon_psi lat_psi s_rho ocean_time');
                               
            case {'average_z_u', 'history_z_u'}
                nccreate(newfile, zvars{ii,2}, ...
                    'Dimensions', {'xi_u', ...
                                   'eta_u', ...
                                   's_rho', ...
                                   'ocean_time'}, ...
                    'Datatype', 'single');
                ncwriteatt(newfile, zvars{ii,2}, 'long_name', 'depth at horizontal u points, rho depth');
                ncwriteatt(newfile, zvars{ii,2}, 'units', 'm');
                ncwriteatt(newfile, zvars{ii,2}, 'coordinates', 'lon_u lat_u s_rho ocean_time');
                               
            case {'average_z_v', 'history_z_v'}
                nccreate(newfile, zvars{ii,2}, ...
                    'Dimensions', {'xi_v', ...
                                   'eta_v', ...
                                   's_rho', ...
                                   'ocean_time'}, ...
                    'Datatype', 'single');
                ncwriteatt(newfile, zvars{ii,2}, 'long_name', 'depth at horizontal v points, rho depth');
                ncwriteatt(newfile, zvars{ii,2}, 'units', 'm');
                ncwriteatt(newfile, zvars{ii,2}, 'coordinates', 'lon_v lat_v s_rho ocean_time');
                               
            case {'average_z_psi_w', 'history_z_psi_w'}
                nccreate(newfile, zvars{ii,2}, ...
                    'Dimensions', {'xi_psi', ...
                                   'eta_psi', ...
                                   's_w', ...
                                   'ocean_time'}, ...
                    'Datatype', 'single');
                ncwriteatt(newfile, zvars{ii,2}, 'long_name', 'depth at horizontal psi points, w depth');
                ncwriteatt(newfile, zvars{ii,2}, 'units', 'm');
                ncwriteatt(newfile, zvars{ii,2}, 'coordinates', 'lon_psi lat_psi s_w ocean_time');
                
            case {'average_z_u_w', 'history_z_u_w'}
                nccreate(newfile, zvars{ii,2}, ...
                    'Dimensions', {'xi_u', ...
                                   'eta_u', ...
                                   's_w', ...
                                   'ocean_time'}, ...
                    'Datatype', 'single');
                ncwriteatt(newfile, zvars{ii,2}, 'long_name', 'depth at horizontal u points, w depth');
                ncwriteatt(newfile, zvars{ii,2}, 'units', 'm');
                ncwriteatt(newfile, zvars{ii,2}, 'coordinates', 'lon_u lat_u s_w ocean_time');
                
            case {'average_z_v_w', 'history_z_v_w'}
                nccreate(newfile, zvars{ii,2}, ...
                    'Dimensions', {'xi_v', ...
                                   'eta_v', ...
                                   's_w', ...
                                   'ocean_time'}, ...
                    'Datatype', 'single');
                ncwriteatt(newfile, zvars{ii,2}, 'long_name', 'depth at horizontal v points, w depth');
                ncwriteatt(newfile, zvars{ii,2}, 'units', 'm');
                ncwriteatt(newfile, zvars{ii,2}, 'coordinates', 'lon_v lat_v s_w ocean_time');
                
            case {'station_z_rho'}
                nccreate(newfile, zvars{ii,2}, ...
                    'Dimensions', {'station', ...
                                   's_rho', ...
                                   'ocean_time'}, ...
                    'Datatype', 'single');
                ncwriteatt(newfile, zvars{ii,2}, 'long_name', 'depth at stations, rho depth');
                ncwriteatt(newfile, zvars{ii,2}, 'units', 'm');
                
            case {'station_z_w'}
                nccreate(newfile, zvars{ii,2}, ...
                    'Dimensions', {'station', ...
                                   's_w', ...
                                   'ocean_time'}, ...
                    'Datatype', 'single');
                ncwriteatt(newfile, zvars{ii,2}, 'long_name', 'depth at stations, w depth');
                ncwriteatt(newfile, zvars{ii,2}, 'units', 'm');
            
        end
        
    end
    
end

%----------------------
% Calculate depths
%----------------------

% Do calculation in chunks to avoid memory problems

for ii = 1:length(Opt.ftype)

    if Ref.(Opt.ftype{ii}).flag
    
        if Opt.verbose
            fprintf('Calculating depth values...%s\n', Opt.ftype{ii});
        end

        % Calculate start/count for each write chunk

        I = ncinfo(Ref.(Opt.ftype{ii}).zetafile, 'ocean_time');
        nt = I.Size;

        st = 1:Opt.ntper:nt; % start index
        endi = min(st+Opt.ntper, nt);
        cn = endi - st + 1;
        nchunk = length(st);

        % Set status printing

        if Opt.verbose
            c = ConsoleProgressBar;
            c.setMaximum(nchunk);
            c.start();
        end

        % File name shortcut

        fl = @(v) fullfile(Opt.outfolder, sprintf('%s_%s_%s.nc', Opt.outbase, Opt.ftype{ii}, v));

        % Primary calculations

        switch Opt.ftype{ii}
            case {'average', 'history'}

                Tmp = ncstruct(Ref.(Opt.ftype{ii}).cfile, ...
                    'Vtransform', 'hc', 's_rho', 's_w', 'Cs_r', 'Cs_w');

                Tmp.h = ncread(Opt.gridfile, 'h');

                [Lp, Mp] = size(Tmp.h);
                L=Lp-1;
                M=Mp-1;

                % Need depths at rho, psi, u, and v points
                % All input on xi x eta x depth x time grid

                Cr = permute(Tmp.Cs_r, [2 3 1]);
                Cw = permute(Tmp.Cs_w, [2 3 1]);
                sr = permute(Tmp.s_rho, [2 3 1]);
                sw = permute(Tmp.s_w, [2 3 1]);

                zwrite = ismember({'z_rho','z_w','z_psi','z_psi_w', 'z_u','z_u_w','zv','z_v_w'}, zvars(:,2));

                for ic = 1:nchunk

                    if Opt.verbose
                        c.setValue(ic);
                        c.setText(sprintf('%d/%d', ic, nchunk));
                    end

                    scin  = {[1 1 st(ic)], [Inf Inf cn(ic)]}; % start, count
                    scout = [1 1 1 st(ic)]; % start

                    zeta = ncread(Ref.(Opt.ftype{ii}).zetafile, 'zeta', scin{:});

                    hr = Tmp.h;
                    zetar = zeta;

                    hp    = 0.25.*(hr(  1:L,1:M  )+hr(  2:Lp,1:M  )+hr(  1:L,2:Mp  )+hr(  2:Lp,2:Mp  ));
                    zetap = 0.25.*(zeta(1:L,1:M,:)+zeta(2:Lp,1:M,:)+zeta(1:L,2:Mp,:)+zeta(2:Lp,2:Mp,:));

                    hu    = 0.5.*(hr(  1:L,1:Mp  )+hr(  2:Lp,1:Mp  ));
                    zetau = 0.5.*(zeta(1:L,1:Mp,:)+zeta(2:Lp,1:Mp,:));

                    hv    = 0.5.*(hr(  1:Lp,1:M  )+hr(  1:Lp,2:Mp  ));
                    zetav = 0.5.*(zeta(1:Lp,1:M,:)+zeta(1:Lp,2:Mp,:));

                    [zrr,zrw] = romsdepth(Tmp.Vtransform, Cr, Cw, sr, sw, Tmp.hc, ...
                                          hr, permute(zetar, [1 2 4 3]));
                    [zpr,zpw] = romsdepth(Tmp.Vtransform, Cr, Cw, sr, sw, Tmp.hc, ...
                                          hp, permute(zetap, [1 2 4 3]));
                    [zur,zuw] = romsdepth(Tmp.Vtransform, Cr, Cw, sr, sw, Tmp.hc, ...
                                          hu, permute(zetau, [1 2 4 3]));
                    [zvr,zvw] = romsdepth(Tmp.Vtransform, Cr, Cw, sr, sw, Tmp.hc, ...
                                          hv, permute(zetav, [1 2 4 3]));
                    
                    
                    if zwrite(1); ncwrite(fl('z_rho'  ), 'z_rho',   zrr, scout); end
                    if zwrite(2); ncwrite(fl('z_w'    ), 'z_w',     zrw, scout); end
                    if zwrite(3); ncwrite(fl('z_psi'  ), 'z_psi',   zpr, scout); end
                    if zwrite(4); ncwrite(fl('z_psi_w'), 'z_psi_w', zpw, scout); end
                    if zwrite(5); ncwrite(fl('z_u'    ), 'z_u',     zur, scout); end
                    if zwrite(6); ncwrite(fl('z_u_w'  ), 'z_u_w',   zuw, scout); end
                    if zwrite(7); ncwrite(fl('z_v'    ), 'z_v',     zvr, scout); end
                    if zwrite(8); ncwrite(fl('z_v_w'  ), 'z_v_w',   zvw, scout); end
                end
                if Opt.verbose
                    c.stop();
                    fprintf('\n');
                end

                tdata = ncread(Ref.(Opt.ftype{ii}).zetafile, 'ocean_time');
                if zwrite(1); ncwrite(fl('z_rho'  ), 'ocean_time', tdata); end
                if zwrite(2); ncwrite(fl('z_w'    ), 'ocean_time', tdata); end
                if zwrite(3); ncwrite(fl('z_psi'  ), 'ocean_time', tdata); end
                if zwrite(4); ncwrite(fl('z_psi_w'), 'ocean_time', tdata); end
                if zwrite(5); ncwrite(fl('z_u'    ), 'ocean_time', tdata); end
                if zwrite(6); ncwrite(fl('z_u_w'  ), 'ocean_time', tdata); end
                if zwrite(7); ncwrite(fl('z_v'    ), 'ocean_time', tdata); end
                if zwrite(8); ncwrite(fl('z_v_w'  ), 'ocean_time', tdata); end

                if zwrite(1); ncaddhis(fl('z_rho'  ), 'z_rho variable added via roms_level2_depths.m'); end
                if zwrite(2); ncaddhis(fl('z_w'    ), 'z_w variable added via roms_level2_depths.m'); end
                if zwrite(3); ncaddhis(fl('z_psi'  ), 'z_psi variable added via roms_level2_depths.m'); end
                if zwrite(4); ncaddhis(fl('z_psi_w'), 'z_psi_w variable added via roms_level2_depths.m'); end
                if zwrite(5); ncaddhis(fl('z_u'    ), 'z_u variable added via roms_level2_depths.m'); end
                if zwrite(6); ncaddhis(fl('z_u_w'  ), 'z_u_w variable added via roms_level2_depths.m'); end
                if zwrite(7); ncaddhis(fl('z_v'    ), 'z_v variable added via roms_level2_depths.m'); end
                if zwrite(8); ncaddhis(fl('z_v_w'  ), 'z_v_w variable added via roms_level2_depths.m'); end
                
            case {'station'}

                Tmp = ncstruct(Ref.(Opt.ftype{ii}).cfile, ...
                    'Vtransform', 'hc', 's_rho', 's_w', 'Cs_r', 'Cs_w', 'h');

                zwrite = ismember({'z_rho','z_w'}, zvars(:,2));

                for ic = 1:nchunk

                    if Opt.verbose
                        c.setValue(ic);
                        c.setText(sprintf('%d/%d', ic, nchunk));
                    end

                    scin  = {[1 st(ic)], [Inf cn(ic)]}; % start, count
                    scout = [1 1 st(ic)]; % start

                    zeta = ncread(Ref.(Opt.ftype{ii}).zetafile, 'zeta', scin{:});

                    [zr,zw] = romsdepth(Tmp.Vtransform, ...
                                    permute(Tmp.Cs_r,  [2 1]), ...
                                    permute(Tmp.Cs_w,  [2 1]), ...
                                    permute(Tmp.s_rho, [2 1]), ...
                                    permute(Tmp.s_w,   [2 1]), ...
                                    Tmp.hc, ...
                                    Tmp.h, ...
                                    permute(zeta, [1 3 2]));

                    if zwrite(1); ncwrite(fl('z_rho'), 'z_rho',   zr, scout); end
                    if zwrite(2); ncwrite(fl('z_w'  ), 'z_w',     zw, scout); end
                end
                if Opt.verbose
                    c.stop();
                    fprintf('\n');
                end

                % Transfer full ocean_time data now

                tdata = ncread(Ref.(Opt.ftype{ii}).zetafile, 'ocean_time');
                if zwrite(1); ncwrite(fl('z_rho'), 'ocean_time', tdata); end
                if zwrite(2); ncwrite(fl('z_w'  ), 'ocean_time', tdata); end
                
                if zwrite(1); ncaddhis(fl('z_rho'  ), 'z_rho variable added via roms_level2_depths.m'); end
                if zwrite(2); ncaddhis(fl('z_w'    ), 'z_w variable added via roms_level2_depths.m'); end

        end
    
    end
    
end



end

%     % Calculate depths
%     
%     if flag(itp)
% 
% %         if ~exist(afol, 'dir')
% %             mkdir(afol);
% %         end
%         if Opt.verbose
%             fprintf('Creating depths file: %s...\n', ftype{itp});
%         end
%         depthvars(fol{itp}, Opt.outfolder, Opt.outbase, Opt.gridfile, Opt.verbose, cfile, zetafile);
%     end
% end


% end

%----------------------
% Create depth file and
% calculate depths
%----------------------


% function depthvars(ftype, outfolder, outbase, gridfile, vflag, cfile,zetafile)
% 
%     % File name (and check if it exists)
% 
%     zfile = fullfile(outfolder, ftype, 'depths.nc');
%     haszfile = exist(zfile, 'file');
% 
%     nt = ncinfo(zetafile, 'ocean_time');
%     nt = nt.Size;
%     
%     if haszfile
%         switch ftype
%             case {'average', 'history'}
%                 V = ncinfo(zfile, 'z_v_w'); % Check last variable written
%             case {'station'}
%                 V = ncinfo(zfile, 'z_w');
%         end
%         ntv = V.Size(strcmp({V.Dimensions.Name}, 'ocean_time'));
%     end
% 
%     % Gather dimensions
%     
%     I(1) = ncinfo(gridfile);
%     I(2) = ncinfo(zetafile);
%     I(3) = ncinfo(cfile);
%     D = [I(1).Dimensions I(2).Dimensions I(3).Dimensions];
%     D = table2struct(unique(struct2table(D)));
%     
%     % First, get depth parameters and create new file, if needed
% 
%     Inew = struct;
%     Inew.Name = '/';
%     Inew.Format = 'classic';
%     Inew.Attributes = I(3).Attributes;
% 
%     switch ftype
%         case {'average', 'history'}
% 
%             Tmp = ncstruct(cfile, ...
%                 'Vtransform', 'hc', 's_rho', 's_w', 'Cs_r', 'Cs_w');
% 
%             Tmp.h = ncread(gridfile, 'h');
% 
%             [Lp, Mp] = size(Tmp.h);
%             L=Lp-1;
%             M=Mp-1;
% 
%             % Need depths at rho, psi, u, and v points
%             % All input on xi x eta x depth x time grid
% 
%             Cr = permute(Tmp.Cs_r, [2 3 1]);
%             Cw = permute(Tmp.Cs_w, [2 3 1]);
%             sr = permute(Tmp.s_rho, [2 3 1]);
%             sw = permute(Tmp.s_w, [2 3 1]);
% 
%             % Create file
% 
%             dname = {'xi_rho', 'eta_rho', ...
%                      'xi_psi', 'eta_psi', ...
%                      'xi_u',   'eta_u', ...
%                      'xi_v',   'eta_v', ...
%                      's_rho',  's_w', 'ocean_time'};
%             [tf,loc] = ismember(dname, {D.Name});
%             if ~all(tf)
%                 error('Not all necessary dimensions found');
%             end
%             Inew.Dimensions = D(loc);
% 
%             Inew.Variables = [...
%                 varstruct('z_rho', ...  
%                           {'xi_rho', 'eta_rho', 's_rho', 'ocean_time'}, ...
%                           {'long_name', 'depth at horizontal rho points, rho depth', ... 
%                            'units', 'm'}, ...
%                           'double', ...
%                           [])
%                 varstruct('z_w', ...  
%                           {'xi_rho', 'eta_rho', 's_w', 'ocean_time'}, ...
%                           {'long_name', 'depth at horizontal rho points, w depth', ... 
%                            'units', 'm'}, ...
%                           'double', ...
%                           [])
%                 varstruct('z_psi', ...  
%                           {'xi_psi', 'eta_psi', 's_rho', 'ocean_time'}, ...
%                           {'long_name', 'depth at horizontal psi points, rho depth', ... 
%                            'units', 'm'}, ...
%                           'double', ...
%                           [])
%                 varstruct('z_u', ...  
%                           {'xi_u', 'eta_u', 's_rho', 'ocean_time'}, ...
%                           {'long_name', 'depth at horizontal u points, rho depth', ... 
%                            'units', 'm'}, ...
%                           'double', ...
%                           [])
%                 varstruct('z_v', ...  
%                           {'xi_v', 'eta_v', 's_rho', 'ocean_time'}, ...
%                           {'long_name', 'depth at horizontal v points, rho depth', ... 
%                            'units', 'm'}, ...
%                           'double', ...
%                           [])
%                 varstruct('z_psi_w', ...  
%                           {'xi_psi', 'eta_psi', 's_w', 'ocean_time'}, ...
%                           {'long_name', 'depth at horizontal psi points, w depth', ... 
%                            'units', 'm'}, ...
%                           'double', ...
%                           [])
%                 varstruct('z_u_w', ...  
%                           {'xi_u', 'eta_u', 's_w', 'ocean_time'}, ...
%                           {'long_name', 'depth at horizontal u points, w depth', ... 
%                            'units', 'm'}, ...
%                           'double', ...
%                           [])
%                 varstruct('z_v_w', ...  
%                           {'xi_v', 'eta_v', 's_w', 'ocean_time'}, ...
%                           {'long_name', 'depth at horizontal v points, w depth', ... 
%                            'units', 'm'}, ...
%                           'double', ...
%                           [])
%                 ];
%             Inew = updatencschema(Inew);
% 
%             if ~haszfile
%                 ncwriteschema(zfile, Inew);
%                 ncaddhis(zfile, 'File created by roms_level2_depths.m');
%             end
% 
%         case {'station'}
% 
%             Tmp = ncstruct(cfile, ...
%                 'Vtransform', 'hc', 's_rho', 's_w', 'Cs_r', 'Cs_w', 'h');
% 
%             % Station files only need depths at station locations
%             % All input on station x depth x time grid
% 
%             % Write to file
% 
%             dname = {'station','s_w','s_rho','ocean_time'};
%             [tf,loc] = ismember(dname, {D.Name});
%             if ~all(tf)
%                 error('Not all necessary dimensions found');
%             end
%             Inew.Dimensions = D(loc);
% 
%             Inew.Variables = [...
%                 varstruct('z_rho', ...  
%                           {'station', 's_rho', 'ocean_time'}, ...
%                           {'long_name', 'depth at stations, rho depth', ... 
%                            'units', 'm'}, ...
%                           'double', ...
%                           [])
%                 varstruct('z_w', ...  
%                           {'station', 's_w', 'ocean_time'}, ...
%                           {'long_name', 'depth stations, w depth', ... 
%                            'units', 'm'}, ...
%                           'double', ...
%                           [])
%                 ];
%             Inew = updatencschema(Inew);
% 
%             zfile = fullfile(outfolder, ftype, 'depths.nc');
%             if ~haszfile
%                 ncwriteschema(zfile, Inew);
%             end
%     end
% 
%     % Next, calculate depths and write to file.  Do this in chunks to
%     % avoid memory problems.
% 
%     if vflag
%         fprintf('Calculating depth values...\n');
%     end
% 
%     dt = 10;
%     st = 1:dt:nt; % start index
%     endi = min(st+dt, nt);
%     cn = endi - st + 1;
%     nchunk = length(st);
% 
%     if vflag
%         c = ConsoleProgressBar;
%         c.setMaximum(nchunk);
%         c.start();
%     end
% 
%     switch ftype
%         case {'average', 'history'}
% 
%             for ic = 1:nchunk
% 
%                 if vflag
%                     c.setValue(ic);
%                     c.setText(sprintf('%d/%d', ic, nchunk));
%                 end
% 
%                 scin  = {[1 1 st(ic)], [Inf Inf cn(ic)]}; % start, count
%                 scout = [1 1 1 st(ic)]; % start
% 
%                 zeta = ncread(zetafile, 'zeta', scin{:});
% 
%                 hr = Tmp.h;
%                 zetar = zeta;
% 
%                 hp    = 0.25.*(hr(  1:L,1:M  )+hr(  2:Lp,1:M  )+hr(  1:L,2:Mp  )+hr(  2:Lp,2:Mp  ));
%                 zetap = 0.25.*(zeta(1:L,1:M,:)+zeta(2:Lp,1:M,:)+zeta(1:L,2:Mp,:)+zeta(2:Lp,2:Mp,:));
% 
%                 hu    = 0.5.*(hr(  1:L,1:Mp  )+hr(  2:Lp,1:Mp  ));
%                 zetau = 0.5.*(zeta(1:L,1:Mp,:)+zeta(2:Lp,1:Mp,:));
% 
%                 hv    = 0.5.*(hr(  1:Lp,1:M  )+hr(  1:Lp,2:Mp  ));
%                 zetav = 0.5.*(zeta(1:Lp,1:M,:)+zeta(1:Lp,2:Mp,:));
% 
%                 [zrr,zrw] = romsdepth(Tmp.Vtransform, Cr, Cw, sr, sw, Tmp.hc, ...
%                                       hr, permute(zetar, [1 2 4 3]));
%                 [zpr,zpw] = romsdepth(Tmp.Vtransform, Cr, Cw, sr, sw, Tmp.hc, ...
%                                       hp, permute(zetap, [1 2 4 3]));
%                 [zur,zuw] = romsdepth(Tmp.Vtransform, Cr, Cw, sr, sw, Tmp.hc, ...
%                                       hu, permute(zetau, [1 2 4 3]));
%                 [zvr,zvw] = romsdepth(Tmp.Vtransform, Cr, Cw, sr, sw, Tmp.hc, ...
%                                       hv, permute(zetav, [1 2 4 3]));
% 
%                 ncwrite(zfile, 'z_rho',   zrr, scout);
%                 ncwrite(zfile, 'z_w',     zrw, scout);
%                 ncwrite(zfile, 'z_psi',   zpr, scout);
%                 ncwrite(zfile, 'z_psi_w', zpw, scout);
%                 ncwrite(zfile, 'z_u',     zur, scout);
%                 ncwrite(zfile, 'z_u_w',   zuw, scout);
%                 ncwrite(zfile, 'z_v',     zvr, scout);
%                 ncwrite(zfile, 'z_v_w',   zvw, scout);
%             end
%             if vflag
%                 c.stop();
%                 fprintf('\n');
%             end
% 
%         case {'station'}
% 
%             for ic = 1:nchunk
% 
%                 if vflag
%                     c.setValue(ic);
%                     c.setText(sprintf('%d/%d', ic, nchunk));
%                 end
% 
%                 scin  = {[1 st(ic)], [Inf cn(ic)]}; % start, count
%                 scout = [1 1 st(ic)]; % start
% 
%                 zeta = ncread(zetafile, 'zeta', scin{:});
% 
%                 [zr,zw] = romsdepth(Tmp.Vtransform, ...
%                                 permute(Tmp.Cs_r,  [2 1]), ...
%                                 permute(Tmp.Cs_w,  [2 1]), ...
%                                 permute(Tmp.s_rho, [2 1]), ...
%                                 permute(Tmp.s_w,   [2 1]), ...
%                                 Tmp.hc, ...
%                                 Tmp.h, ...
%                                 permute(zeta, [1 3 2]));
%                 ncwrite(zfile, 'z_rho',   zr, scout);
%                 ncwrite(zfile, 'z_w',     zw, scout);
%             end
%             if vflag
%                 c.stop();
%                 fprintf('\n');
%             end
%     end
% end


%----------------------
% New output file stuff
%----------------------

% Variable structure

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


