function roms_level2_bottomsurf(varargin)
%ROMS_LEVEL2_BOTTOMSURF Extract bottom/surface average of a variable
%
% roms_level2_bottomsurf(infile, outfolder, outbase, zfilebase, var)
% roms_level2_bottomsurf(infile, outfolder, outbase, zfilebase, var, ...)
%
% This function calculates average surface and/or bottom properties based
% on 3D ROMS output.  It is designed to be used with Level 1-processed ROMS
% output (see roms_level1.m).
%
% Input variables:
%
%   infile:     path to Level 1 or 2 file holding 3D ROMS output
%
%   outfolder:  path to output folder
%
%   outbase:    base name to use for output file.  File will be named
%               <outbase>_<var>_<dz>m.nc (see below for <var> and <dz>
%               input descriptions).  Note that in this case, unlike level
%               1 postprocessing, no file type designation is used in the
%               automatic naming scheme, so this should be included as part
%               of <outbase> if desired.
%
%   zfilebase:  base name for precalculated z-value files (see
%               roms_level2_depths).  This function will look for a file
%               matching the pattern <zfilebase>_z_<grid>.nc, where <grid>
%               is the appropriate horizontal grid for the specified
%               variable (w, u_w, or v_w).
%
%   var:        variable to average.  Should match the primary variable in
%               the specified input file.
%
% Optional input variables (passed as parameter/value pairs):
%
%   varbose:    scalar logical, true to print progress to screen [true]
%
%   dz:         numeric scalar, thickness (m) of the surface or bottom
%               layer used to calculate averages.  Bottom layers are
%               measured from the bottom, and surface layers from the free
%               surface. [5]
%
%   location:   'bottom', 'surface', or 'both', indicating whether to
%               calculate a bottom-layer average, a surface-layer average,
%               or both.

% Copyright 2020 Kelly Kearney

p = inputParser;
p.addRequired('infile',    @(x) validateattributes(x, {'char'}, {'scalartext'}));
p.addRequired('outfolder', @(x) validateattributes(x, {'char'}, {'scalartext'}));
p.addRequired('outbase',   @(x) validateattributes(x, {'char'}, {'scalartext'}));
p.addRequired('zfilebase', @(x) validateattributes(x, {'char'}, {'scalartext'}));
p.addRequired('var',       @(x) validateattributes(x, {'char'}, {'scalartext'}));
p.addParameter('verbose',    true, @(x) validateattributes(x, {'logical'}, {'scalar'}));
p.addParameter('dz',            5, @(x) validateattributes(x, {'numeric'}, {'scalar'}));
p.addParameter('location', 'both', @(x) validateattributes(x, {'char'}, {'scalartext'}));

p.parse(varargin{:});
Opt = p.Results;

validatestring(Opt.location, {'top', 'bottom', 'both'});

% Make sure file exists and has expected dimensions

if ~exist(Opt.infile, 'file')
    error('Cannot find file: %s', Opt.infile);
end

I = ncinfo(Opt.infile, Opt.var);

grdflag = ismember({'xi_rho', 'xi_u', 'xi_v'}, {I.Dimensions.Name}); 
if ~any(grdflag)
    error('Expected grid dimensions not found');
end
if sum(grdflag)>1
    error('More than one grid dimension?');
end

if grdflag(1) % rho variable
    [~,loc] = ismember({I.Dimensions.Name}, {'xi_rho','eta_rho'});
    zfield = 'z_w';
elseif grdflag(2) % u variable
    [~,loc] = ismember({I.Dimensions.Name}, {'xi_u','eta_u'});
    zfield = 'z_u_w';
elseif grdflag(3) % v variable
    [~,loc] = ismember({I.Dimensions.Name}, {'xi_v','eta_v'});
    zfield = 'z_v_w';
end
n = I.Dimensions(loc(1)).Length;
m = I.Dimensions(loc(2)).Length;

zfile = sprintf('%s_%s.nc', Opt.zfilebase, zfield);
if ~exist(zfile, 'file')
    error('Cannot find necessary depths file: %s', zfile);
end

% Time subsetting

nt = ncinfo(Opt.infile, 'ocean_time');
nt = nt.Size;

nchunk = 10;
sta = 1:nchunk:nt;
cnt = min(nchunk, nt-sta+1);

% Read data

[valbot, valtop] = deal(nan(n,m,nt));

if Opt.verbose
    fprintf('Bottom/surface averages, %s...\n', Opt.var);
    c = ConsoleProgressBar;
    c.setMaximum(length(sta));
    c.start();
end

for ii = 1:length(sta)
    
    if Opt.verbose
        c.setValue(ii);
        c.setText(sprintf('%d/%d', ii, length(sta)));
    end
    
    Scs = struct('ocean_time', [sta(ii) cnt(ii) 1]);
 
    Data = ncstruct(Opt.infile, Opt.var, Scs);
    C = ncstruct(zfile, zfield, Scs);
     
    [valbot(:,:,sta(ii):(cnt(ii)+sta(ii)-1)), ...
     valtop(:,:,sta(ii):(cnt(ii)+sta(ii)-1))] = calcmean(C.(zfield), Data.(Opt.var), Opt.dz);
end
if Opt.verbose
    c.stop();
    fprintf('\n');
end

% Create new files

switch Opt.location
    case 'bottom'
        sliceadddata(Opt.infile, Opt.var, permute(valbot, [1 2 4 3]), 'bottom',  Opt.dz, Opt.outfolder, Opt.outbase);
    case 'top'
        sliceadddata(Opt.infile, Opt.var, permute(valtop, [1 2 4 3]), 'surface', Opt.dz, Opt.outfolder, Opt.outbase);
    case 'both'
        sliceadddata(Opt.infile, Opt.var, permute(valbot, [1 2 4 3]), 'bottom',  Opt.dz, Opt.outfolder, Opt.outbase);
        sliceadddata(Opt.infile, Opt.var, permute(valtop, [1 2 4 3]), 'surface', Opt.dz, Opt.outfolder, Opt.outbase);
end


% [moxdir, Grd, nxi, neta] = b10kdata;
% [nu,mu] = size(Grd.lat_u);
% [nv,mv] = size(Grd.lat_v);
% 
% dztopbot = 5;

% % Files
% 
% tfile = '../Postprocessed/average/temp.nc';
% ufile = '../Postprocessed/average/u.nc';
% vfile = '../Postprocessed/average/v.nc';
% cfile = '../Postprocessed/average/depths.nc';
% 
% % Time subsetting
% 
% I = ncinfo(tfile, 'ocean_time');
% nt = I.Size;
% 
% nchunk = 20;
% sta = 1:nchunk:nt;
% cnt = min(nchunk, nt-sta+1);

end

% Subfunctions

function [bot, top] = calcmean(zw, data, dztopbot)

    % Mean, bottom X m
    
    zwrel = min(zw - zw(:,:,1,:), dztopbot);
    dz = diff(zwrel, 1, 3);
    frac = dz./sum(dz,3);
        
    bot = sum(frac .* data, 3);
    
    % Mean, top X m
    
    zwrel = max(zw - zw(:,:,end,:), -dztopbot);
    dz = diff(zwrel, 1, 3);
    frac = dz./sum(dz,3);
        
    top = sum(frac .* data, 3);

end

function sliceadddata(oldfile, var, data, topbot, dz, outfolder, outbase)

    newfile = fullfile(outfolder, sprintf('%s_%s_%s%dm.nc', outbase, var, topbot, dz));
    
    if ~exist(newfile, 'file')
        cmd = sprintf('ncks -F -O -d s_rho,1,1 %s %s', oldfile, newfile);
        system(cmd);
    end

    ncwrite(newfile, var, data);
    lname = ncreadatt(oldfile, var, 'long_name');
    ncwriteatt(newfile, var, 'long_name', sprintf('%s, %s %dm mean', lname, topbot, dz));  

    ncaddhis(newfile, sprintf('%s average calculated via roms_level2_bottomsurf.m', [upper(topbot(1)) topbot(2:end)]));
    
    cmd = sprintf('ncwa -O -a s_rho %s %s', newfile, newfile);
    system(cmd);
    
end


