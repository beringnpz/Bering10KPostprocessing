function roms_level2_depthintegrate(varargin)
%ROMS_LEVEL2_DEPTHINTEGRATE Integrate over depth for a ROMS variable
%
% roms_level2_depthintegrate(infile, outfolder, outbase, zfilebase, var)
% roms_level2_depthintegrate(infile, outfolder, outbase, zfilebase, var, ...)
%
% This function calculates depth-integrated properties based on 3D ROMS
% output.  It is designed to be used with Level 1-processed ROMS output
% (see roms_level1.m) or files with similar format and naming scheme. 
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
%   var:        variable to integrate. Should match the primary variable in
%               the specified input file.
%
% Optional input variables (passed as parameter/value pairs):
%
%   verbose:    scalar logical, true to print progress to screen [true]
%
%   simplesum:  scalar logical, true if original file values already
%               represent integrated-over-layer values, and therefore
%               layers just need to be summed together.  If true, zfilebase
%               input is ignored.

% Copyright 2020 Kelly Kearney

p = inputParser;
p.addRequired('infile',    @(x) validateattributes(x, {'char'}, {'scalartext'}));
p.addRequired('outfolder', @(x) validateattributes(x, {'char'}, {'scalartext'}));
p.addRequired('outbase',   @(x) validateattributes(x, {'char'}, {'scalartext'}));
p.addRequired('zfilebase', @(x) validateattributes(x, {'char'}, {'scalartext'}));
p.addRequired('var',       @(x) validateattributes(x, {'char'}, {'scalartext'}));
p.addParameter('verbose',    true, @(x) validateattributes(x, {'logical'}, {'scalar'}));
p.addParameter('simplesum',  false, @(x) validateattributes(x, {'logical'}, {'scalar'}));


p.parse(varargin{:});
Opt = p.Results;

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

if ~Opt.simplesum
    if contains(Opt.zfilebase, 'diagnos')
        zfile = sprintf('%s_%s.nc', strrep(Opt.zfilebase, '_diagnos_', '_average_'), zfield);
    else
        zfile = sprintf('%s_%s.nc', Opt.zfilebase, zfield);
    end
    if ~exist(zfile, 'file')
        error('Cannot find necessary depths file: %s', zfile);
    end
end

% Time subsetting

nt = ncinfo(Opt.infile, 'ocean_time');
nt = nt.Size;

nchunk = 10;
sta = 1:nchunk:nt;
cnt = min(nchunk, nt-sta+1);

% Read data

valint = nan(n,m,nt);

if Opt.verbose
    fprintf('Depth integration...%s\n', Opt.var);
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
    if Opt.simplesum
        valint(:,:,sta(ii):(cnt(ii)+sta(ii)-1)) = sum(Data.(Opt.var),3);
    else
        C = ncstruct(zfile, zfield, Scs);
        valint(:,:,sta(ii):(cnt(ii)+sta(ii)-1)) = calcsum(C.(zfield), Data.(Opt.var));
    end
end
if Opt.verbose
    c.stop();
    fprintf('\n');
end

% Create new files

sliceadddata(Opt.infile, Opt.var, permute(valint, [1 2 4 3]), Opt.outfolder, Opt.outbase, Opt.simplesum);

end

% Subfunctions

function val = calcsum(zw, data)
    dz = diff(zw, 1, 3);
    val = sum(data .* dz, 3);
end

function sliceadddata(oldfile, var, data, outfolder, outbase, sumflag)

    newfile = fullfile(outfolder, sprintf('%s_%s_integrated.nc', outbase, var));
    
    if ~exist(newfile, 'file')
        cmd = sprintf('ncks -F -O -d s_rho,1,1 %s %s', oldfile, newfile);
        system(cmd);
    end

    ncwrite(newfile, var, data);
    
    I  = ncinfo(oldfile, var);
    
    if ismember('long_name', {I.Attributes.Name})
        lname = ncreadatt(oldfile, var, 'long_name');
        if sumflag
            ncwriteatt(newfile, var, 'long_name', sprintf('%s, summed over depth', lname)); 
        else
            ncwriteatt(newfile, var, 'long_name', sprintf('%s, integrated over depth', lname));  
        end
    end
    if ~sumflag && ismember('units', {I.Attributes.Name})
        unit = ncreadatt(oldfile, var, 'units');
        ncwriteatt(newfile, var, 'units', sprintf('(%s)*m', unit));  
    end
    
    ncaddhis(newfile, 'integrated value calculated via roms_level2_depthintegrate.m');
    
    cmd = sprintf('ncwa -O -a s_rho %s %s', newfile, newfile);
    system(cmd);
    
end


