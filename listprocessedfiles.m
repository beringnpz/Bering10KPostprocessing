function Ftbl = listprocessedfiles(ppfolder)
%LISTPROCESSEDFILES Lists all post-processed files
%
% This function lists all the netCDF files it finds in the indicated folder
% under a LevelX subfolder.  It also parses out some relevant information
% from the file name.
%
% Input variables: 
%
%   ppfolder:   post-processing parent folder.  If not included, defaults
%               to the hyak-mox /gscratch/bumblereem/roms_for_public folder
%
% Output variables:
%
%   Ftbl:       table with file info.  Includes all fields associated with
%               a call to the dir function, as well as the following
%   
%               level:      processing level (1, 2, or 3)
%
%               simulation: name of simulation
%
%               type:       type of file (average, history, station, or
%                           other)
%
%               variable:   primary variable in file ('' if not
%                           relevant)
%
%               years:      year span in file ('' if not processed and
%                           named with time blocks)

% Copyright 2021 Kelly Kearney

if nargin<1
    moxdir = b10kdata;
    ppfolder = fullfile(moxdir, 'roms_for_public');
end

% List of simulations that have been processed

F = dir(ppfolder);
issim = [F.isdir] & ~startsWith({F.name}, '.');

simnames = {F(issim).name}';

simfolder = fullfile(ppfolder, simnames);

% Gather list of post-processed files

fol = unique(simfolder);
tmp = cell(size(fol));
for ii = 1:length(fol)
    tmp{ii} = dir(fullfile(fol{ii}, 'Lev*','*.nc'));
end
npersim = cellfun(@length, tmp);

tmp = cat(1, tmp{:});
Ftbl = struct2table(tmp);
nfile = height(Ftbl);

% Mark simulation and level

[pth, levfol] = fileparts(Ftbl.folder);
level = str2double(strrep(levfol, 'Level', ''));
[pth, simfol] = fileparts(pth);

Ftbl.level = level;
Ftbl.simulation = simfol;

% Mark each file as average, history, station, or other (other includes
% Level 3 stuff)

filters = {...
    contains(Ftbl.name, '_average_')
    contains(Ftbl.name, '_history_')
    contains(Ftbl.name, '_station_')
    };

type = cell(nfile, 1);
[type{:}] = deal('other');
[type{filters{1}}] = deal('average');
[type{filters{2}}] = deal('history');
[type{filters{3}}] = deal('station');

Ftbl.type = type;

% Mark variable, if relevant

tok = regexp(Ftbl.name, '(average|history|station)_(\w*)(\.nc)', 'tokens', 'once');
isin = ~cellfun(@isempty, tok);
tok = cat(1, tok{:});

var = cell(nfile,1);
[var{~isin}] = deal('');
var(isin) = tok(:,2);

Ftbl.variable = var;

% Mark time block, if relevant

tok = regexp(Ftbl.name, '_(\d\d\d\d-\d\d\d\d)_', 'tokens', 'once');
isin = ~cellfun(@isempty, tok);
tok = cat(1, tok{:});

yrs = cell(nfile,1);
[yrs{~isin}] = deal('');
yrs(isin) = tok;

Ftbl.years = yrs;

