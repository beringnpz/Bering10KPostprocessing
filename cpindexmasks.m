function M = cpindexmasks(Grd)
%CPINDEXMASKS Returns mask data used to calculate cold pool indices 
%
% M = cpindexmasks
% M = cpindexmasks(Grd)
%
% This small wrapper function coordinates the masks used to calculate cold
% pool indices across a number of scripts.
%
% Input variables:
% 
%   Grd:    ROMS grid file structure, with added survey strata fields
%           (currently custom to Bering10K domain).  If not included,
%           output will only include the name field.
%
% Output variables:
%
%   M:      structure with the following fields:
%
%           mask:   nxi x neta x 4 logical array, indicating which model
%                   grid cells fall within each masking polygon
%
%           name:   4 x 1 cell array of character arrays, names associated
%                   with each mask

% Copyright 2024 Kelly Kearney

M.name = {...
    'SEBS-northwest'
    'SEBS+northwest'
    'SEBS-northwest, model shelf'
    'SEBS+northwest, model shelf'};

if nargin > 1
    M.mask = cell(4,1);
    M.mask{1} = Grd.surveystrata_comboeast >= 10 & Grd.surveystrata_comboeast <= 62; % SEBS, no NW extension
    M.mask{2} = M.mask{1} | ismember(Grd.surveystrata_comboeast, [82 90]); % SEBS+NW
    M.mask{3} = M.mask{1} & Grd.h <= 200; % SEBS, no shelf mismatch
    M.mask{4} = M.mask{2} & Grd.h <= 200; % SEBS+NW, no shelf mismatch
    M.mask = cat(3, M.mask{:});
end
