function template = akfinwrite(varargin)
%AKFINWRITE Create AKFIN Database-formatted ESP file
%
% text = akfinwrite(param, val, ...)
%
% Input variables (passed as parameter/value pairs):
%
%   submission_year:    Current year of contribution submission 
%
%   indicator_name:     Composite key (meaning this must be unique to the
%                       indicator) based on the ESP naming convention and
%                       used for joining ESP data tables.  
%
%   description:        Brief description of the indicator and why it is
%                       important to groundfish fishery management. Please
%                       make sure this description includes information on
%                       the spatial distribution of the indicator and how
%                       the data for the indicator are collected. The
%                       spatial resolution can be a cell size dimension
%                       (e.g., 5 km resolution for gridded data) or area of
%                       sampling for a survey (e.g., Shelikof Strait). The
%                       data collection method can be brief (e.g., survey
%                       name and gear type, satellite sensor and version,
%                       stock assessment model output, fishery observer
%                       data, community reports, etc.) and can include a
%                       reference to methods detailed elswhere. (Limit to
%                       4000 characters)    
%
%   status_trends:      Information on the current status of the indicator
%                       in the context of historical trends. This is
%                       similar to the Ecosystem Status Report contribution
%                       (Limit to 4000 characters)   
%
%   factors:            Information on the potential causes for observed
%                       trends and current status (focus on most recent
%                       year). This is similar to the Ecosystem Status
%                       Report contribution. (Limit to 4000 characters)  
%   
%   implications:       Information that briefly answers these questions:
%                       What are the implications or impacts of the
%                       observed trends on the ecosystem or ecosystem
%                       components? What do the trends mean? Why are they
%                       important? How can this information be used to
%                       inform groundfish management decisions? This is
%                       similar to the Ecosystem Status Report
%                       contribution. (Limit to 4000 characters)
%
%   references:         Include any full references that are associated
%                       with the indicator. This may include data
%                       references such as from an ERDDAP webpage or
%                       literature cited (plus the DOI where possible).
%                       Please also provide the full reference if there is
%                       an associated Ecosystem Status Report contribution
%                       that provides more details on the indicator.
%
%   year:               List of years for the indicator contribution
%
%   data:               List of data values for the indicator contribuion
%                       (must match the year list length) 
%
%   file:               string,  file name where text will be written
%
% Output variables:
%
%   text:               string array, text of AKFIN-formatted file (will
%                       also be written to file if the file input is not
%                       empty)  

% Copyright 2023 Kelly Kearney

p = inputParser;
p.addParameter('submission_year',[], @(x) validateattributes(x, {'numeric'}, {'scalar','integer'}));
p.addParameter('indicator_name', '', @(x) validateattributes(x, {'string','char'}, {'scalartext'}));
p.addParameter('description',    '', @(x) validateattributes(x, {'string','char'}, {'scalartext'}));
p.addParameter('status_trends',  '', @(x) validateattributes(x, {'string','char'}, {'scalartext'}));
p.addParameter('factors',        '', @(x) validateattributes(x, {'string','char'}, {'scalartext'}));
p.addParameter('implications',   '', @(x) validateattributes(x, {'string','char'}, {'scalartext'}));
p.addParameter('references',     '', @(x) validateattributes(x, {'string','char'}, {'scalartext'}));
p.addParameter('year',           [], @(x) validateattributes(x, {'numeric'}, {'integer', 'vector'}));
p.addParameter('data',           [], @(x) validateattributes(x, {'numeric'}, {'vector'}));
p.addParameter('file',           '', @(x) validateattributes(x, {'string','char'}, {'scalartext'}));

p.parse(varargin{:});
Opt = p.Results;

% Fill in template with user details

template = [...
"#Ecosystem and Socioeconomic Profile (ESP) indicator contribution for stocks managed under the North Pacific Fisheries Management Council"
"#This template is required for updating ESP indicator contribution information"
"#There are two required sections to check or update (see below): Indicator Review and Indicator Data"
"#Please fill in the text (surrounded by "" "") or data as values in the line after each field marked with a # and capitalized name (e.g., #INDICATOR_NAME, the next line should be the name of your indicator, for example ""Annual_Arrowtooth_Biomass_GOA_Model"")"
"#Note that all fields are described in the Alaska ESP User Guide, please see pdf attached to indicator assignment email for more details"
"#INDICATOR_REVIEW ----------------------------------------------------------------------------------------"
"#SUBMISSION_YEAR - Current year of contribution submission"
string(Opt.submission_year)
"#INDICATOR_NAME - Composite key (meaning this must be unique to the indicator) based on the ESP naming convention and used for joining ESP data tables. Please see email with your indicator names, and copy/paste name to this location. Note: this name must match the ESP records provided in the email, please do not change. Questions, contact kalei.shotwell@noaa.gov"
"""" + Opt.indicator_name + """"
"#DESCRIPTION - Brief description of the indicator and why it is important to groundfish fishery management. Please make sure this description includes information on the spatial distribution of the indicator and how the data for the indicator are collected. The spatial resolution can be a cell size dimension (e.g., 5 km resolution for gridded data) or area of sampling for a survey (e.g., Shelikof Strait). The data collection method can be brief (e.g., survey name and gear type, satellite sensor and version, stock assessment model output, fishery observer data, community reports, etc.) and can include a reference to methods detailed elswhere. (Limit to 4000 characters)"  
"""" + Opt.description + """"
"#STATUS_TRENDS - Information on the current status of the indicator in the context of historical trends. This is similar to the Ecosystem Status Report contribution (Limit to 4000 characters)"
"""" + Opt.status_trends + """"
"#FACTORS - Information on the potential causes for observed trends and current status (focus on most recent year). This is similar to the Ecosystem Status Report contribution. (Limit to 4000 characters)"
"""" + Opt.factors + """"
"#IMPLICATIONS - Information that briefly answers these questions: What are the implications or impacts of the observed trends on the ecosystem or ecosystem components? What do the trends mean? Why are they important? How can this information be used to inform groundfish management decisions? This is similar to the Ecosystem Status Report contribution. (Limit to 4000 characters)"
"""" + Opt.implications + """"
"#REFERENCES - Include any full references that are associated with the indicator. This may include data references such as from an ERDDAP webpage or literature cited (plus the DOI where possible). Please also provide the full reference if there is an associated Ecosystem Status Report contribution that provides more details on the indicator."
"""" + Opt.references + """"
"#INDICATOR_DATA ----------------------------------------------------------------------------------------------"
"#YEAR - List of years for the indicator contribution"
join(compose("%d",Opt.year), " ")
"#INDICATOR_VALUE - List of data values for the indicator contribuion (must match the year list length)"
join(compose("%g",Opt.data), " ")
];

% Write to file

if ~isempty(Opt.file)
    if exist(Opt.file, 'file')
        warning('File %s already exists, exiting without writing to file', Opt.file);
    else
        fid = fopen(Opt.file, 'wt');
        fprintf(fid, '%s\n', template{:});
        fclose(fid);
    end
end
