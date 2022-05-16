function R = raw_variable_summary

moxdir = b10kdata;
b10koutfol = fullfile(moxdir, 'bering10k', 'output');

simdetails = {
'B10K-K20P19_CMIP6_cesm_historical'  , fullfile(b10koutfol, 'forecasts/cmip6/cesm_historical/Out')        ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
'B10K-K20P19_CMIP6_cesm_ssp126'      , fullfile(b10koutfol, 'forecasts/cmip6/cesm_ssp126/Out')            ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
'B10K-K20P19_CMIP6_cesm_ssp585'      , fullfile(b10koutfol, 'forecasts/cmip6/cesm_ssp585/Out')            ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
'B10K-K20P19_CMIP6_gfdl_historical'  , fullfile(b10koutfol, 'forecasts/cmip6/gfdl_historical/Out')        ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
'B10K-K20P19_CMIP6_gfdl_ssp126'      , fullfile(b10koutfol, 'forecasts/cmip6/gfdl_ssp126/Out')            ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
'B10K-K20P19_CMIP6_gfdl_ssp585'      , fullfile(b10koutfol, 'forecasts/cmip6/gfdl_ssp585/Out')            ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
'B10K-K20P19_CMIP6_miroc_historical' , fullfile(b10koutfol, 'forecasts/cmip6/miroc_historical/Out')       ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
'B10K-K20P19_CMIP6_miroc_ssp126'     , fullfile(b10koutfol, 'forecasts/cmip6/miroc_ssp126/Out')           ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
'B10K-K20P19_CMIP6_miroc_ssp585'     , fullfile(b10koutfol, 'forecasts/cmip6/miroc_ssp585/Out')           ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
'B10K-K20_CORECFS'                   , fullfile(b10koutfol, 'hindcasts/npz_201904_delta30/Out')           ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
'B10K-H16_CORECFS'                   , fullfile(b10koutfol, 'hindcasts/npz_201904_aclim/Out/')            ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
'B10K-K20P19_CMIP5_CESM_rcp45'       , fullfile(b10koutfol, 'cmip5-avg-carbon/CESM_rcp45/')               ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
'B10K-K20P19_CMIP5_CESM_rcp85'       , fullfile(b10koutfol, 'cmip5-avg-carbon/CESM_rcp85/')               ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
'B10K-K20P19_CMIP5_GFDL_rcp45'       , fullfile(b10koutfol, 'cmip5-avg-carbon/GFDL_rcp45/')               ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
'B10K-K20P19_CMIP5_GFDL_rcp85'       , fullfile(b10koutfol, 'cmip5-avg-carbon/GFDL_rcp85/')               ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
'B10K-K20P19_CMIP5_MIROC_rcp45'      , fullfile(b10koutfol, 'cmip5-avg-carbon/MIROC_rcp45/')              ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
'B10K-K20P19_CMIP5_MIROC_rcp85'      , fullfile(b10koutfol, 'cmip5-avg-carbon/MIROC_rcp85/')              ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
'B10K-H16_CMIP5_CESM_BIO_rcp85'      , fullfile(b10koutfol, 'cmip5-avg/CESM_BIO_rcp85/')                  ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
'B10K-H16_CMIP5_CESM_rcp45'          , fullfile(b10koutfol, 'cmip5-avg/CESM_rcp45/')                      ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
'B10K-H16_CMIP5_CESM_rcp85'          , fullfile(b10koutfol, 'cmip5-avg/CESM_rcp85/')                      ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
'B10K-H16_CMIP5_GFDL_BIO_rcp85'      , fullfile(b10koutfol, 'cmip5-avg/GFDL_BIO_rcp85/')                  ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
'B10K-H16_CMIP5_GFDL_rcp45'          , fullfile(b10koutfol, 'cmip5-avg/GFDL_rcp45/')                      ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
'B10K-H16_CMIP5_GFDL_rcp85'          , fullfile(b10koutfol, 'cmip5-avg/GFDL_rcp85/')                      ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
'B10K-H16_CMIP5_MIROC_rcp45'         , fullfile(b10koutfol, 'cmip5-avg/MIROC_rcp45/')                     ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
'B10K-H16_CMIP5_MIROC_rcp85'         , fullfile(b10koutfol, 'cmip5-avg/MIROC_rcp85/')                     ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc'            
'B10K-K20nobio_CFS'                  , fullfile(b10koutfol, 'hindcasts/phys_201910_30_singlesource/Out/') ,'physhccfs*_avg_*.nc'   ,'physhccfs*_his_*.nc'  ,'physhccfs*_sta_*.nc' 
'B10K-K20nobio_CORE'                 , fullfile(b10koutfol, 'hindcasts/phys_201910_30_singlesource/Out/') ,'physhccore*_avg_*.nc'  ,'physhccore*_his_*.nc' ,'physhccore*_sta_*.nc'
'B10K-K20P19_CORECFS'                , '/gscratch/jisao/pilchd/Bering10k/output/hindcast_aug2021/OUT/'    ,'*avg*.nc'              ,'*his*.nc'             ,'*sta*.nc' 
};

simdetails = cell2table(simdetails, 'VariableNames', {'base','path', 'afilter', 'hfilter', 'sfilter'});
nsim = height(simdetails);

for isim = 1:nsim
    simdetails.afiles{isim} = dirfull(fullfile(simdetails.path{isim}, simdetails.afilter{isim}));
    simdetails.hfiles{isim} = dirfull(fullfile(simdetails.path{isim}, simdetails.hfilter{isim}));
    simdetails.sfiles{isim} = dirfull(fullfile(simdetails.path{isim}, simdetails.sfilter{isim}));
end

% New postprocessing parent folders

simdetails.ppfol = fullfile(moxdir, 'roms_for_public', simdetails.base);

% Scan the raw files

for is = nsim:-1:1
    if ~isempty(simdetails.afiles{is})
        S(is,1) = simvariableinfo({simdetails.afiles{is}.name}, 'average');
    end
    if ~isempty(simdetails.hfiles{is})
        S(is,2) = simvariableinfo({simdetails.hfiles{is}.name}, 'history');
    end
    if ~isempty(simdetails.sfiles{is})
        S(is,3) = simvariableinfo({simdetails.sfiles{is}.name}, 'station');
    end
end

% [tf,loc] = ismember(V.sname, simdetails.base);
% S = S(loc);

refidx = find(strcmp(simdetails.base, 'B10K-K20_CORECFS')); % for long name and unit clashes
for is = 1:length(S)
    S(is,1).priority = repmat(is==refidx, length(S(is,1).name), 1);
    S(is,2).priority = repmat(is==refidx, length(S(is,2).name), 1);
    S(is,3).priority = repmat(is==refidx, length(S(is,3).name), 1);
end

% Parse out unique variables and how they appear across types of files
% Note: if a variable is only available in the station files, we can't
% automatically tell what horizontal grid is falls on.

Sall = arrayfun(@(x) struct2table(x), S, 'uni', 0);

Sah = cat(1, Sall{:,1:2});
Sah = sortrows(Sah, 'priority', 'descend');
Sah = unique(removevars(Sah, 'priority'));

Ss = cat(1, Sall{:,3});
Ss = sortrows(Ss, 'priority', 'descend');
Ss = unique(removevars(Ss, 'priority'));

R.vname = unique([Sah.name; Ss.name]);
R.sname = simdetails.base;
R.tname = {'average', 'history', 'station'};

[tf1,loc1] = ismember(R.vname, Sah.name);
[tf2,loc2] = ismember(R.vname, Ss.name);

nv = length(R.vname);
R.lname = cell(nv,1);
R.grid_ah = cell(nv,1);
R.grid_s = cell(nv,1);

for iv = 1:length(R.vname)
    if tf1(iv)
        R.grid_ah{iv} = Sah.grid{loc1(iv)};
        R.lname{iv} = Sah.lname{loc1(iv)};
        R.units{iv} = Sah.units{loc1(iv)};
    end
    if tf2(iv)
        R.grid_s{iv} = Ss.grid{loc2(iv)};
        if ~tf1(iv)
            R.lname{iv} = Ss.lname{loc2(iv)};
            R.units{iv} = Ss.units{loc2(iv)};
        end
    end
end

% Mark which variables were saved for each simulation

R.varsaved = false(nv,nsim,3);
for is = 1:nsim
    for it = 1:3
        if ~isempty(S(is,it).name)
            R.varsaved(:,is,it) = ismember(R.vname, S(is,it).name);
        end
    end
end




