function Out = coldpoolcluster(Idx, B, Grd, varargin)
%COLDPOOLCLUSTER Bottom temp/cold pool clustering exercise
%
% Out = coldpoolcluster(Idx, B, Grd, nxi, neta, plot)
%
% This function runs the cold pool clustering exercise I usually do in prep
% for the Spring PEEC and/or September pre-assessment meetings.  It
% analyzes various metrics related to EBS bottom temperature and cold pool
% extent in the most recent year, and determines which prior years are most
% similar.  
%
% This analysis is a bit art and a bit science, since there are
% plenty of different parameters to play around with in the clustering
% routines.  But these particular metrics have proven useful thus far.
%
% This function is called by the peec_workflow.m script to prep one of the
% slides I've typically presented at the PEEC.
%
% Input variables:
%
%   Idx:    table array with at least the following variables:
%
%           t:          datetimes
%
%           fracbelow0: fraction SEBS below 0 deg C
%
%           fracbelow2: fraction SEBS below 2 deg C 
%
%   B:      structure with the following fields:
%
%           t:      datetimes corresponding to btemp field
%
%           btemp:  nxi x neta x time array of bottom temperature values on
%                   or near July 1 of each year
%
%   Grd:    ROMS grid structure with extended grid fields
%
% Optional input variables (passed as parameter/value pairs)
%
%   plot:   logical scalar, true to create diagnostic plots showing
%           dendrograms and patterns representing each cluster, one figure
%           per metric [flase]
%
%   anom:   logical scalar, true to plot anomalies maps instead of value
%           [false]
%
%   years:  vector of (summer) years to include in cluster analysis
%           [1971:present] 
%
% Output variables:
%
%   Out:    structure with the following fields:
%
%           yr:     1 x nyr array, years within the dataset, with each year
%                   identified by its summer year (i.e. year XXXX runs from
%                   Oct XXXX-1 to Sep XXXX). 
%
%           cp2:    52 x nyr, weekly 2-deg cold pool index
%
%           cp1:    52 x nyr, weekly 1-deg cold pool index
%
%           cp0:    52 x nyr, weekly 2-deg cold pool index
%
%           tmid:   52 x 1 datetime array, weekly time midpoints
%
%           mask:   1 x 3 cell array of ROMS logical masks indicating the
%                   regions used for clustering metrics 2-4
%
%           c:      nyr x 4 array, cluster number corresponding to each
%                   year using 4 different metrics:
%                   1:  0-deg and 2-deg cold pool indices from Oct 1 of the
%                       previous year through Aug 5 of the current one
%                   2:  July 1 bottom temperature across the entire EBS
%                       (Grd.strata_comboeast) 
%                   3:  July 1 bottom temperature in the SEBS region
%                       (strata 10-62) 
%                   4:  July 1 bottom temperature in the NBS region
%                       (non-SEBS comboeast strata) 

% Copyright 2021 Kelly Kearney

p = inputParser;
p.addParameter('plot', false, @(x) validateattributes(x, {'logical'}, {'scalar'}));
p.addParameter('anom', false, @(x) validateattributes(x, {'logical'}, {'scalar'}));
p.addParameter('years', 1971:year(datetime('today')), @(x) validateattributes(x, {'numeric'}, {'integer'}));
p.parse(varargin{:});
Opt = p.Results;

[nxi, neta] = size(Grd.h);

% Selected years

tmask = ismember(year(B.t), Opt.years);
B.t = B.t(tmask);
B.btemp = B.btemp(:,:,tmask);

btemp = reshape(B.btemp, nxi*neta, []);

% Cold pool index timeseries

[~, Aligned2] = romsavgclimatology(Idx.fracbelow2, Idx.t, 'nbin', 52, ...
    'realign', true, 'pivotmonth', 10);
[~, Aligned0] = romsavgclimatology(Idx.fracbelow0, Idx.t, 'nbin', 52, ...
    'realign', true, 'pivotmonth', 10);
cp2 = Aligned2.y';
cp0 = Aligned0.y';
tmid = Aligned2.t(1,1:end-1)'+diff(Aligned2.t(1,:))'./2;
yrcp = year(Aligned2.t(:,end))';


% [cp2, yrcp, tmid] = reshapetimeseries(Idx.t, Idx.fracbelow2, ...
%     'bin', 52, 'pivotdate', [10 1]);
% 
% 
% % [cp1, yrcp, tmid] = reshapetimeseries(Idx.t, Idx.fracbelow1, ...
% %     'bin', 52, 'pivotdate', [10 1]);
% [cp0, yrcp, tmid] = reshapetimeseries(Idx.t, Idx.fracbelow0, ...
%     'bin', 52, 'pivotdate', [10 1]);

seasonmask = days(tmid - tmid(1)) <= 315;
cp2 = cp2(seasonmask,:);
cp0 = cp0(seasonmask,:);
tmid = tmid(seasonmask);

% yrcp = yrcp+1; % label based on summer month, not start

tmask = ismember(yrcp,Opt.years);
cp0 = cp0(:,tmask);
cp2 = cp2(:,tmask);
yrcp = yrcp(tmask);

Out.yr = yrcp;
Out.cp2 = cp2;
% Out.cp1 = cp1;
Out.cp0 = cp0;
Out.tmid = tmid;

for ii = 1:size(cp2,2) % fill NaNs that result from the pivot choice
    isn = isnan(cp2(:,ii));
    cp2(isn,ii) = interp1(tmid(~isn), cp2(~isn,ii), tmid(isn), 'linear');
    
%     isn = isnan(cp1(:,ii));
%     cp1(isn,ii) = interp1(tmid(~isn), cp1(~isn,ii), tmid(isn), 'linear');
    
    isn = isnan(cp0(:,ii));
    cp0(isn,ii) = interp1(tmid(~isn), cp0(~isn,ii), tmid(isn), 'linear');
end

% Clustering metrics

isprior = tmid <= datetime(yrcp(1),8,5);
tmid = tmid(isprior);
cp2 = cp2(isprior,:);
% cp1 = cp1(isprior,:);
cp0 = cp0(isprior,:);

mask{1} = ~isnan(Grd.surveystrata_comboeast) & Grd.mask_rho==1; % EBS
mask{2} = Grd.surveystrata_comboeast <= 62   & Grd.mask_rho==1; % SEBS
mask{3} = Grd.surveystrata_comboeast > 62    & Grd.mask_rho==1; % NEBS

Out.mask = mask;

metric = {...
    [1-cp2; cp0]'  % 0- and 2-deg cold pool index timeseries, Oct-midAug
    btemp(mask{1},:)' % summer EBS temp
    btemp(mask{2},:)' % summer SEBS temp
    btemp(mask{3},:)'}; % summer NEBS temp
ttl = {'1) 0/2-deg Cold pool timeseries', ...
       '2) EBS temp in July', ...
       '3) SEBS temp in July', ...
       '4) NEBS temp in July'};
   
%-----------------   
% Cluster and plot
%-----------------

c = zeros(length(yrcp), length(metric));
for im = 1:length(metric)

    z = linkage(metric{im}, 'ward');

    ctmp = clusterdata(metric{im}, ...
        'linkage', 'ward', 'distance', 'euclidean', 'criterion', ...
        'distance', 'cutoff', 0.2 * max(z(:,3)));
    
    % Cluster group numbers are arbitrary... rescale so based on first year
    % with it (note: this will also keep clusters numbered the same when
    % new years are added if the new year gets added to an existing
    % cluster)  
    
    minyr = splitapply(@(x) min(x,[],1), yrcp', ctmp);
    ctmp = minyr(ctmp);
    
    c(:,im) = findgroups(ctmp);
    
    if Opt.plot
    
        % Set up plot

        h = plotgrid('size', [1 2], 'mar', 0.05);
        h.axr = subgridaxes(h.ax(2), [0.2 0.8], 1); % Dendrogram, timeseries

        setpos(h.fig, '# # 10in 10in');
        set(h.fig, 'color', 'w');

        % Plot dendrogram

        axes(h.axr(1));
        [hh,t,p] = dendrogram(z,0, 'label', compose('%d',yrcp'), ...
            'orientation', 'right', 'colorthreshold', 0.2 * max(z(:,3)));

        cols = dendrocolor(hh, length(yrcp));
        isunq = [true; any(diff(cols,1,1),2)];
        colunq = flipud(cols(isunq,:));
    %     colunq = flipud(unique(cols, 'rows', 'stable')); % unique colors, top to bottom


        yrperm = yrcp(p);

        cplt = c(:,im);
        cplt(isnan(cplt)) = 0;

        [~,yloc] = ismember(yrcp, yrperm);
        cmap = unique([cplt cols(yloc,:)], 'rows');

        plot(h.axr(2), yrcp, cp2(end-5,:), 'k');
        hold(h.axr(2), 'on');
        scatter(h.axr(2), yrcp, cp2(end-5,:), [], cplt, 'filled');
        set(h.axr(2), 'clim', minmax(cplt)+[-0.5 0.5], ...
            'colormap', cmap(:,2:end), 'ylim', [0 1.2], 'box', 'off', ...
            'fontsize', 8);
        title(h.axr(1), ttl{im});

        cmapdetails{im} = cmap;

        set(h.axr(1), 'xaxisloc', 'top', 'fontsize', 8);
        set(h.fig, 'color', 'w');
        ylabel(h.axr(2), '2\circ cold pool index');

        % Plot clusters

        corder = flipud(unique(cplt(p), 'stable'));
        ncluster = length(corder);

        if im == 1 % Line plots for indices

            h.axl = subgridaxes(h.ax(1), ncluster, 2);

            for ii = 1:ncluster
                isin = c(:,im) == corder(ii);
                if ~any(isin)
                    isin = isnan(c(:,im));
                end
                coltmp = cmap(corder(ii)==cmap(:,1),2:4);

                plot(h.axl(ii,1), tmid, cp2(:,~isin), 'color', rgb('light gray'));
                hold(h.axl(ii,1), 'on');
                plot(h.axl(ii,1), tmid, cp2(:,isin), 'color', colunq(ii,:));

                plot(h.axl(ii,2), tmid, cp0(:,~isin), 'color', rgb('light gray'));
                hold(h.axl(ii,2), 'on');
                plot(h.axl(ii,2), tmid, cp0(:,isin), 'color', colunq(ii,:));

            end
            set(h.axl, 'ylim', [0 1], 'xlim', minmax(tmid));
            arrayfun(@(x) set(x.XAxis, 'TickLabelFormat', 'M'), h.axl);


        else % maps

            latlim = minmax(Grd.lat_rho(mask{im-1}));
            lonlim = minmax(Grd.lon_rho(mask{im-1}));

            h.axl = subgridaxes(h.ax(1), ceil(ncluster/2), 2);
            for ii = 1:ncluster

                isin = c(:,im) == corder(ii);
                if ~any(isin)
                    isin = isnan(c(:,im));
                end

                axes(h.axl(ii));
                worldmap(latlim, lonlim);
                if Opt.anom
                    plotromsrho(Grd, nanmean(B.btemp(:,:,isin),3)-nanmean(B.btemp,3), false);
                else
                    plotromsrho(Grd, nanmean(B.btemp(:,:,isin),3), false);
                end
                setm(h.axl(ii), 'fedgecolor', colunq(ii,:), 'flinewidth', 1, 'meridianlabel', 'off', 'parallellabel', 'off');
            end

            if Opt.anom
                set(h.axl, 'clim', [-3 3], 'colormap', cmocean('balance'));
            else
                set(h.axl, 'clim', [-2 3], 'colormap', vivid(cmocean('-dense',5),[0.3 0.7]));
%                 set(h.axl, 'clim', [-2 3], 'colormap', cmocean('-dense',5));
            end         
            
            set(h.axl(ncluster+1:end), 'visible', 'off');
            h.cb = colorbar(h.axl(end,1), 'south');
            setpos(h.cb, '# 0.01 # #');
            set(h.cb, 'tickdir', 'out');

        end
        Out.h{im} = h;
    
    end
    
end

Out.c = c;


