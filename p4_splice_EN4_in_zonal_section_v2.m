% p4_splice_EN4_in_zonal_section_v2
% SJones 04/21

% Add EN4 data across 47N to existing boundary dataset.  Recompute
% along-boundary distances.

clear; close all;

addpath(genpath('D:\Work_computer_sync\MATLAB_functions'));

minlon = -65;
maxlon = 5;
minlat = 46;
maxlat = 67;



%% Gridding parameters 
pres_grid = 20:20:5500;
% EN4_pres_threshold = 1900; % What depth to start splicing in EN4 data

GRID.season_description = ['JFM'; 'AMJ'; 'JAS' ; 'OND'];
season_months = [1 2 3; 4 5 6; 7 8 9; 10 11 12];

%% Load boundary dataset

load D:\Work_computer_sync\OSNAP_postdoc\PAPERS_NEW\N_Atlantic_boundary\matlab\V3_050321\intermediate_saves\contour_data_1000
load D:\Work_computer_sync\OSNAP_postdoc\PAPERS_NEW\N_Atlantic_boundary\matlab\V3_050321\intermediate_saves\pressure_gridded_1000

%% Load EN4, reduce to zonal transect
EN4 = load('E:\EN4\EN4_SPG_2000_2020.mat');
% Throw out all except 47 N (EN4.lat(3))
EN4.CT = squeeze(EN4.CT(:,3,:,:));
EN4.S = squeeze(EN4.S(:,3,:,:));
EN4.SA = squeeze(EN4.SA(:,3,:,:));
EN4.S_weight = squeeze(EN4.S_weight(:,3,:,:));
EN4.T_weight = squeeze(EN4.T_weight(:,3,:,:));
EN4.Tptemp = squeeze(EN4.Tptemp(:,3,:,:));
EN4.p = squeeze(EN4.p(:,3,:));
EN4.sigma0all = squeeze(EN4.sigma0(:,3,:,:));

EN4.month = nan*ones(length(EN4.time),1);
for aa = 1:length(EN4.time)
    EN4.month(aa,1) = str2num(datestr(EN4.time(aa),'mm'));
    % EN4.year(1,aa) = str2num(datestr(EN4.time(aa),'yyyy'));
end

% Average by season
EN4S.CT = nan*ones(length(EN4.lon),length(EN4.depth),4);
EN4S.SA = EN4S.CT;
EN4S.sigma0 = EN4S.CT;
EN4S.lon = EN4.lon;
EN4S.p = EN4.p(1,:)';

for season = 1:4
    ind = find(EN4.month == season_months(season,1) | EN4.month == season_months(season,2) | EN4.month == season_months(season,3));
    
    EN4S.CT(:,:,season) = nanmean(EN4.CT(:,:,ind),3);
    EN4S.SA(:,:,season) = nanmean(EN4.SA(:,:,ind),3);
    EN4S.sigma0(:,:,season) = nanmean(EN4.sigma0all(:,:,ind),3);
end

%% Prepare input and output arrays for interpolation
% in this version just extract grid cells in the zonal lon range
lon_start = cont_lon(end);
lon_end = cont_lon(1);
ind = find(EN4S.lon > lon_start & EN4S.lon < lon_end);
EN4S.CT = EN4S.CT(ind,:,:);
EN4S.SA = EN4S.SA(ind,:,:);
EN4S.sigma0 = EN4S.sigma0(ind,:,:);
EN4S.lon = EN4S.lon(ind);

% Expand lat_grid and lon_grid to capture EN4 extension
transect_mask_grid = 0 * lon_grid;
lon_grid = [lon_grid EN4S.lon'];
lat_grid = [lat_grid repmat(47,1,length(EN4S.lon))];
dist_grid = [dist_grid nan*ones(1,length(EN4S.lon))];
dist_grid_accurate = [dist_grid_accurate nan*ones(1,length(EN4S.lon))];
transect_mask_grid = [transect_mask_grid repmat(1,1,length(EN4S.lon))];
GRID.SA = cat(2,GRID.SA,nan*ones(length(pres_grid),length(EN4S.lon),4));
GRID.CT = cat(2,GRID.CT,nan*ones(length(pres_grid),length(EN4S.lon),4));
GRID.dens = cat(2,GRID.dens,nan*ones(length(pres_grid),length(EN4S.lon),4));

% also expand contour data
zonal_transect_mask = zeros(1,length(cont_lon));
cont_lon = [cont_lon EN4S.lon'];
cont_lat = [cont_lat repmat(47,1,length(EN4S.lon))];
zonal_transect_mask = [zonal_transect_mask repmat(1,1,length(EN4S.lon))];
cont_distance = [cont_distance nan*ones(1,length(EN4S.lon))];

% preallocate
CT_insert = nan*ones(length(pres_grid),length(EN4S.lon),4);
SA_insert = CT_insert;



for season = 1:4 % for each season
    
    %% perform interpolation
    
    for aa = 1:length(EN4S.lon) % for each profile
        % Vq = interp1(X,V,Xq)
        CT_insert(:,aa,season) = interp1(EN4S.p,squeeze(EN4S.CT(aa,:,season))',pres_grid);
        SA_insert(:,aa,season) = interp1(EN4S.p,squeeze(EN4S.SA(aa,:,season))',pres_grid);
    end

    % sigma0 = gsw_sigma0(SA,CT)
    dens_insert = gsw_sigma0(SA_insert,CT_insert);
    
    %% Now splice in EN4 profiles
    ind = find(transect_mask_grid == 1);
    GRID.SA(:,ind,season) = SA_insert(:,:,season);
    GRID.CT(:,ind,season) = CT_insert(:,:,season);
    GRID.dens(:,ind,season) = dens_insert(:,:,season);

end % End season loop

% Add on bathy from EN4
startpoint = length(bathy_depth);
for aa = 1:length(EN4S.lon)
    maxdepthind = max(find(~isnan(EN4S.CT(aa,:,1))));
    bathy_depth(1,aa+startpoint) = EN4.depth(maxdepthind);
    bathy_pres(1,aa+startpoint) = EN4S.p(maxdepthind);
end


% %% Trim profiles based on pressure of bed 
% 
% for ee = 1:length(bathy_depth) % for each profile
%     ind = find(pres_grid > bathy_pres(ee));
%     GRID.SA(ind,ee,:) = nan;
%     GRID.CT(ind,ee,:) = nan;
%     GRID.dens(ind,ee,:) = nan;
% end

% Recalculate cont_distance
ind = find(zonal_transect_mask == 1);
for aa = 1:length(ind)
    lon1 = cont_lon(ind(aa)-1); lon2 = cont_lon(ind(aa));
    lat1 = cont_lat(ind(aa)-1); lat2 = cont_lat(ind(aa));
    lon_in = [lon1 lon2]; lat_in = [lat1 lat2];
    cont_distance(1,ind(aa)-1) = (gsw_distance(lon_in,lat_in))/1000;
end

% And loop back round to first location
lon1 = cont_lon(end); lon2 = cont_lon(1);
lat1 = cont_lat(end); lat2 = cont_lat(1);
lon_in = [lon1 lon2]; lat_in = [lat1 lat2];
cont_distance(1,end) = (gsw_distance(lon_in,lat_in))/1000;

% Redo summed distance.
cont_dist_cumul = cumsum(cont_distance);

%insert these values into the gridded distance variable too...
contind = find(zonal_transect_mask == 1);
gridind = find(transect_mask_grid == 1);

% Need to add a nasty 
% little fudge to account for the fact that dist_grid has the same number of cells 
% as lon_grid etc, not n-1 due to the way it was calculated in p1. 
dist_grid(1,gridind(1:end-1)) = cont_dist_cumul(contind(1:end-1));
dx = dist_grid(end-1) - dist_grid(end-2); 
dist_grid(1,gridind(end)) = dist_grid(end-1) + dx;

dist_grid_accurate(1,gridind(1:end-1)) = cont_dist_cumul(contind(1:end-1));
dx = dist_grid_accurate(end-1) - dist_grid_accurate(end-2); 
dist_grid_accurate(1,gridind(end)) = dist_grid_accurate(end-1) + dx;

%% Add the final distance value on manually - last cell to first cell.
dx2 = cont_distance(end) + dist_grid(1);
dist_grid = [dist_grid dist_grid(end)+dx2];

dx2 = cont_distance(end) + dist_grid_accurate(1);
dist_grid_accurate = [dist_grid_accurate dist_grid_accurate(end)+dx2];

lon_grid = [lon_grid lon_grid(1)];
lat_grid = [lat_grid lat_grid(1)];
transect_mask_grid = [transect_mask_grid transect_mask_grid(1)];
bathy_depth = [bathy_depth bathy_depth(1)];
bathy_pres = [bathy_pres bathy_pres(1)];
% Concatenate GRID variables as well
GRID.SA = cat(2,GRID.SA,GRID.SA(:,1,:));
GRID.CT = cat(2,GRID.CT,GRID.CT(:,1,:));
GRID.dens = cat(2,GRID.dens,GRID.dens(:,1,:));



%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Fig 3: T,S, dens contour plots

figure(3);
for season = 1:4
subplot(4,1,season);
hold on;
% dist_grid,CT_wint
levels = [0:0.5:14];
caxis([2 14]);
contourf(dist_grid,pres_grid,GRID.CT(:,:,season),levels); shading flat;
% plot(dist_grid,bathy_pres,'k','linewidth',2); 

set(gca,'ydir','reverse');
%ylim([0 4800])
ylim([0 1000])
ylabel('Pres (db)','fontsize',14);
xlabel('Distance along contour (km)','fontsize',14);
title(GRID.season_description(season,:),'fontsize',14)
set(gca,'fontsize',14);
colours = (lbmap(24,'RedBlue'));
colormap(flipud(colours));
cb = colorbar;
ylabel(cb,'Conservative temperature (^oC)','fontsize',14);
end

%% print figure
width  = 1500;  % frame width
height = 2000;  % frame height
pngname = ('plots/p3_CT_EN4added');

% set background color (outside axes)
set(gcf,'color',[1 1 1]);

% don't change background color when printing
set(gcf,'inverthardcopy','off');

% set size of frame to be written
resolution=150;
set(gcf,'paperunits','inches');
set(gcf,'paperposition',[0 0 width height]./resolution);

% write .png file
% the 'zbuffer' method is likely to look similar to the plot window
print('-dpng', ['-r' num2str(resolution)], '-opengl', pngname);


%% Sal
figure(4)
for season = 1:4
subplot(4,1,season);
hold on;
% dist_grid,CT_wint
levels = [33:0.05:37];
caxis([34.5 35.7]);
%pcolor(dist_grid,-depth_grid,SA_wint); shading flat;
contourf(dist_grid,pres_grid,GRID.SA(:,:,season),levels); shading flat;
%plot(dist_grid,bathy_pres,'k','linewidth',2); 

set(gca,'ydir','reverse');
%ylim([0 4800])
 ylim([0 1000])
ylabel('Pres (db)','fontsize',14);
xlabel('Distance along contour (km)','fontsize',14);
title(GRID.season_description(season,:),'fontsize',14)
set(gca,'fontsize',14);
colours = (lbmap(24,'RedBlue'));
colormap(flipud(colours));
cb = colorbar;
ylabel(cb,'Absolute salinity (g/kg)','fontsize',14);
end

%% print figure
width  = 1500;  % frame width
height = 2000;  % frame height
pngname = ('plots/p4_SA_EN4added');

% set background color (outside axes)
set(gcf,'color',[1 1 1]);

% don't change background color when printing
set(gcf,'inverthardcopy','off');

% set size of frame to be written
resolution=150;
set(gcf,'paperunits','inches');
set(gcf,'paperposition',[0 0 width height]./resolution);

% write .png file
% the 'zbuffer' method is likely to look similar to the plot window
print('-dpng', ['-r' num2str(resolution)], '-opengl', pngname);

%% Dens
figure(5)
for season = 1:4
subplot(4,1,season);
hold on;
% dist_grid,CT_wint
levels = [25:0.05:28];
caxis([26.8 27.8]);
%pcolor(dist_grid,-depth_grid,SIG_0_wint); shading flat;
contourf(dist_grid,pres_grid,GRID.dens(:,:,season),levels); shading flat;
% plot(dist_grid,bathy_pres,'k','linewidth',2); 

set(gca,'ydir','reverse');
%ylim([0 4800])
 ylim([0 1000])
ylabel('Pres (db)','fontsize',14);
xlabel('Distance along contour (km)','fontsize',14);
title(GRID.season_description(season,:),'fontsize',14)
set(gca,'fontsize',14);
colours = (lbmap(20,'RedBlue'));
colormap(flipud(colours));
cb = colorbar;
ylabel(cb,'Density Sig\theta (kg m^-^3 - 1000)','fontsize',14);
end


%% print figure
width  = 1500;  % frame width
height = 2000;  % frame height
pngname = ('plots/p5_dens_season_EN4added');

% set background color (outside axes)
set(gcf,'color',[1 1 1]);

% don't change background color when printing
set(gcf,'inverthardcopy','off');

% set size of frame to be written
resolution=150;
set(gcf,'paperunits','inches');
set(gcf,'paperposition',[0 0 width height]./resolution);

% write .png file
% the 'zbuffer' method is likely to look similar to the plot window
print('-dpng', ['-r' num2str(resolution)], '-opengl', pngname);







% %% Save variables
save D:\Work_computer_sync\OSNAP_postdoc\PAPERS_NEW\N_Atlantic_boundary\matlab\V3_050321\intermediate_saves\pressure_gridded_1000_EN4inserted dist_grid pres_grid lon_grid lat_grid transect_mask_grid GRID bathy_depth bathy_pres
% save('D:\Work_computer_sync\OSNAP_postdoc\PAPERS_NEW\N_Atlantic_boundary\matlab\V3_050321\intermediate_saves\contour_data_1000_EN4inserted',...
%    'cont_distance','cont_dist_cumul','cont_lon','cont_lat','zonal_transect_mask',...
%    'dist_grid','lon_grid','lat_grid','dist_grid_accurate','transect_mask_grid');
