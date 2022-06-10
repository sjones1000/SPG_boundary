% p5_calc_velocities_and_transports
% S Jones 140421 based on SC's code.
% Note I don't think the transports calc'd here are actually used.  We
% recalculate them as needed when working on fluxes, overturning, etc.

clearvars ; close all ;
addpath(genpath('D:\Work_computer_sync\MATLAB_functions'));
addpath subfunctions

%% Load

% Load boundary product
load('D:\Work_computer_sync\OSNAP_postdoc\PAPERS_NEW\N_Atlantic_boundary\matlab\V3_050321\intermediate_saves/pressure_gridded_1000_EN4inserted.mat')
% Load contour data
contour = load('D:\Work_computer_sync\OSNAP_postdoc\PAPERS_NEW\N_Atlantic_boundary\matlab\V3_050321\intermediate_saves/contour_data_1000_EN4inserted.mat');
% Load satellite altimetry (already formatted for boundary product)
load('D:\Work_computer_sync\OSNAP_postdoc\PAPERS_NEW\N_Atlantic_boundary\matlab\V3_050321\PLOT6_altimetry_along_grid/boundary_ssh.mat')
% load bathymetry

%% Prep datasets

g = 9.82 ; % acceleration due to gravity

% Smooth SSH using approximately the same radius as hydrography.  Note the
% smooth distance varies according to the grid spacing, but this seems
% appropriate to boundary / EN4 data origins.
SSH.adt_JFM = smooth(SSH.adt_JFM,5)';
SSH.adt_AMJ = smooth(SSH.adt_AMJ,5)';
SSH.adt_JAS = smooth(SSH.adt_JAS,5)';
SSH.adt_OND = smooth(SSH.adt_OND,5)';

% Distance grid for plotting velocities
for aa = 1:length(dist_grid)-1
    %dx2(aa) = dist_grid(aa+1) - dist_grid(aa);
    dist_grid_vel(aa) = dist_grid(aa) + ((dist_grid(aa+1) - dist_grid(aa)) / 2);
end

% Distance between stations [NOTE MUST AGREE WITH DISTANCE METHOD USED FOR GEOSTROPHIC TRANSPORT CALC]
dx=diff(dist_grid).*1000 ; % dx in m <- distance following contour
% dx = gsw_distance(lon_grid,lat_grid); <- straight line distance

% Add zero db level to pres_grid for transport calcs, copy 20 db values to
% surface. SA, CT have dimensions (x,pres,season)
pres_grid = [0 pres_grid];
GRID.SA = cat(1,GRID.SA(1,:,:),GRID.SA(:,:,:));
GRID.CT = cat(1,GRID.CT(1,:,:),GRID.CT(:,:,:));
% sigma0 = gsw_sigma0(SA,CT)
GRID.dens = gsw_sigma0(GRID.SA,GRID.CT);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1.  Compute main product... surface-referenced geostrophic velocities, summing
% dynamic height and satellite-adt then calculating geo vels from resultant gradients
% Note 3D variables SA, CT have dimensions (pres,x,season) where season is 1-4

% TOTAL_above_1000m geo vel (ref surface): JFM
dy_ht_JFM = gsw_geo_strf_dyn_height(GRID.SA(:,:,1),GRID.CT(:,:,1),pres_grid',0) ; % [m^2/s^2]
dy_ht_surf_JFM = dy_ht_JFM + (g * SSH.adt_JFM);
% [TOTAL_above_1000m.gvel_JFM, mid_lat, mid_long] = gsw_geostrophic_velocity(dy_ht_surf_JFM,lon_grid,lat_grid); % <- straight line distance between stations
[TOTAL_above_1000m.gvel_JFM, mid_lat, mid_long] = gsw_geostrophic_velocity_DISTANCE_SUPPLIED(dy_ht_surf_JFM,lon_grid,lat_grid,dx); % <- distance around contour

% TOTAL_above_1000m geo vel (ref surface): AMJ
dy_ht_AMJ = gsw_geo_strf_dyn_height(GRID.SA(:,:,2),GRID.CT(:,:,2),pres_grid',0) ; % [m^2/s^2]
dy_ht_surf_AMJ = dy_ht_AMJ + (g * SSH.adt_AMJ);
%[TOTAL_above_1000m.gvel_AMJ, mid_lat, mid_long] = gsw_geostrophic_velocity(dy_ht_surf_AMJ,lon_grid,lat_grid); % <- straight line distance between stations
[TOTAL_above_1000m.gvel_AMJ, mid_lat, mid_long] = gsw_geostrophic_velocity_DISTANCE_SUPPLIED(dy_ht_surf_AMJ,lon_grid,lat_grid,dx); % <- distance around contour

% TOTAL_above_1000m geo vel (ref surface): JAS
dy_ht_JAS = gsw_geo_strf_dyn_height(GRID.SA(:,:,3),GRID.CT(:,:,3),pres_grid',0) ; % [m^2/s^2]
dy_ht_surf_JAS = dy_ht_JAS + (g * SSH.adt_JAS);
%[TOTAL_above_1000m.gvel_JAS, mid_lat, mid_long] = gsw_geostrophic_velocity(dy_ht_surf_JAS,lon_grid,lat_grid); % <- straight line distance between stations
[TOTAL_above_1000m.gvel_JAS, mid_lat, mid_long] = gsw_geostrophic_velocity_DISTANCE_SUPPLIED(dy_ht_surf_JAS,lon_grid,lat_grid,dx); % <- distance around contour

% TOTAL_above_1000m geo vel (ref surface): OND
dy_ht_OND = gsw_geo_strf_dyn_height(GRID.SA(:,:,4),GRID.CT(:,:,4),pres_grid',0) ; % [m^2/s^2]
dy_ht_surf_OND = dy_ht_OND + (g * SSH.adt_OND);
%[TOTAL_above_1000m.gvel_OND, mid_lat, mid_long] = gsw_geostrophic_velocity(dy_ht_surf_OND,lon_grid,lat_grid); % <- straight line distance between stations
[TOTAL_above_1000m.gvel_OND, mid_lat, mid_long] = gsw_geostrophic_velocity_DISTANCE_SUPPLIED(dy_ht_surf_OND,lon_grid,lat_grid,dx); % <- distance around contour

% Create a bathy mask while we have the full cross-section of velocities (the deep portion is removed later)
bathy_mask = 0* TOTAL_above_1000m.gvel_JFM;
bathy_mask = bathy_mask+1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1.5.  Compute correct cross-sectional areas before transport calcs %%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx_metres = dx;
%[lat2d,p2d] = meshgrid(mid_lat,pres_grid');
%z_grid = gsw_z_from_p(p2d,lat2d); % if interested in the difference caused
%by working with depth instead of pressure...
% create dz (basic 'canvas', with 0-10 m half the dz of the other cells)
for aa = 1:length(pres_grid)
    % top clause
    if aa == 1
        dz_metres(aa,1) = 10;
    else
        dz_metres(aa,1) = 20;
    end
end
dz_metres = repmat(dz_metres,1,length(mid_lat));

% Next trim using 'velocity bathymetry'
dz_metres = dz_metres.*bathy_mask;

% The bottom cell is only half dz (i.e. 990-1000 dbar).
for aa = 1:length(dz_metres(1,:)) % for each column
   ind = max(find(~isnan(dz_metres(:,aa)))); 
   dz_metres(ind,aa) = 10;
end

% Compute area
area_metres2 = dx_metres .* dz_metres;

% Now deal with edge area under-representation due to geo vel process
edgemask = double(~isnan(area_metres2));
edgemask(edgemask==0) = nan;
for aa = 1:length(edgemask(:,1)) % for each row
    for bb = 2:(length(edgemask(1,:)))-1 % for each column [note missing start and end column]      
        if ~isnan(area_metres2(aa,bb)) % if the cell is in the water column
            if isnan(edgemask(aa,bb-1)) | isnan(edgemask(aa,bb+1)) % if next to a blank cell
                edgemask(aa,bb) = 1.5;
            end
        end
    end
end

figure(1);
subplot(2,1,1);
pcolor(edgemask); shading flat;
title('Edge multiplier check!');
colorbar
set(gca,'ydir','reverse');

% This is the important bit!
area_metres2 = area_metres2.*edgemask;

subplot(2,1,2);
pcolor(area_metres2); shading flat;
title('Area metres^2');
colorbar
set(gca,'ydir','reverse');

% Also want an area 'map' for above and below 1000 m. [Note this is specific to this number of grid cells and will break if they change] 
area_metres2_above1000m = area_metres2;
area_metres2_below1000m = area_metres2;
% Trim off bottom, halve 'exposed' cells...
area_metres2_above1000m(52:end,:) = nan;
area_metres2_above1000m(51,64:100) = area_metres2_above1000m(51,64:100)./2;
dz_metres_above1000m = dz_metres; 
dz_metres_above1000m(52:end,:) = nan;
dz_metres_above1000m(51,64:100) = dz_metres_above1000m(51,64:100)./2;
% Trim off top, halve 'exposed' cells...
area_metres2_below1000m(1:51,:) = nan;
area_metres2_below1000m(52,64:100) = area_metres2_below1000m(52,64:100)./2;
dz_metres_below1000m = dz_metres; 
dz_metres_below1000m(1:51,:) = nan;
dz_metres_below1000m(52,64:100) = dz_metres_below1000m(52,64:100)./2;

figure(2)

subplot(3,1,1);
pcolor(area_metres2); shading flat;
title('Area metres^2 - TOTAL');
colorbar
set(gca,'ydir','reverse');

subplot(3,1,2);
pcolor(area_metres2_above1000m); shading flat;
title('Above 1000m');
colorbar
set(gca,'ydir','reverse');

subplot(3,1,3);
pcolor(area_metres2_below1000m); shading flat;
title('Below 1000m');
colorbar
set(gca,'ydir','reverse');





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1.75 Now compute transports using the cell areas we calculated above

% Compute approx depth average velocity to 1000 m (not needed for tranport calcs, just plot)
DA_total_JFM = nanmean(TOTAL_above_1000m.gvel_JFM(1:51,:),1);

% Compute volume transports using more accurate 'per cell' method - first
% transport per cell
TOTAL_above_1000m.tran_JFM_2D = (TOTAL_above_1000m.gvel_JFM .* area_metres2_above1000m) / 1E6 ; % TRANSPORT IN Sv
TOTAL_above_1000m.tran_AMJ_2D = (TOTAL_above_1000m.gvel_AMJ .* area_metres2_above1000m) / 1E6 ;
TOTAL_above_1000m.tran_JAS_2D = (TOTAL_above_1000m.gvel_JAS .* area_metres2_above1000m) / 1E6 ;
TOTAL_above_1000m.tran_OND_2D = (TOTAL_above_1000m.gvel_OND .* area_metres2_above1000m) / 1E6 ;

% now depth integral
TOTAL_above_1000m.tran_JFM = sum(TOTAL_above_1000m.tran_JFM_2D,1,'omitnan'); % TRANSPORT IN Sv
TOTAL_above_1000m.tran_AMJ = sum(TOTAL_above_1000m.tran_AMJ_2D,1,'omitnan');
TOTAL_above_1000m.tran_JAS = sum(TOTAL_above_1000m.tran_JAS_2D,1,'omitnan');
TOTAL_above_1000m.tran_OND = sum(TOTAL_above_1000m.tran_OND_2D,1,'omitnan');

% Compute cumulative transports - NEW method
TOTAL_above_1000m.tran_cum_JFM = cumsum(TOTAL_above_1000m.tran_JFM) ;
TOTAL_above_1000m.tran_cum_AMJ = cumsum(TOTAL_above_1000m.tran_AMJ) ;
TOTAL_above_1000m.tran_cum_JAS = cumsum(TOTAL_above_1000m.tran_JAS) ;
TOTAL_above_1000m.tran_cum_OND = cumsum(TOTAL_above_1000m.tran_OND) ;


% Quick figure to plot the above data products
figure(3)
subplot(2,1,1)
hold on
title('Transport above 1000 m (Sv, out of volume negative)')
plot(dist_grid_vel,TOTAL_above_1000m.tran_JFM,'b');
plot(dist_grid_vel,TOTAL_above_1000m.tran_AMJ,'g');
plot(dist_grid_vel,TOTAL_above_1000m.tran_JAS,'r');
plot(dist_grid_vel,TOTAL_above_1000m.tran_OND,'y');
grid on

subplot(2,1,2)
title('Cumulative transport above 1000 m (Sv)')
hold on; 
plot(dist_grid_vel,TOTAL_above_1000m.tran_cum_JFM,'b');
plot(dist_grid_vel,TOTAL_above_1000m.tran_cum_AMJ,'g');
plot(dist_grid_vel,TOTAL_above_1000m.tran_cum_JAS,'r');
plot(dist_grid_vel,TOTAL_above_1000m.tran_cum_OND,'y');
legend('JFM','AMJ','JAS','OND','Location','best');
grid on



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Compute surface vels only from altimetry %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SSH first difference (m) 
dadt_JFM = diff(SSH.adt_JFM) ; 
dadt_AMJ = diff(SSH.adt_AMJ) ; 
dadt_JAS = diff(SSH.adt_JAS) ; 
dadt_OND = diff(SSH.adt_OND) ; 

% compute f at the velocity midpoints
for i = 1 : length(SSH.lon) - 1
    fmn(i) = sw_f(mean([SSH.lat(i) SSH.lat(i+1)])) ;
end

gfmn = g ./ fmn ;

dadtdx_JFM = dadt_JFM ./ dx ;
dadtdx_AMJ = dadt_AMJ ./ dx ;
dadtdx_JAS = dadt_JAS ./ dx ;
dadtdx_OND = dadt_OND ./ dx ;

BT.adt_vel_JFM = ( gfmn .* dadtdx_JFM ) ;
BT.adt_vel_AMJ = ( gfmn .* dadtdx_AMJ ) ;
BT.adt_vel_JAS = ( gfmn .* dadtdx_JAS ) ;
BT.adt_vel_OND = ( gfmn .* dadtdx_OND ) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Compute Baroclinic vels only, referenced to 1000 m %%%%%%%%%%%%%%%%%%

% Also calc baroclinic only.
dy_ht_bc_JFM = gsw_geo_strf_dyn_height(GRID.SA(:,:,1),GRID.CT(:,:,1),pres_grid',1000) ; % [m^2/s^2]
% [BC.gvel_bc_JFM, mid_lat, mid_long] =
% gsw_geostrophic_velocity(dy_ht_bc_JFM,lon_grid,lat_grid);  % <- straight line distance between stations
[BC.gvel_bc_JFM, mid_lat, mid_long] = gsw_geostrophic_velocity_DISTANCE_SUPPLIED(dy_ht_bc_JFM,lon_grid,lat_grid,dx);

% Also calc baroclinic only.
dy_ht_bc_AMJ = gsw_geo_strf_dyn_height(GRID.SA(:,:,2),GRID.CT(:,:,2),pres_grid',1000) ; % [m^2/s^2]
%[BC.gvel_bc_AMJ, mid_lat, mid_long] = gsw_geostrophic_velocity(dy_ht_bc_AMJ,lon_grid,lat_grid); % <- straight line distance between stations
[BC.gvel_bc_AMJ, mid_lat, mid_long] = gsw_geostrophic_velocity_DISTANCE_SUPPLIED(dy_ht_bc_AMJ,lon_grid,lat_grid,dx);

% Also calc baroclinic only.
dy_ht_bc_JAS = gsw_geo_strf_dyn_height(GRID.SA(:,:,3),GRID.CT(:,:,3),pres_grid',1000) ; % [m^2/s^2]
%[BC.gvel_bc_JAS, mid_lat, mid_long] = gsw_geostrophic_velocity(dy_ht_bc_JAS,lon_grid,lat_grid); % <- straight line distance between stations
[BC.gvel_bc_JAS, mid_lat, mid_long] = gsw_geostrophic_velocity_DISTANCE_SUPPLIED(dy_ht_bc_JAS,lon_grid,lat_grid,dx);

% Also calc baroclinic only.
dy_ht_bc_OND = gsw_geo_strf_dyn_height(GRID.SA(:,:,4),GRID.CT(:,:,4),pres_grid',1000) ; % [m^2/s^2]
%[BC.gvel_bc_OND, mid_lat, mid_long] = gsw_geostrophic_velocity(dy_ht_bc_OND,lon_grid,lat_grid); % <- straight line distance between stations
[BC.gvel_bc_OND, mid_lat, mid_long] = gsw_geostrophic_velocity_DISTANCE_SUPPLIED(dy_ht_bc_OND,lon_grid,lat_grid,dx);

% Compute depth average velocity to 1000 m [not used for tranport calc, just plot]
DA_baroclinic_JFM = nanmean(BC.gvel_bc_JFM(1:51,:),1);

% Plot of baroclinic, barotropic, total
figure(4);
title('Comparison between total, baroclinic and barotropic velocities [JFM] (ms-1)');
hold on;
plot(dist_grid_vel,DA_total_JFM,'r','linewidth',2);
plot(dist_grid_vel,DA_baroclinic_JFM,'b');
plot(dist_grid_vel,BT.adt_vel_JFM,'k');
legend('Total','Baroclinic','Barotropic','Location','best');
grid on;

% contour of baroclinic vs barotropic vels >  1000 db
figure(5);
subplot(2,1,1)
hold on
title('Total vels - JFM (m/s)')
levels = [-0.3:0.01:0.3];
contourf(dist_grid_vel,-pres_grid,TOTAL_above_1000m.gvel_JFM,levels);
ylim([-1000 0]);
caxis([-0.05 0.05]);
colorbar;

subplot(2,1,2)
hold on
title('Baroclinic vels - JFM (m/s)')
contourf(dist_grid_vel,-pres_grid,BC.gvel_bc_JFM,levels);
ylim([-1000 0]);
caxis([-0.05 0.05]);
colorbar;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Compute velocities and transports for sub - 1000 m
% Start with total velocity product as calc'd above
% Set all velocities above 1000 m to nan. 

DEEP.gvel_JFM = TOTAL_above_1000m.gvel_JFM;
DEEP.gvel_AMJ = TOTAL_above_1000m.gvel_AMJ;
DEEP.gvel_JAS = TOTAL_above_1000m.gvel_JAS;
DEEP.gvel_OND = TOTAL_above_1000m.gvel_OND;

DEEP.gvel_JFM(1:51,:) = nan;
DEEP.gvel_AMJ(1:51,:) = nan;
DEEP.gvel_JAS(1:51,:) = nan;
DEEP.gvel_OND(1:51,:) = nan;

% Compute volume transports using more accurate 'per cell' method - first
% transport per cell
DEEP.tran_JFM_2D = (DEEP.gvel_JFM .* area_metres2_below1000m) / 1E6 ; % TRANSPORT IN Sv
DEEP.tran_AMJ_2D = (DEEP.gvel_AMJ .* area_metres2_below1000m) / 1E6 ;
DEEP.tran_JAS_2D = (DEEP.gvel_JAS .* area_metres2_below1000m) / 1E6 ;
DEEP.tran_OND_2D = (DEEP.gvel_OND .* area_metres2_below1000m) / 1E6 ;

% now depth integral
DEEP.tran_JFM = sum(DEEP.tran_JFM_2D,1,'omitnan'); % TRANSPORT IN Sv
DEEP.tran_AMJ = sum(DEEP.tran_AMJ_2D,1,'omitnan');
DEEP.tran_JAS = sum(DEEP.tran_JAS_2D,1,'omitnan');
DEEP.tran_OND = sum(DEEP.tran_OND_2D,1,'omitnan');

% Compute cumulative transports
DEEP.tran_cum_JFM = cumsum(DEEP.tran_JFM,'omitnan') ;
DEEP.tran_cum_AMJ = cumsum(DEEP.tran_AMJ,'omitnan') ;
DEEP.tran_cum_JAS = cumsum(DEEP.tran_JAS,'omitnan') ;
DEEP.tran_cum_OND = cumsum(DEEP.tran_OND,'omitnan') ;

% Couple of diagnostic plots for deep component
figure(6);
subplot(3,1,1)
hold on
title('Deep vels - JFM (m/s)')
levels = [-0.3:0.01:0.3];
contourf(dist_grid_vel,-pres_grid,DEEP.gvel_JFM,levels);
xlim([9500 max(dist_grid_vel)]);
ylim([-5000 0]);
caxis([-0.05 0.05]);
colorbar;

subplot(3,1,2)
hold on
title('Transport below 1000 m (Sv, out of volume negative)')
plot(dist_grid_vel,DEEP.tran_JFM,'b');
plot(dist_grid_vel,DEEP.tran_AMJ,'g');
plot(dist_grid_vel,DEEP.tran_JAS,'r');
plot(dist_grid_vel,DEEP.tran_OND,'y');
xlim([9500 max(dist_grid_vel)]);
grid on

subplot(3,1,3)
hold on
title('Cumulative transport below 1000 m (Sv, out of volume negative)')
plot(dist_grid_vel,DEEP.tran_cum_JFM,'b');
plot(dist_grid_vel,DEEP.tran_cum_AMJ,'g');
plot(dist_grid_vel,DEEP.tran_cum_JAS,'r');
plot(dist_grid_vel,DEEP.tran_cum_OND,'y');
xlim([9500 max(dist_grid_vel)]);
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. Redo mid-lons, mid-lats to be along the contour rather than halfway along the straight line between two grid points
mid_long2 = mid_long;
mid_lat2 = mid_lat;
boundind = find(transect_mask_grid == 0);
% ONLY OPERATING ON BOUNDARY BIT, AS 'NEAREST NEIGHBOR DOESN'T WORK ON EN4
% CONTOUR DATA.  AND REMOVING LAST 'WRAP-AROUND' CELL.
for bb = 1:length(boundind)-2 % for each along-slope cell
    % locate lon and lat of each grid point by minimisation
    dist_delete = abs(dist_grid_vel(bb) - contour.cont_dist_cumul); distminind = find(dist_delete == min(dist_delete)); distind = distminind(1);
    mid_long2(1,bb) = contour.cont_lon(distind);
    mid_lat2(1,bb) = contour.cont_lat(distind); 
end

figure(7);
hold on;
plot(contour.cont_lon,contour.cont_lat,'color',[0.5 0.5 0.5]);
plot(lon_grid,lat_grid,'.b');
plot(mid_long,mid_lat,'.r');
plot(mid_long2,mid_lat2,'xg')
legend('contour data','lat grid','mid lat OLD','mid lat modified');

% Happy with Fig 5 so keep the new mid_lons and mid_lats
mid_long = mid_long2;
mid_lat = mid_lat2;

%% 6. Interpolate T and S onto velocity 'mid-point' coordinates. Recalculate density.

% prep dimensions
[dist_in,pres_in] = meshgrid(dist_grid,pres_grid);
[dist_out,pres_out] = meshgrid(dist_grid_vel,pres_grid);

% Vq = interp2(X,Y,V,Xq,Yq)
WATER_PROPERTIES.CT_JFM = interp2(dist_in,pres_in,squeeze(GRID.CT(:,:,1)),dist_out,pres_out);
WATER_PROPERTIES.CT_AMJ = interp2(dist_in,pres_in,squeeze(GRID.CT(:,:,2)),dist_out,pres_out);
WATER_PROPERTIES.CT_JAS = interp2(dist_in,pres_in,squeeze(GRID.CT(:,:,3)),dist_out,pres_out);
WATER_PROPERTIES.CT_OND = interp2(dist_in,pres_in,squeeze(GRID.CT(:,:,4)),dist_out,pres_out);

WATER_PROPERTIES.SA_JFM = interp2(dist_in,pres_in,squeeze(GRID.SA(:,:,1)),dist_out,pres_out);
WATER_PROPERTIES.SA_AMJ = interp2(dist_in,pres_in,squeeze(GRID.SA(:,:,2)),dist_out,pres_out);
WATER_PROPERTIES.SA_JAS = interp2(dist_in,pres_in,squeeze(GRID.SA(:,:,3)),dist_out,pres_out);
WATER_PROPERTIES.SA_OND = interp2(dist_in,pres_in,squeeze(GRID.SA(:,:,4)),dist_out,pres_out);

% sigma0 = gsw_sigma0(SA,CT)
WATER_PROPERTIES.sig0_JFM = gsw_sigma0(WATER_PROPERTIES.SA_JFM,WATER_PROPERTIES.CT_JFM);
WATER_PROPERTIES.sig0_AMJ = gsw_sigma0(WATER_PROPERTIES.SA_AMJ,WATER_PROPERTIES.CT_AMJ);
WATER_PROPERTIES.sig0_JAS = gsw_sigma0(WATER_PROPERTIES.SA_JAS,WATER_PROPERTIES.CT_JAS);
WATER_PROPERTIES.sig0_OND = gsw_sigma0(WATER_PROPERTIES.SA_OND,WATER_PROPERTIES.CT_OND);

save intermediate_saves/SPG_geovels dist_grid_vel mid_long mid_lat pres_grid dx_metres ...
    dz_metres dz_metres_above1000m dz_metres_below1000m area_metres2 ...
    area_metres2_above1000m area_metres2_below1000m TOTAL_above_1000m BT BC DEEP WATER_PROPERTIES

