% p3e_gridding_P_grid_4seasons
% SJ 03/21
%
% Gridding of the 1000 m-adjacent component of the dataset.

clear; close all;

addpath(genpath('D:\Work_computer_sync\MATLAB_functions'));


minlon = -65;
maxlon = 5;
minlat = 46;
maxlat = 67;

%% Gridding parameters 
pres_grid = 20:20:5000;
% density_grid_size = 0.015; % kg/m3
preferred_profile_count = 75;
min_search_radius = 150;
max_search_radius = 300;
% dist_grid = 150:150:11900;

%% Load QC'd raw profiles

load D:\Work_computer_sync\OSNAP_postdoc\PAPERS_NEW\N_Atlantic_boundary\matlab\V3_050321\intermediate_saves\contour_data_1000
load D:\Work_computer_sync\OSNAP_postdoc\PAPERS_NEW\N_Atlantic_boundary\matlab\V3_050321\intermediate_saves\raw_profs_1000m_P2_QC

ocean.month = nan*ones(1,length(ocean.juld));
for aa = 1:length(ocean.juld)
    ocean.month(1,aa) = str2num(datestr(ocean.juld(aa),'mm'));
end

%ind = min(find(zonal_transect_mask));
%start_of_zonal = cont_dist_cumul(ind);

% preallocate main gridded arrays
GRID.season_description = ['Winter JFM'; 'Spring AMJ'; 'Summer JAS'; 'Autumn OND'];
GRID.SA = nan*ones(length(pres_grid),length(dist_grid),4);
GRID.CT = GRID.SA;
GRID.dens = GRID.SA;
prof_derived_bathy = nan*ones(1,length(dist_grid),4);

season_months = [1 2 3; 4 5 6; 7 8 9; 10 11 12];

for season = 1:4 % for each season
    
    for bb = 1:length(dist_grid) % for each along-slope profile
        
        %% Gradually expand our search radius until we get the requisite number
        % of observations
        search_radius = min_search_radius;
        num_prof = 0;
        
        while num_prof < preferred_profile_count & search_radius < max_search_radius
            profind = find((ocean.cont_distance > (dist_grid(bb)-search_radius)& ocean.cont_distance < (dist_grid(bb)+search_radius)) ...
                & (ocean.month == season_months(season,1) | ocean.month == season_months(season,2) | ocean.month == season_months(season,3)));
            
            % S = sum(X,DIM)
            goodprof = nanmean(ocean.sigma0(:,profind),1);
            num_prof = length(find(~isnan(goodprof)));
            
            
            
            if num_prof >= preferred_profile_count
                break
            else
                search_radius = search_radius + 25;
            end
        end
        
        % Add a couple of diagnostics to the output
        % GRID.number_of_profiles(1,bb,season) = num_prof;
        % GRID.search_radius(1,bb,season) = search_radius;
        
        %% average bathy depth if covering the zonal part
%         if dist_grid(bb) < start_of_zonal % if in the 1000 m contour region
%             prof_derived_bathy(1,bb) = 1000;
%         else % if traversing the zonal transect
%             prof_derived_bathy(1,bb) = -nanmean(ocean.profdepth(profind));
%         end
        
        prof_derived_bathy(1,bb) = 1000;
        
        
        %% This bit actually does the averaging.
        
        if ~isempty(profind)
            sig0_temp = ocean.sigma0(:,profind);
            CT_temp = ocean.CT(:,profind);
            SA_temp = ocean.SA(:,profind);
            pres_temp = ocean.pres(:,profind);
            
            %% QC figure 100: plot of all raw profiles, with average on top
            figure(100);
            clf
            hold on;
            for ff = 1:length(profind)
                A = plot(sig0_temp(:,ff),pres_temp(:,ff),'g','linewidth',1);
                plot(sig0_temp(:,ff),pres_temp(:,ff),'g.','markersize',10);
            end % end of plotting loop
            set(gca,'ydir','reverse');
            ylim([0 prof_derived_bathy(1,bb)]);
            xlim([26 27.9]);
            xlabel('Density');
            ylabel('pressure');
            title(['Gridding - raw vs. gridded for cell ' num2str(bb) ', Distance: ' num2str(dist_grid(bb)) ' km']);
            
            % preallocate empty arrays to take temporary profile group
            dens_pgridded_UNAVERAGED = nan*ones(length(pres_grid),length(profind));
            CT_pgridded_UNAVERAGED = nan*ones(length(pres_grid),length(profind));
            SA_pgridded_UNAVERAGED = nan*ones(length(pres_grid),length(profind));
            pres_pgridded_UNAVERAGED = nan*ones(length(pres_grid),length(profind));
            
            
            %% interpolate all the raw profiles onto the pressure grid
            for dd = 1:length(profind) % for each profile in this grab...
                ind = find(~isnan(pres_temp(:,dd)));
                if ~isempty(ind)
                    % Vq = interp1(X,V,Xq)
                    CT_pgridded_UNAVERAGED(:,dd) = interp1(pres_temp(ind,dd),CT_temp(ind,dd),pres_grid);
                    SA_pgridded_UNAVERAGED(:,dd) = interp1(pres_temp(ind,dd),SA_temp(ind,dd),pres_grid);
                end
            end % end of 'for each profile'
            
            % Slightly redundant 'mean' step, but didn't want any indexing
            % mistakes!
            CT_pgridded_MEAN = nan*ones(length(pres_grid),1);
            SA_pgridded_MEAN = nan*ones(length(pres_grid),1);
            sig0_pgridded_MEAN = nan*ones(length(pres_grid),1);
            
            for ee = 1:length(pres_grid)
                %             CT_pgridded_MEAN(ee,1) = median(CT_pgridded(ee,:),'omitnan');
                %             SA_pgridded_MEAN(ee,1) = median(SA_pgridded(ee,:),'omitnan');
                CT_pgridded_MEAN(ee,1) = nanmean(CT_pgridded_UNAVERAGED(ee,:));
                SA_pgridded_MEAN(ee,1) = nanmean(SA_pgridded_UNAVERAGED(ee,:));
            end
            % sigma0 = gsw_sigma0(SA,CT)
            sig0_pgridded_MEAN = gsw_sigma0(SA_pgridded_MEAN,CT_pgridded_MEAN);
            
            %% populate main arrays
            GRID.SA(:,bb,season) = SA_pgridded_MEAN;
            GRID.CT(:,bb,season) = CT_pgridded_MEAN;
            GRID.dens(:,bb,season) = sig0_pgridded_MEAN;
            
            
            %% Add to gridding plot...
            
            % this is interpolated, sorted composite
            E = plot(GRID.dens(:,bb,season),pres_grid,'k','linewidth',1.5);
            plot(GRID.dens(:,bb,season),pres_grid,'k.','markersize',12);
            
            legend([A,E],'Raw profile data','Pressure gridded composite','Location','SouthWest');
            grid on
            
            % print figure
            width  = 1500;  % frame width
            height = 1500;  % frame height
            pngname = (['plots/gridding_intermediate/season_' num2str(season) '/Gridding_sanity_check_for_cell_' num2str(bb) '.png']);
            
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
            
        end % end of 'if there are any profiles in this along slope bin...'
        
        bb
    end % 'for each along-slope profile...'
    
end % End season loop

%% Trim profiles based on pressure of bed 
% prof_derived_bathy = nan*ones(1,length(dist_grid),4);
bathy_depth = nanmean(prof_derived_bathy,3);

bathy_pres = gsw_p_from_z(-bathy_depth,lat_grid);
for ee = 1:length(bathy_depth) % for each profile
    ind = find(pres_grid > bathy_pres(ee));
    GRID.SA(ind,ee,:) = nan;
    GRID.CT(ind,ee,:) = nan;
    GRID.dens(ind,ee,:) = nan;
end







%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Plot 1: Scatter of data locations

figure(1); hold on;

bath = load('GEBCO_world_1D.mat');
m_proj('mercator','longitudes',[minlon maxlon],'latitudes',[minlat maxlat]);

bath.bathy = bath.bathy';

% bathymetry
ind = find(bath.longitude >= minlon & bath.longitude < maxlon); bath.longitude = bath.longitude(ind); bath.bathy = bath.bathy(:,ind);
ind = find(bath.latitude >= minlat & bath.latitude < maxlat); bath.latitude = bath.latitude(ind); bath.bathy = bath.bathy(ind,:);

%bath.bathy(bath.bathy > 0) = 0;
bathinv = bath.bathy.*-1;

m_contour(bath.longitude,bath.latitude,bathinv,[1000 2000 3000 4000 5000],'color',[0.6 0.6 0.6],'linewidth',1.5);
m_contour(bath.longitude,bath.latitude,bathinv,[0 0],'color','k','linewidth',2);

%caxis([6200 7500])
m_plot(cont_lon,cont_lat,'.r','linewidth',4);
m_plot(lon_grid,lat_grid,'+k','linewidth',2);

% for aa = 1:4:length(dist_grid)
%    m_text(lon_grid(aa)+0.1,lat_grid(aa)+0.1,[num2str(dist_grid(aa)) ' km'],'fontsize',12); 
% end

% plot lon and lat of profiles
m_plot(ocean.lon,ocean.lat,'.','color',[0.4 0.4 0.4]);

% Distance labels
for aa = 5:5:length(dist_grid)
   A = m_plot(lon_grid(aa),lat_grid(aa),'+k','markersize',20,'linewidth',3);
   uistack(A,'top');
   m_text(lon_grid(aa)+0.8,lat_grid(aa)+1.3,[num2str(dist_grid(aa)) ' km'],'fontsize',12,'color','k','BackgroundColor','w'); 
end


% tidying uup
m_grid('box','fancy','fontsize',20,'color','k','tickdir','in');


%% print figure
width  = 2000;  % frame width
height = 1000;  % frame height
pngname = ('plots/p2_map');

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






%% Fig 3: T,S, dens contour plots

figure(3);
for season = 1:4
subplot(4,1,season);
hold on;
% dist_grid,CT_wint
levels = [0:0.5:14];
caxis([2 14]);
contourf(dist_grid,pres_grid,GRID.CT(:,:,season),levels); shading flat;
plot(dist_grid,bathy_pres,'k','linewidth',2); 
set(gca,'ydir','reverse');
%ylim([0 4700])
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
pngname = ('plots/p3_CT');

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
levels = [33:0.05:35.7];
caxis([34.5 35.7]);
%pcolor(dist_grid,-depth_grid,SA_wint); shading flat;
contourf(dist_grid,pres_grid,GRID.SA(:,:,season),levels); shading flat;
plot(dist_grid,bathy_pres,'k','linewidth',2); 
set(gca,'ydir','reverse');
%ylim([0 4700])
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
pngname = ('plots/p4_SA');

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
plot(dist_grid,bathy_pres,'k','linewidth',2); 
set(gca,'ydir','reverse');

%ylim([0 4700])
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
pngname = ('plots/p5_dens');

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
% NEW SAVE VARIABLES BELOW
save D:\Work_computer_sync\OSNAP_postdoc\PAPERS_NEW\N_Atlantic_boundary\matlab\V3_050321\intermediate_saves\pressure_gridded_1000 dist_grid pres_grid lon_grid lat_grid GRID bathy_depth bathy_pres

