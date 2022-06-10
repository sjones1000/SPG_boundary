% p2_QC_WOD
%
% Perform quality control on WOD profiles.
% SJ 03/21
% QC subsetted profile data
%%%%% Takes about 20 mins to run on full dataset %%%%%%%%%

clear; close all;

addpath(genpath('D:\Work_computer_sync\MATLAB_functions'));
%addpath(genpath('D:\Work_computer_sync\OSNAP_postdoc\PAPERS_NEW\N_Atlantic_boundary\matlab\V2_250220\intermediate_saves'));

min_profile_length = 50;

% load bathy
%load('GEBCO_world_1D.mat')

minlon = -65;
maxlon = 5;
minlat = 46;
maxlat = 67;

% Load raw profiles
load D:\Work_computer_sync\OSNAP_postdoc\PAPERS_NEW\N_Atlantic_boundary\matlab\V3_050321\intermediate_saves\raw_profs_1000m
load D:\Work_computer_sync\OSNAP_postdoc\PAPERS_NEW\N_Atlantic_boundary\matlab\V3_050321\intermediate_saves\contour_data_1000


% First things first; find and destroy data from Faroe Shetland Channel!
ind = find(ocean.lon > -8 & ocean.lon < -5.5 & ocean.lat > 59.95 & ocean.lat < 60.6);
ocean.temp(:,ind) = nan;
ocean.sal(:,ind) = nan;
ocean.pres(:,ind) = nan;

% And region separating Flemish Cap from shelf
ind = find(ocean.lon > -47.5 & ocean.lon < -45.2 & ocean.lat > 46.5 & ocean.lat < 48.1);
ocean.temp(:,ind) = nan;
ocean.sal(:,ind) = nan;
ocean.pres(:,ind) = nan;

%% Sift 0.5 ... We want full depth profiles, with some represenatation right 
% through the water column.  Throw out shallow or incomplete profiles.

rough_pres_bins = [0:200:1000];
for aa = 1:length(ocean.lon)

        tmp_prof = ocean.pres(:,aa);
        for bb = 1:length(rough_pres_bins)-1
            datacount(bb,1) = length(find(tmp_prof > rough_pres_bins(bb) & tmp_prof <= rough_pres_bins(bb+1)));
        end
        
        if min(datacount) < 2
            ocean.temp(:,aa) = nan;
            ocean.sal(:,aa) = nan;
            ocean.pres(:,aa) = nan;
        end
end


%% Sift 1 ... throw out data with < 50 (currently) observations
goodprof = zeros(1,length(ocean.lon));
for aa = 1:length(ocean.lon)
    pres_length(1,aa)  = length(find(~isnan(ocean.pres(:,aa))));
    
    if pres_length(1,aa) >= min_profile_length
        goodprof(1,aa) = 1;    
    end  
end






% Trim all variables given this criteria
ind = find(goodprof==1);
ocean.lon = ocean.lon(1,ind);
ocean.lat = ocean.lat(1,ind);
ocean.profdepth = ocean.profdepth(1,ind);
ocean.juld = ocean.juld(1,ind);
ocean.temp = ocean.temp(:,ind);
ocean.sal = ocean.sal(:,ind);
ocean.pres = ocean.pres(:,ind);
ocean.depth = ocean.depth(:,ind);
ocean.cont_vicinity = ocean.cont_vicinity(1,ind);
ocean.cont_distance = ocean.cont_distance(1,ind);
ocean.contour_lon = ocean.contour_lon(1,ind);
ocean.contour_lat = ocean.contour_lat(1,ind);
ocean.dist_from_contour = ocean.dist_from_contour(1,ind);

%% Convert to TEOS-10

disp('TEOS 10 conversions (1 min)')
% Convert to TEOS-10
% [SA, in_ocean] = gsw_SA_from_SP(SP,p,long,lat)
ocean.SA = gsw_SA_from_SP(ocean.sal,ocean.pres,ocean.lon,ocean.lat);
% CT = gsw_CT_from_t(SA,t,p)
ocean.CT = gsw_CT_from_t(ocean.SA,ocean.temp,ocean.pres);
% sigma0 = gsw_sigma0(SA,CT)
ocean.sigma0 = gsw_sigma0(ocean.SA,ocean.CT);



% Create a 'backup' temp and sal to compare; make sure the QC has worked as
% expected.
ocean.temp_NOQC = ocean.temp;
ocean.sal_NOQC = ocean.sal;
ocean.SA_NOQC = ocean.SA;
ocean.CT_NOQC = ocean.CT;
ocean.sigma0_NOQC = ocean.sigma0;



%% 1. LOCAL RANGE TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Break down profiles by distance along contour...
cont_QC_dist = 0:500:max(ocean.cont_distance); cont_QC_dist = [cont_QC_dist max(cont_QC_dist)+500];
pres_QC_bracket = 0:50:5000; pres_QC_mids = 25:50:5000;
SD_accept = 4;

% The outer QC loop
for QC_set = 1:length(cont_QC_dist)-1
    ind = find(ocean.cont_distance > cont_QC_dist(QC_set) & ocean.cont_distance <= cont_QC_dist(QC_set+1));
    
    TEMP.lon = ocean.lon(ind);
    TEMP.lat = ocean.lat(ind);
    TEMP.pres = ocean.pres(:,ind);
    TEMP.temp = ocean.temp(:,ind);
    TEMP.sal = ocean.sal(:,ind);
    TEMP.SA = ocean.SA(:,ind);
    TEMP.CT = ocean.CT(:,ind);
    TEMP.dens = ocean.sigma0(:,ind);
    
    %% Initialise some plots
    figure(1);
    clf;
    
    % TS diagram...
    ax1 = subplot(1,3,1);
    hold on;
    for aa = 1:length(TEMP.lon)
        plot(ax1,TEMP.dens(:,aa),TEMP.pres(:,aa),'.r');
    end
    set(gca,'Ydir','reverse');
    xlabel('Density - Sigma0');
    ylabel('Pres (dbar)');    
    
    ax2 = subplot(1,3,2);
    hold on;
    for aa = 1:length(TEMP.lon)
        plot(ax2,TEMP.temp(:,aa),TEMP.pres(:,aa),'.r');
    end
    set(gca,'Ydir','reverse');
    xlabel('Temperature (^oC)');
    ylabel('Pres (dbar)');    
    
    ax3 = subplot(1,3,3);
    hold on;
    for aa = 1:length(TEMP.lon)
        plot(ax3,TEMP.sal(:,aa),TEMP.pres(:,aa),'.r');
    end
    set(gca,'Ydir','reverse');
    xlabel('Salinity');
    ylabel('Pres (dbar)'); 
    
    % Set up empty arrays to handle 'good' profiles
    %goodprofs = nan * TEMP.temp;
    TEMPGOOD.temp = nan * TEMP.temp;
    TEMPGOOD.sal = nan * TEMP.temp;
    TEMPGOOD.dens = nan * TEMP.temp;
    TEMPGOOD.sal = nan * TEMP.temp;
    TEMPGOOD.SA = nan * TEMP.temp;
    TEMPGOOD.CT = nan * TEMP.temp;
    
    %% Commence pressure loop.  Step down thru pressure brackets, calculating T, S and dens SD for each.
    % Establish which profiles to keep / throw out, update goodprofs
    for pres_loop = 1:length(pres_QC_bracket)-1
        presind = find(TEMP.pres > pres_QC_bracket(pres_loop) & TEMP.pres <= pres_QC_bracket(pres_loop+1));
        T_mean = nanmean(TEMP.temp(presind)); T_mean_plot(pres_loop) = T_mean;
        T_std = nanstd(TEMP.temp(presind));
        S_mean = nanmean(TEMP.sal(presind)); S_mean_plot(pres_loop) = S_mean;
        S_std = nanstd(TEMP.sal(presind));
        dens_mean = nanmean(TEMP.dens(presind)); dens_mean_plot(pres_loop) = dens_mean;
        dens_std = nanstd(TEMP.dens(presind));
        
        % addition for case where there's only one deep profile in the QC
        % bracket
        T_std(T_std < 0.05) = 0.05;
        S_std(S_std < 0.01) = 0.01;
        dens_std(dens_std < 0.005) = 0.005;
        
        T_min = T_mean-(T_std*SD_accept); T_min_plot(pres_loop) = T_min;
        T_max = T_mean+(T_std*SD_accept); T_max_plot(pres_loop) = T_max;
        S_min = S_mean-(S_std*SD_accept); S_min_plot(pres_loop) = S_min;
        S_max = S_mean+(S_std*SD_accept); S_max_plot(pres_loop) = S_max;
        dens_min = dens_mean-(dens_std*SD_accept); dens_min_plot(pres_loop) = dens_min;
        dens_max = dens_mean+(dens_std*SD_accept); dens_max_plot(pres_loop) = dens_max;
        
        
        % Nasty fudge to throw out some scraps of profiles which otherwise
        % creep through the filter.  Will need to be modified for 1800 m
        % contour...
        if cont_QC_dist(QC_set) <= 3000 % km 
        goodind = find(TEMP.temp(presind) > T_min & TEMP.temp(presind) < T_max & ...
            TEMP.sal(presind) > S_min & TEMP.sal(presind) < S_max & ...
            TEMP.dens(presind) > dens_min & TEMP.dens(presind) < dens_max & ...
            TEMP.sal(presind) > 34.75);
        else
            goodind = find(TEMP.temp(presind) > T_min & TEMP.temp(presind) < T_max ...
                & TEMP.sal(presind) > S_min & TEMP.sal(presind) < S_max & ... 
                TEMP.dens(presind) > dens_min & TEMP.dens(presind) < dens_max);
        end
        
        % Update goodprofs
%       goodprofs(presind(goodind)) = 1;
        TEMPGOOD.temp(presind(goodind)) = TEMP.temp(presind(goodind));
        TEMPGOOD.sal(presind(goodind)) = TEMP.sal(presind(goodind));
        TEMPGOOD.SA(presind(goodind)) = TEMP.SA(presind(goodind));      
        TEMPGOOD.CT(presind(goodind)) = TEMP.CT(presind(goodind));
        TEMPGOOD.dens(presind(goodind)) = TEMP.dens(presind(goodind));
        
        
    end % end pressure loop
    
    % Now append the accepted values onto plots.
    plot(ax1,TEMPGOOD.dens,TEMP.pres,'.g');
    plot(ax2,TEMPGOOD.temp,TEMP.pres,'.g');
    plot(ax3,TEMPGOOD.sal,TEMP.pres,'.g');
    
    
    % Overlay mean line over the plots now the scatter is finished.
    plot(ax1,dens_mean_plot,pres_QC_mids,'b','linewidth',1);
    plot(ax2,T_mean_plot,pres_QC_mids,'b','linewidth',1);
    plot(ax3,S_mean_plot,pres_QC_mids,'b','linewidth',1);
    
    %plot(S_min_plot,T_min_plot,'.k','linewidth',1.5);
    %plot(S_max_plot,T_max_plot,'.k','linewidth',1.5);    

    %% print figure
    width  = 1750;  % frame width
    height = 600;  % frame height
    pngname = (['QC_plots_1000/QC_' num2str(cont_QC_dist(QC_set)) ' km to ' num2str(cont_QC_dist(QC_set+1)) 'km.png']);
    
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
    
    clear T_mean_plot S_mean_plot T_min_plot T_max_plot S_min_plot S_max_plot
    
    %% Most important bit: Replace ocean.temp with a QC'd version.  Keep the old for comparison
    ocean.temp(:,ind) = TEMPGOOD.temp;
    ocean.sal(:,ind) = TEMPGOOD.sal;
    ocean.SA(:,ind) = TEMPGOOD.SA;
    ocean.CT(:,ind) = TEMPGOOD.CT;
    ocean.sigma0(:,ind) = TEMPGOOD.dens;
    
    
    clear TEMP TEMPGOOD
    disp(['Completed QC set number: ' num2str(QC_set)])
end




%% plot T and S for ALL DATA

figure(2);
clf;


ax1 = subplot(1,3,1);
hold on;
for aa = 1:length(ocean.lon)
    plot(ocean.sigma0_NOQC(:,aa),ocean.pres(:,aa),'.r');
end

for aa = 1:length(ocean.lon)
    plot(ocean.sigma0(:,aa),ocean.pres(:,aa),'.g');
end

set(gca,'Ydir','reverse');
xlabel('Density - Sigma0');
ylabel('Pres (dbar)');

ax2 = subplot(1,3,2);
hold on;
for aa = 1:length(ocean.lon)
    plot(ocean.temp_NOQC(:,aa),ocean.pres(:,aa),'.r');
end

for aa = 1:length(ocean.lon)
    plot(ocean.temp(:,aa),ocean.pres(:,aa),'.g');
end

set(gca,'Ydir','reverse');
xlabel('Temperature (^oC)');
ylabel('Pres (dbar)');

ax3 = subplot(1,3,3);
hold on;
for aa = 1:length(ocean.lon)
    plot(ocean.sal_NOQC(:,aa),ocean.pres(:,aa),'.r');
end

for aa = 1:length(ocean.lon)
    plot(ocean.sal(:,aa),ocean.pres(:,aa),'.g');
end

set(gca,'Ydir','reverse');
xlabel('Salinity');
ylabel('Pres (dbar)');


%% print figure
width  = 1750;  % frame width
height = 600;  % frame height
pngname = (['QC_plots_1000/QC_ALL.png']);

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





%% 2: DENSITY INVERSION TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pres_grid_rough = 20:20:5000;
cc = 1; dd = 1;
badprof = 0;

% density inversion figure
figure(40);
hold on;

for aa = 1:length(ocean.lat) % for each profile
    
    pres_temp = ocean.pres(:,aa);
    sigma0_temp = ocean.sigma0(:,aa);
    
    ind = find(~isnan(sigma0_temp));
    if length(ind) > min_profile_length-10
        % Vq = interp1(X,V,Xq)
        dens_interp = interp1(pres_temp(ind),sigma0_temp(ind),pres_grid_rough);
        
        % are there inversions in this profile?
        for bb = 1:length(dens_interp)-1
            thisobs = dens_interp(bb);
            
            minaccept = thisobs - 0.03; % kg/m3 maximum density inversion permitted
            nextobs = dens_interp(bb+1);
            
            if nextobs <= minaccept % if there is an inversion
                invers(1,cc) = aa; % counter
                %baddata(:,cc) = dens_interp';
                

                
                badprof = 1; % 
                cc = cc+1;
                
            end
        end
    else % if there are more than a few gaps in the profile we don't want it anyway...
        toogappy(1,cc) = aa;
        ocean.temp(:,aa) = nan;
        ocean.sal(:,aa) = nan;
        ocean.pres(:,aa) = nan;
        ocean.depth(:,aa) = nan;
        ocean.CT(:,aa) = nan;
        ocean.SA(:,aa) = nan;
        ocean.sigma0(:,aa) = nan;
        dd = dd+1;
    end
    
    if badprof==1
      
        % Add to density inversion figure
        plot(ocean.sigma0(:,aa),ocean.pres(:,aa),'r');
        plot(ocean.sigma0(:,aa),ocean.pres(:,aa),'.r');    
        
        ocean.temp(:,aa) = nan;
        ocean.sal(:,aa) = nan;
        ocean.pres(:,aa) = nan;
        ocean.depth(:,aa) = nan;
        ocean.CT(:,aa) = nan;
        ocean.SA(:,aa) = nan;
        ocean.sigma0(:,aa) = nan;
        
        badprof = 0;
    end
end

% Tidy up and print figure
set(gca,'Ydir','reverse');
xlabel('Density - Sigma0');
ylabel('Pres (dbar)');
title('Profiles rejected due to density inversions');

%% print figure
width  = 1000;  % frame width
height = 1000;  % frame height
pngname = (['QC_plots_1000/Dens_inversion_fails.png']);

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


disp(['Profiles thrown out due to gappiness: ' num2str(length(toogappy))]);
disp(['Profiles thrown out due to dens inversions: ' num2str(length(invers))]);




%% Temporary save for sanity purposes.
save('D:\Work_computer_sync\OSNAP_postdoc\PAPERS_NEW\N_Atlantic_boundary\matlab\V3_050321\intermediate_saves\raw_profs_1000m_P2_QC','ocean','-v7.3');
