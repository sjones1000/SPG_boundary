% p22c_calc_surf_PminusE_ERA5
% SJones 04/22
clear; close all;

addpath(genpath('D:\Work_computer_sync\MATLAB_functions'));


% Load 1000m contour data
C = load('D:\Work_computer_sync\OSNAP_postdoc\PAPERS_NEW\N_Atlantic_boundary\matlab\V3_050321\intermediate_saves/contour_data_1000_EN4inserted.mat');

year_range = 2000:2019;


%% Load precip, evap ERA5
% https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=form
% Load global vars
filename = 'ERA5/PE_260422.nc';
lon = double(ncread(filename,'longitude'));
lat = double(ncread(filename,'latitude'));

basetime = (datenum('01-01-1900'));
time = double(ncread(filename,'time'));
time = (time/24)+basetime;


precip = ncread(filename,'tp'); precip = squeeze(precip(:,:,1,:));
evap = ncread(filename,'e'); evap = squeeze(evap(:,:,1,:));
% precip(precip < 0 | precip > 200) = nan;

%% Subset the 2000-2020 period
timeind = find(time >= datenum(['01-01-' num2str(year_range(1))]) & time <= datenum(['12-31-' num2str(year_range(end))]));

time = time(timeind);
precip = precip(:,:,timeind);
evap = evap(:,:,timeind);


%% Get mean periods
month = nan*time;
for aa = 1:length(time)
month(aa) = str2num(datestr(time(aa),'mm'));
end

prec_mean = nanmean(precip,3); 
prec_mean = prec_mean.*365.25;

evap_mean = nanmean(evap,3); 
evap_mean = evap_mean.*365.25;


% Do the same for each season
ind = find(month == 1 | month == 2 | month == 3);
prec_mean_JFM = nanmean(precip(:,:,ind),3); prec_mean_JFM = prec_mean_JFM.*365.25;
evap_mean_JFM = nanmean(evap(:,:,ind),3); evap_mean_JFM = evap_mean_JFM.*365.25;

ind = find(month == 4 | month == 5 | month == 6);
prec_mean_AMJ = nanmean(precip(:,:,ind),3); prec_mean_AMJ = prec_mean_AMJ.*365.25;
evap_mean_AMJ = nanmean(evap(:,:,ind),3); evap_mean_AMJ = evap_mean_AMJ.*365.25;

ind = find(month == 7 | month == 8 | month == 9);
prec_mean_JAS = nanmean(precip(:,:,ind),3); prec_mean_JAS = prec_mean_JAS.*365.25;
evap_mean_JAS = nanmean(evap(:,:,ind),3); evap_mean_JAS = evap_mean_JAS.*365.25;

ind = find(month == 10 | month == 11 | month == 12);
prec_mean_OND = nanmean(precip(:,:,ind),3); prec_mean_OND = prec_mean_OND.*365.25;
evap_mean_OND = nanmean(evap(:,:,ind),3); evap_mean_OND = evap_mean_OND.*365.25;




%% Get area associated with each cell
dlon = diff(lon); dlon = dlon(1);
dlat = diff(lat); dlat = -dlat(1);

% metres in degree of latitude and longitude
for aa = 1:length(lon)
    for bb = 1:length(lat)

    dlat_m(aa,bb) = dlat .* 111139;
    dlon_m(aa,bb) = dlon .* (111139*cosd(lat(bb))) ;

    end
end

area_m2 = dlon_m.*dlat_m;



%% only interested in region within SPG

[lonm,latm] = meshgrid(lon,lat);
NX = size(lonm,1);
NY = size(lonm,2);
veclength = NX*NY;

lonvec = reshape(lonm,1,veclength);
latvec = reshape(latm,1,veclength);

% IN = inpolygon(X,Y,X_boundary,Y_boundary)
IN = inpolygon(lonvec,latvec,C.cont_lon,C.cont_lat);

SPG_mask = double(reshape(IN,NX,NY))';







prec_mean(SPG_mask == 0) = nan;
evap_mean(SPG_mask == 0) = nan;

prec_mean_JFM(SPG_mask == 0) = nan;
evap_mean_JFM(SPG_mask == 0) = nan;
prec_mean_AMJ(SPG_mask == 0) = nan;
evap_mean_AMJ(SPG_mask == 0) = nan;
prec_mean_JAS(SPG_mask == 0) = nan;
evap_mean_JAS(SPG_mask == 0) = nan;
prec_mean_OND(SPG_mask == 0) = nan;
evap_mean_OND(SPG_mask == 0) = nan;

area_m2(SPG_mask == 0) = nan;

P_E = prec_mean+evap_mean;

P_E_JFM = prec_mean_JFM+evap_mean_JFM;
P_E_AMJ = prec_mean_AMJ+evap_mean_AMJ;
P_E_JAS = prec_mean_JAS+evap_mean_JAS;
P_E_OND = prec_mean_OND+evap_mean_OND;






%% Precip figure
figure;
subplot(2,2,1);
hold on;
lonmin = -80; lonmax = 0; latmin = 45; latmax = 70;
pcolor(lon,lat,prec_mean_JFM'); shading flat
xlim([lonmin lonmax]);
ylim([latmin latmax]);
plot(C.cont_lon,C.cont_lat,'r','linewidth',2);
title('prec (m) JFM');
colorbar;
caxis([0 2]);

subplot(2,2,2);
hold on;
lonmin = -80; lonmax = 0; latmin = 45; latmax = 70;
pcolor(lon,lat,prec_mean_AMJ'); shading flat
xlim([lonmin lonmax]);
ylim([latmin latmax]);
plot(C.cont_lon,C.cont_lat,'r','linewidth',2);
title('prec (m) AMJ');
colorbar;
caxis([0 2]);

subplot(2,2,3);
hold on;
lonmin = -80; lonmax = 0; latmin = 45; latmax = 70;
pcolor(lon,lat,prec_mean_JAS'); shading flat
xlim([lonmin lonmax]);
ylim([latmin latmax]);
plot(C.cont_lon,C.cont_lat,'r','linewidth',2);
title('prec (m) JAS');
colorbar;
caxis([0 2]);

subplot(2,2,4);
hold on;
lonmin = -80; lonmax = 0; latmin = 45; latmax = 70;
pcolor(lon,lat,prec_mean_OND'); shading flat
xlim([lonmin lonmax]);
ylim([latmin latmax]);
plot(C.cont_lon,C.cont_lat,'r','linewidth',2);
title('prec (m) OND');
colorbar;
caxis([0 2]);






%% Evap figure
figure;
subplot(2,2,1);
hold on;
lonmin = -80; lonmax = 0; latmin = 45; latmax = 70;
pcolor(lon,lat,evap_mean_JFM'); shading flat
xlim([lonmin lonmax]);
ylim([latmin latmax]);
plot(C.cont_lon,C.cont_lat,'r','linewidth',2);
title('evap (m) JFM');
colorbar;
caxis([-2 0]);

subplot(2,2,2);
hold on;
lonmin = -80; lonmax = 0; latmin = 45; latmax = 70;
pcolor(lon,lat,evap_mean_AMJ'); shading flat
xlim([lonmin lonmax]);
ylim([latmin latmax]);
plot(C.cont_lon,C.cont_lat,'r','linewidth',2);
title('evap (m) AMJ');
colorbar;
caxis([-2 0]);

subplot(2,2,3);
hold on;
lonmin = -80; lonmax = 0; latmin = 45; latmax = 70;
pcolor(lon,lat,evap_mean_JAS'); shading flat
xlim([lonmin lonmax]);
ylim([latmin latmax]);
plot(C.cont_lon,C.cont_lat,'r','linewidth',2);
title('evap (m) JAS');
colorbar;
caxis([-2 0]);

subplot(2,2,4);
hold on;
lonmin = -80; lonmax = 0; latmin = 45; latmax = 70;
pcolor(lon,lat,evap_mean_OND'); shading flat
xlim([lonmin lonmax]);
ylim([latmin latmax]);
plot(C.cont_lon,C.cont_lat,'r','linewidth',2);
title('evap (m) OND');
colorbar;
caxis([-2 0]);








%% P-E figure
figure;
subplot(2,2,1);
hold on;
lonmin = -80; lonmax = 0; latmin = 45; latmax = 70;
pcolor(lon,lat,P_E_JFM'); shading flat
xlim([lonmin lonmax]);
ylim([latmin latmax]);
plot(C.cont_lon,C.cont_lat,'r','linewidth',2);
title('P+E (m) JFM');
colorbar;
caxis([0 1]);

subplot(2,2,2);
hold on;
lonmin = -80; lonmax = 0; latmin = 45; latmax = 70;
pcolor(lon,lat,P_E_AMJ'); shading flat
xlim([lonmin lonmax]);
ylim([latmin latmax]);
plot(C.cont_lon,C.cont_lat,'r','linewidth',2);
title('P+E (m) AMJ');
colorbar;
caxis([0 1]);

subplot(2,2,3);
hold on;
lonmin = -80; lonmax = 0; latmin = 45; latmax = 70;
pcolor(lon,lat,P_E_JAS'); shading flat
xlim([lonmin lonmax]);
ylim([latmin latmax]);
plot(C.cont_lon,C.cont_lat,'r','linewidth',2);
title('P+E (m) JAS');
colorbar;
caxis([0 1]);

subplot(2,2,4);
hold on;
lonmin = -80; lonmax = 0; latmin = 45; latmax = 70;
pcolor(lon,lat,P_E_OND'); shading flat
xlim([lonmin lonmax]);
ylim([latmin latmax]);
plot(C.cont_lon,C.cont_lat,'r','linewidth',2);
title('P+E (m) OND');
colorbar;
caxis([0 1]);



%% Multiply up to get totals

area_tot_m = nansum(nansum(area_m2));

% Get PE in m2/s for each cell
PE_m2s = (P_E.*area_m2) ./ (365.25 * 3600 * 24);

PE_m2s_JFM = (P_E_JFM.*area_m2) ./ (365.25 * 3600 * 24);
PE_m2s_AMJ = (P_E_AMJ.*area_m2) ./ (365.25 * 3600 * 24);
PE_m2s_JAS = (P_E_JAS.*area_m2) ./ (365.25 * 3600 * 24);
PE_m2s_OND = (P_E_OND.*area_m2) ./ (365.25 * 3600 * 24);

% total
PE_sv_annual =  (nansum(nansum(PE_m2s)))./1e6;

PE_sv_JFM =  (nansum(nansum(PE_m2s_JFM)))./1e6;
PE_sv_AMJ =  (nansum(nansum(PE_m2s_AMJ)))./1e6;
PE_sv_JAS =  (nansum(nansum(PE_m2s_JAS)))./1e6;
PE_sv_OND =  (nansum(nansum(PE_m2s_OND)))./1e6;


%% Rainfall and evap totals for text
% P_E_m = (nansum(nansum(P_E.*area_m2))) ./ area_tot_m;
P_m = (nansum(nansum(prec_mean.*area_m2))) ./ area_tot_m;
E_m = (nansum(nansum(evap_mean.*area_m2))) ./ area_tot_m;



disp(['Precipitation mean= ' num2str(P_m)]);
disp(['Evaporation mean= ' num2str(E_m)]);
disp(' ');
disp('Surface flux (Sv) *************');
disp(['Annual mean: ' num2str(PE_sv_annual)]);
disp(['JFM: ' num2str(PE_sv_JFM)]);
disp(['AMJ: ' num2str(PE_sv_AMJ)]);
disp(['JAS: ' num2str(PE_sv_JAS)]);
disp(['OND: ' num2str(PE_sv_OND)]);



