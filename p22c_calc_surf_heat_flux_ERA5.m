% Calc_surf_heat_flux
% SJones 0322
clear; close all;

addpath(genpath('D:\Work_computer_sync\MATLAB_functions'));


% Load 1000m contour data
C = load('D:\Work_computer_sync\OSNAP_postdoc\PAPERS_NEW\N_Atlantic_boundary\matlab\V3_050321\intermediate_saves/contour_data_1000_EN4inserted.mat');

year_range = 2000:2019;


%% Load ERA5
% https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=form
% Load global vars
filename = 'ERA5/HF_160522.nc';
lon = double(ncread(filename,'longitude'));
lat = double(ncread(filename,'latitude'));

basetime = (datenum('01-01-1900'));
time = double(ncread(filename,'time'));
time = (time/24)+basetime;


latent = ncread(filename,'slhf'); latent = squeeze(latent(:,:,1,:));
sensible = ncread(filename,'sshf'); sensible = squeeze(sensible(:,:,1,:));
sw = ncread(filename,'ssr'); sw = squeeze(sw(:,:,1,:));
lw = ncread(filename,'str'); lw = squeeze(lw(:,:,1,:));





% precip(precip < 0 | precip > 200) = nan;

%% Subset the 2000-2020 period
timeind = find(time >= datenum(['01-01-' num2str(year_range(1))]) & time <= datenum(['12-31-' num2str(year_range(end))]));

time = time(timeind);
latent = latent(:,:,timeind);
sensible = sensible(:,:,timeind);
sw = sw(:,:,timeind);
lw = lw(:,:,timeind);







%% One day accumulation period.  We want in seconds to convert to Watts.
latent = latent./(3600*24);
sensible = sensible./(3600*24);
sw = sw./(3600*24);
lw = lw./(3600*24);

%% Get mean periods
month = nan*time;
year = nan*time;
for aa = 1:length(time)
month(aa) = str2num(datestr(time(aa),'mm'));
year(aa) = str2num(datestr(time(aa),'yyyy'));
end

latent_mean = nanmean(latent,3); 
sensible_mean = nanmean(sensible,3); 
sw_mean = nanmean(sw,3);
lw_mean = nanmean(lw,3);


%% Do a separate subset to get annual means for error estimation

for aa = 1:length(year_range)
    ind = find(year == year_range(aa));
    latent_mean_y(:,:,aa) = nanmean(latent(:,:,ind),3);
    sensible_mean_y(:,:,aa) = nanmean(sensible(:,:,ind),3);
    sw_mean_y(:,:,aa) = nanmean(sw(:,:,ind),3);
    lw_mean_y(:,:,aa) = nanmean(lw(:,:,ind),3);

    ind = find(year == year_range(aa) & (month == 1 | month == 2 | month == 3));
    latent_mean_y_JFM(:,:,aa) = nanmean(latent(:,:,ind),3);
    sensible_mean_y_JFM(:,:,aa) = nanmean(sensible(:,:,ind),3);
    sw_mean_y_JFM(:,:,aa) = nanmean(sw(:,:,ind),3);
    lw_mean_y_JFM(:,:,aa) = nanmean(lw(:,:,ind),3);

    ind = find(year == year_range(aa) & (month == 4 | month == 5 | month == 6));
    latent_mean_y_AMJ(:,:,aa) = nanmean(latent(:,:,ind),3);
    sensible_mean_y_AMJ(:,:,aa) = nanmean(sensible(:,:,ind),3);
    sw_mean_y_AMJ(:,:,aa) = nanmean(sw(:,:,ind),3);
    lw_mean_y_AMJ(:,:,aa) = nanmean(lw(:,:,ind),3);

        ind = find(year == year_range(aa) & (month == 7 | month == 8 | month == 9));
    latent_mean_y_JAS(:,:,aa) = nanmean(latent(:,:,ind),3);
    sensible_mean_y_JAS(:,:,aa) = nanmean(sensible(:,:,ind),3);
    sw_mean_y_JAS(:,:,aa) = nanmean(sw(:,:,ind),3);
    lw_mean_y_JAS(:,:,aa) = nanmean(lw(:,:,ind),3);

        ind = find(year == year_range(aa) & (month == 10 | month == 11 | month == 12));
    latent_mean_y_OND(:,:,aa) = nanmean(latent(:,:,ind),3);
    sensible_mean_y_OND(:,:,aa) = nanmean(sensible(:,:,ind),3);
    sw_mean_y_OND(:,:,aa) = nanmean(sw(:,:,ind),3);
    lw_mean_y_OND(:,:,aa) = nanmean(lw(:,:,ind),3);

end

Total_y = latent_mean_y + sensible_mean_y + sw_mean_y + lw_mean_y; Total_stE = (nanstd(Total_y,0,3) ./ sqrt(length(Total_y(1,1,:))));

Total_y_JFM = latent_mean_y_JFM + sensible_mean_y_JFM + sw_mean_y_JFM + lw_mean_y_JFM; Total_stE_JFM = (nanstd(Total_y_JFM,0,3) ./ sqrt(length(Total_y_JFM(1,1,:))));
Total_y_AMJ = latent_mean_y_AMJ + sensible_mean_y_AMJ + sw_mean_y_AMJ + lw_mean_y_AMJ; Total_stE_AMJ = (nanstd(Total_y_AMJ,0,3) ./ sqrt(length(Total_y_AMJ(1,1,:))));
Total_y_JAS = latent_mean_y_JAS + sensible_mean_y_JAS + sw_mean_y_JAS + lw_mean_y_JAS; Total_stE_JAS = (nanstd(Total_y_JAS,0,3) ./ sqrt(length(Total_y_JAS(1,1,:))));
Total_y_OND = latent_mean_y_OND + sensible_mean_y_OND + sw_mean_y_OND + lw_mean_y_OND; Total_stE_OND = (nanstd(Total_y_OND,0,3) ./ sqrt(length(Total_y_OND(1,1,:))));

% latent_std = nanstd(latent,0,3); 
% sensible_std = nanstd(sensible,0,3); 
% sw_std = nanstd(sw,0,3);
% lw_std = nanstd(lw,0,3);


% Do the same for each season
ind = find(month == 1 | month == 2 | month == 3);
latent_mean_JFM = nanmean(latent(:,:,ind),3); latent_std_JFM = nanstd(latent(:,:,ind),0,3);
sensible_mean_JFM = nanmean(sensible(:,:,ind),3); sensible_std_JFM = nanstd(sensible(:,:,ind),0,3); 
sw_mean_JFM = nanmean(sw(:,:,ind),3); sw_std_JFM = nanstd(sw(:,:,ind),0,3);
lw_mean_JFM = nanmean(lw(:,:,ind),3); lw_std_JFM = nanstd(lw(:,:,ind),0,3); 

ind = find(month == 4 | month == 5 | month == 6);
latent_mean_AMJ = nanmean(latent(:,:,ind),3); latent_std_AMJ = nanstd(latent(:,:,ind),0,3);
sensible_mean_AMJ = nanmean(sensible(:,:,ind),3); sensible_std_AMJ = nanstd(sensible(:,:,ind),0,3);
sw_mean_AMJ = nanmean(sw(:,:,ind),3); sw_std_AMJ = nanstd(sw(:,:,ind),0,3);
lw_mean_AMJ = nanmean(lw(:,:,ind),3); lw_std_AMJ = nanstd(lw(:,:,ind),0,3);

ind = find(month == 7 | month == 8 | month == 9);
latent_mean_JAS = nanmean(latent(:,:,ind),3); latent_std_JAS = nanstd(latent(:,:,ind),0,3);
sensible_mean_JAS = nanmean(sensible(:,:,ind),3); sensible_std_JAS = nanstd(sensible(:,:,ind),0,3);
sw_mean_JAS = nanmean(sw(:,:,ind),3); sw_std_JAS = nanstd(sw(:,:,ind),0,3);
lw_mean_JAS = nanmean(lw(:,:,ind),3); lw_std_JAS = nanstd(lw(:,:,ind),0,3); 

ind = find(month == 10 | month == 11 | month == 12);
latent_mean_OND = nanmean(latent(:,:,ind),3); latent_std_OND = nanstd(latent(:,:,ind),0,3);
sensible_mean_OND = nanmean(sensible(:,:,ind),3); sensible_std_OND = nanstd(sensible(:,:,ind),0,3);
sw_mean_OND = nanmean(sw(:,:,ind),3); sw_std_OND = nanstd(sw(:,:,ind),0,3);
lw_mean_OND = nanmean(lw(:,:,ind),3); lw_std_OND = nanstd(lw(:,:,ind),0,3); 




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







latent_mean(SPG_mask == 0) = nan;
sensible_mean(SPG_mask == 0) = nan;
sw_mean(SPG_mask == 0) = nan;
lw_mean(SPG_mask == 0) = nan;
latent_std(SPG_mask == 0) = nan;
sensible_std(SPG_mask == 0) = nan;
sw_std(SPG_mask == 0) = nan;
lw_std(SPG_mask == 0) = nan;

latent_mean_JFM(SPG_mask == 0) = nan;
sensible_mean_JFM(SPG_mask == 0) = nan;
sw_mean_JFM(SPG_mask == 0) = nan;
lw_mean_JFM(SPG_mask == 0) = nan;
latent_std_JFM(SPG_mask == 0) = nan;
sensible_std_JFM(SPG_mask == 0) = nan;
sw_std_JFM(SPG_mask == 0) = nan;
lw_std_JFM(SPG_mask == 0) = nan;

latent_mean_AMJ(SPG_mask == 0) = nan;
sensible_mean_AMJ(SPG_mask == 0) = nan;
sw_mean_AMJ(SPG_mask == 0) = nan;
lw_mean_AMJ(SPG_mask == 0) = nan;
latent_std_AMJ(SPG_mask == 0) = nan;
sensible_std_AMJ(SPG_mask == 0) = nan;
sw_std_AMJ(SPG_mask == 0) = nan;
lw_std_AMJ(SPG_mask == 0) = nan;

latent_mean_JAS(SPG_mask == 0) = nan;
sensible_mean_JAS(SPG_mask == 0) = nan;
sw_mean_JAS(SPG_mask == 0) = nan;
lw_mean_JAS(SPG_mask == 0) = nan;
latent_std_JAS(SPG_mask == 0) = nan;
sensible_std_JAS(SPG_mask == 0) = nan;
sw_std_JAS(SPG_mask == 0) = nan;
lw_std_JAS(SPG_mask == 0) = nan;

latent_mean_OND(SPG_mask == 0) = nan;
sensible_mean_OND(SPG_mask == 0) = nan;
sw_mean_OND(SPG_mask == 0) = nan;
lw_mean_OND(SPG_mask == 0) = nan;
latent_std_OND(SPG_mask == 0) = nan;
sensible_std_OND(SPG_mask == 0) = nan;
sw_std_OND(SPG_mask == 0) = nan;
lw_std_OND(SPG_mask == 0) = nan;

area_m2(SPG_mask == 0) = nan;

Total_stE(SPG_mask == 0) = nan;
Total_stE_JFM(SPG_mask == 0) = nan;
Total_stE_AMJ(SPG_mask == 0) = nan;
Total_stE_JAS(SPG_mask == 0) = nan;
Total_stE_OND(SPG_mask == 0) = nan;







TOTAL = latent_mean+sensible_mean+sw_mean+lw_mean;
% TOTALstd = latent_std+sensible_std+sw_std+lw_std;

TOTAL_JFM = latent_mean_JFM+sensible_mean_JFM+sw_mean_JFM+lw_mean_JFM;
TOTAL_AMJ = latent_mean_AMJ+sensible_mean_AMJ+sw_mean_AMJ+lw_mean_AMJ;
TOTAL_JAS = latent_mean_JAS+sensible_mean_JAS+sw_mean_JAS+lw_mean_JAS;
TOTAL_OND = latent_mean_OND+sensible_mean_OND+sw_mean_OND+lw_mean_OND;

% TOTALstd_JFM = latent_std_JFM+sensible_std_JFM+sw_std_JFM+lw_std_JFM;
% TOTALstd_AMJ = latent_std_AMJ+sensible_std_AMJ+sw_std_AMJ+lw_std_AMJ;
% TOTALstd_JAS = latent_std_JAS+sensible_std_JAS+sw_std_JAS+lw_std_JAS;
% TOTALstd_OND = latent_std_OND+sensible_std_OND+sw_std_OND+lw_std_OND;






%% TOTAL figure
figure;
subplot(2,2,1);
hold on;
lonmin = -80; lonmax = 0; latmin = 45; latmax = 70;
pcolor(lon,lat,TOTAL_JFM'); shading flat
xlim([lonmin lonmax]);
ylim([latmin latmax]);
plot(C.cont_lon,C.cont_lat,'r','linewidth',2);
title('J m-2 JFM');
colorbar;
caxis([-150 150]);

subplot(2,2,2);
hold on;
lonmin = -80; lonmax = 0; latmin = 45; latmax = 70;
pcolor(lon,lat,TOTAL_AMJ'); shading flat
xlim([lonmin lonmax]);
ylim([latmin latmax]);
plot(C.cont_lon,C.cont_lat,'r','linewidth',2);
title('J m-2 AMJ');
colorbar;
caxis([-150 150]);

subplot(2,2,3);
hold on;
lonmin = -80; lonmax = 0; latmin = 45; latmax = 70;
pcolor(lon,lat,TOTAL_JAS'); shading flat
xlim([lonmin lonmax]);
ylim([latmin latmax]);
plot(C.cont_lon,C.cont_lat,'r','linewidth',2);
title('J m-2 JAS');
colorbar;
caxis([-150 150]);

subplot(2,2,4);
hold on;
lonmin = -80; lonmax = 0; latmin = 45; latmax = 70;
pcolor(lon,lat,TOTAL_OND'); shading flat
xlim([lonmin lonmax]);
ylim([latmin latmax]);
plot(C.cont_lon,C.cont_lat,'r','linewidth',2);
title('J m-2 OND');
colorbar;
caxis([-150 150]);




%% Annual means figure
% latent_mean(SPG_mask == 0) = nan;
% sensible_mean(SPG_mask == 0) = nan;
% sw_mean(SPG_mask == 0) = nan;
% lw_mean(SPG_mask == 0) = nan;

figure;
subplot(2,2,1);
hold on;
lonmin = -80; lonmax = 0; latmin = 45; latmax = 70;
pcolor(lon,lat,latent_mean'); shading flat
xlim([lonmin lonmax]);
ylim([latmin latmax]);
plot(C.cont_lon,C.cont_lat,'r','linewidth',2);
title('J m-2 Latent annual mean');
colorbar;
caxis([-120 0]);

subplot(2,2,2);
hold on;
lonmin = -80; lonmax = 0; latmin = 45; latmax = 70;
pcolor(lon,lat,sensible_mean'); shading flat
xlim([lonmin lonmax]);
ylim([latmin latmax]);
plot(C.cont_lon,C.cont_lat,'r','linewidth',2);
title('J m-2 sensible annual mean');
colorbar;
caxis([-60 0]);

subplot(2,2,3);
hold on;
lonmin = -80; lonmax = 0; latmin = 45; latmax = 70;
pcolor(lon,lat,sw_mean'); shading flat
xlim([lonmin lonmax]);
ylim([latmin latmax]);
plot(C.cont_lon,C.cont_lat,'r','linewidth',2);
title('J m-2 SW annual mean');
colorbar;
caxis([50 120]);

subplot(2,2,4);
hold on;
lonmin = -80; lonmax = 0; latmin = 45; latmax = 70;
pcolor(lon,lat,lw_mean'); shading flat
xlim([lonmin lonmax]);
ylim([latmin latmax]);
plot(C.cont_lon,C.cont_lat,'r','linewidth',2);
title('J m-2 LW annual mean');
colorbar;
caxis([-60 0]);









%% Multiply up to get totals

area_tot_m = nansum(nansum(area_m2));


TOTAL_net = (nansum(nansum(TOTAL.*area_m2))) ./ 1e15;

JFM_net = (nansum(nansum(TOTAL_JFM.*area_m2))) ./ 1e15;
AMJ_net = (nansum(nansum(TOTAL_AMJ.*area_m2))) ./ 1e15;
JAS_net = (nansum(nansum(TOTAL_JAS.*area_m2))) ./ 1e15;
OND_net = (nansum(nansum(TOTAL_OND.*area_m2))) ./ 1e15;



%% Uncertainties
Total_stE = (nansum(nansum(Total_stE.*area_m2))) ./ 1e15;
Total_stE_JFM = (nansum(nansum(Total_stE_JFM.*area_m2))) ./ 1e15;
Total_stE_AMJ = (nansum(nansum(Total_stE_AMJ.*area_m2))) ./ 1e15;
Total_stE_JAS = (nansum(nansum(Total_stE_JAS.*area_m2))) ./ 1e15;
Total_stE_OND = (nansum(nansum(Total_stE_OND.*area_m2))) ./ 1e15;


% %% Rainfall and evap totals for text
% % P_E_m = (nansum(nansum(P_E.*area_m2))) ./ area_tot_m;
% P_m = (nansum(nansum(prec_mean.*area_m2))) ./ area_tot_m;
% E_m = (nansum(nansum(evap_mean.*area_m2))) ./ area_tot_m;

disp('Surface flux (PW) *************');
disp(['Annual mean: ' num2str(TOTAL_net) ' +/- ' num2str(Total_stE)]);
disp(['JFM: ' num2str(JFM_net) ' +/- ' num2str(Total_stE_JFM)]);
disp(['AMJ: ' num2str(AMJ_net) ' +/- ' num2str(Total_stE_AMJ)]);
disp(['JAS: ' num2str(JAS_net) ' +/- ' num2str(Total_stE_JAS)]);
disp(['OND: ' num2str(OND_net) ' +/- ' num2str(Total_stE_OND)]);







