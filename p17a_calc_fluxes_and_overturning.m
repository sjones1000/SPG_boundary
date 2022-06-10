% p17a_calc_fluxes_and_overturning
% Compute fluxes and overturning
% SJones 06/21
clear; close all;

addpath(genpath('D:\Work_computer_sync\MATLAB_functions'));

load D:\Work_computer_sync\OSNAP_postdoc\PAPERS_NEW\N_Atlantic_boundary\matlab\V3_050321\intermediate_saves\SPG_geovels_inc_EKMAN.mat

% Specific heat capacity
% cp0 = gsw_cp0
% take rho.cp as 4.1*10^6, e.g. Lozier et al., 2019 
rhoCp = 4.1e6;

%% set positions of distance ticks for locations of interest
% 2300==FSC, 3350=Iceland, 4100=Rekjarnes 6000=Greenland 7050=Lab Sea,
% 9900=Gulf Stream
distance_ticks = [2300 3350 4100 6000 7050 9900];
OSNAP_crossings = [1958 5820 6105 8528];

area_metres2_3D = repmat(area_metres2,1,1,4);


%% Compute means and anomolies

% Mean velocity
vel = cat(3,TOTAL_above_1000m_inc_EKMAN.vel_JFM,TOTAL_above_1000m_inc_EKMAN.vel_AMJ,TOTAL_above_1000m_inc_EKMAN.vel_JAS,TOTAL_above_1000m_inc_EKMAN.vel_OND);


%% REDUCE GULF STREAM LINEARLY TO ZERO AT BED?
depth_range = length(52:243); % depth index to operate on
vel_multiplier = linspace(1,0,depth_range)';
dist_range = length(65:71); % distance index to operate on
vel_multiplier = repmat(vel_multiplier,1,dist_range);
vel(52:243,65:71,:) = vel(52:243,65:71,:) .* vel_multiplier;


vel_mean = nanmean(vel,3);

% Full depth transport from velocity
tran = vel .* area_metres2_3D;
tran_mean = vel_mean .* area_metres2;

% Temp
T = cat(3,WATER_PROPERTIES.CT_JFM,WATER_PROPERTIES.CT_AMJ,WATER_PROPERTIES.CT_JAS,WATER_PROPERTIES.CT_OND);
T_mean = nanmean(T,3);

% Sal
S = cat(3,WATER_PROPERTIES.SA_JFM,WATER_PROPERTIES.SA_AMJ,WATER_PROPERTIES.SA_JAS,WATER_PROPERTIES.SA_OND);
S_mean = nanmean(S,3);

%% Add mass balanced velocity [Needs to be done on whole (sub-1000 m section)] 
TRAN = nansum(nansum(tran_mean));
AREA = nansum(nansum(area_metres2(52:end,:)));
TRAN_seas = squeeze(nansum(nansum(tran)));


% mean velocity on section
vref = TRAN / AREA;
vref_seas = TRAN_seas ./ AREA;

vel_mean_ADJ = vel_mean; % preallocate
vel_seas_ADJ = vel;

vel_mean_ADJ(52:end,:) = vel_mean_ADJ(52:end,:) - vref;
% Make vel adjustments in a loop as indexing looks a bit sketchy here
for aa = 1:4 % for each season
vel_seas_ADJ(52:end,:,aa) = vel_seas_ADJ(52:end,:,aa) - vref_seas(aa);
end

tran_mean_ADJ = vel_mean_ADJ .* area_metres2;
tran_seas_ADJ = vel_seas_ADJ .* area_metres2_3D;
% To Check! 
% nansum(nansum(tran_mean_ADJ));

%% Heat flux in Watts
% V * temp * rho.Cp0 * area
% SJ 010222: Modified to include a T_REF: 3.9970    3.9855    4.0581    4.0740 (MEAN: 4.0286)
Tref = 4.0286;
HF_mean = vel_mean_ADJ .* (T_mean - Tref) .* rhoCp .* area_metres2;
% Now seasonal
HF_seas = vel_seas_ADJ .* (T - Tref) .* rhoCp .* area_metres2_3D;
% HF_anom = HF_seas - HF_mean; 

%% Freshwater flux
% Define Sbar (mean salinity on section in this case)
%Sref = nanmean(nanmean(S_mean)); 
% Sref = nanmax(nanmax(nanmax(S)));
% Seasonal Sref = 35.1372   35.1359   35.1334   35.1344
Sref = 35.1352;
%Sref = nanmean(nanmean(S_mean(1:51,:)));
% freshwater flux = - ((vel * S) / S) * area
%FF_mean = (((vel_mean_ADJ .* S_mean) ./ Sref) .* area_metres2) * -1;
%FF_seas = (((vel .* S) ./ Sref) .* area_metres2) * -1;
FF_mean = (((vel_mean_ADJ .* (S_mean-Sref)) / Sref) .* area_metres2) .* -1;

FF_seas = (((vel_seas_ADJ .* (S-Sref)) / Sref) .* area_metres2_3D) .* -1;
%FF_anom = FF_seas - FF_mean;

%% Bin into density space...
% recalculate density for mean T and S
% sigma0 = gsw_sigma0(SA,CT)
dens = gsw_sigma0(S,T);
dens_mean = gsw_sigma0(S_mean,T_mean);
dens_grid = 25.5:0.04:28; % 26.2:0.04:28;

% Useful to have T*area and S*area for later calcs
Tarea = T_mean.*area_metres2;
Sarea = S_mean.*area_metres2;
Tarea_seas = T.*area_metres2_3D;
Sarea_seas = S.*area_metres2_3D;

% preallocate
tran_mean_ADJ_DENS = nan*ones(length(dens_grid),length(dist_grid_vel));
HF_mean_DENS = tran_mean_ADJ_DENS;
FF_mean_DENS = tran_mean_ADJ_DENS;
area_DENS = tran_mean_ADJ_DENS;
Tarea_DENS =  tran_mean_ADJ_DENS;
Sarea_DENS =  tran_mean_ADJ_DENS;
tran_mean_ADJ_BOUND_DENS = tran_mean_ADJ_DENS; 
tran_mean_ADJ_47N_DENS = tran_mean_ADJ_DENS;

% and preallocate seasonal...
tran_seas_ADJ_DENS = nan*ones(length(dens_grid),length(dist_grid_vel),4);
HF_seas_DENS = tran_seas_ADJ_DENS;
FF_seas_DENS = tran_seas_ADJ_DENS;
area_seas_DENS = tran_seas_ADJ_DENS;
Tarea_seas_DENS =  tran_seas_ADJ_DENS;
Sarea_seas_DENS =  tran_seas_ADJ_DENS;
tran_seas_ADJ_BOUND_DENS =  tran_seas_ADJ_DENS;
tran_seas_ADJ_47N_DENS = tran_seas_ADJ_DENS;

% Want a boundary only and 47N section only transport version...
tran_mean_ADJ_BOUND = tran_mean_ADJ; tran_mean_ADJ_BOUND(:,63:end) = nan; 
tran_mean_ADJ_47N = tran_mean_ADJ; tran_mean_ADJ_47N(:,1:62) = nan;
tran_seas_ADJ_BOUND = tran_seas_ADJ; tran_seas_ADJ_BOUND(:,63:end,:) = nan; 
tran_seas_ADJ_47N = tran_seas_ADJ; tran_seas_ADJ_47N(:,1:62,:) = nan;

for aa = 1:length(dist_grid_vel) % for each profile
    for bb = 1:length(dens_grid)-1 % for each density bin
        ind = find(dens_mean(:,aa) >= dens_grid(bb) & dens_mean(:,aa) < dens_grid(bb+1));
        
        % for annual mean properties
        if ~isempty(ind)
            tran_mean_ADJ_DENS(bb,aa) = sum(tran_mean_ADJ(ind,aa));
            tran_mean_ADJ_BOUND_DENS(bb,aa) = sum(tran_mean_ADJ_BOUND(ind,aa));
            tran_mean_ADJ_47N_DENS(bb,aa) = sum(tran_mean_ADJ_47N(ind,aa));
            HF_mean_DENS(bb,aa) = sum(HF_mean(ind,aa));
            FF_mean_DENS(bb,aa) = sum(FF_mean(ind,aa));
            
            % Also compute cross-sectional area of density space at this
            % stage (thanks Neil)
            area_DENS(bb,aa) = sum(area_metres2(ind,aa));
            
            % Map across T*area and S* area into density space
            Tarea_DENS(bb,aa) = sum(Tarea(ind,aa));
            Sarea_DENS(bb,aa) = sum(Sarea(ind,aa));
        end
        
        % Do the same thing, but for individual seasons
        for cc = 1:4
            ind = find(dens(:,aa,cc) >= dens_grid(bb) & dens(:,aa,cc) < dens_grid(bb+1));
            if ~isempty(ind)
                tran_seas_ADJ_DENS(bb,aa,cc) = sum(tran_seas_ADJ(ind,aa,cc));
                tran_seas_ADJ_BOUND_DENS(bb,aa,cc) = sum(tran_seas_ADJ_BOUND(ind,aa,cc));
                tran_seas_ADJ_47N_DENS(bb,aa,cc) = sum(tran_seas_ADJ_47N(ind,aa,cc));
                HF_seas_DENS(bb,aa,cc) = sum(HF_seas(ind,aa,cc));
                FF_seas_DENS(bb,aa,cc) = sum(FF_seas(ind,aa,cc));
                
                % Also compute cross-sectional area of density space at this
                % stage (thanks Neil)
                area_seas_DENS(bb,aa,cc) = sum(area_metres2_3D(ind,aa,cc));
                
                % Map across T*area and S* area into density space
                Tarea_seas_DENS(bb,aa,cc) = sum(Tarea_seas(ind,aa,cc));
                Sarea_seas_DENS(bb,aa,cc) = sum(Sarea_seas(ind,aa,cc));
            end
        end
    end
end

%% Overturning streamfunction

tran_sigma = nansum(tran_mean_ADJ_DENS,2);
OTSF = cumsum(tran_sigma,'omitnan');
[OVERTURNING,sigmaind] = max(OTSF);
iso_interest = dens_grid(sigmaind); % isopycnal of interest.  27.28 on > 1000 m data only. About 9 Sv.
% plot(OTSF,-dens_grid);

% boundary bit only
tran_sigma_BOUND = nansum(tran_mean_ADJ_BOUND_DENS,2);
OTSF_BOUND = cumsum(tran_sigma_BOUND,'omitnan');

% 47N bit only
tran_sigma_47N = nansum(tran_mean_ADJ_47N_DENS,2);
OTSF_47N = cumsum(tran_sigma_47N,'omitnan');




% seasonal version.
tran_sigma_seas = squeeze(nansum(tran_seas_ADJ_DENS,2));
OTSF_seas = cumsum(tran_sigma_seas,1,'omitnan');
[OVERTURNING_seas,sigmaind_seas] = max(OTSF_seas);
iso_interest_seas = dens_grid(sigmaind_seas); % isopycnal of interest.  27.28 on > 1000 m data only. About 9 Sv.

% boundary bit only
tran_sigma_BOUND_seas = squeeze(nansum(tran_seas_ADJ_BOUND_DENS,2));
OTSF_BOUND_seas = cumsum(tran_sigma_BOUND_seas,1,'omitnan');

% 47N bit only
tran_sigma_47N_seas = squeeze(nansum(tran_seas_ADJ_47N_DENS,2));
OTSF_47N_seas = cumsum(tran_sigma_47N_seas,1,'omitnan');


%% Decompose into overturning and along isopycnal components of heat flux and FW flux
% From Neil email 240821

% in density space, find the horizontally-averaged velocity and temperature fields 
% ("horizontally" really meaning "along isopycnals"). 

v_DENS = (tran_mean_ADJ_DENS)./area_DENS ;
v_seas_DENS = (tran_seas_ADJ_DENS)./area_seas_DENS ;

% To get the temperature field in density space, you can divide the heat 
% transport for each cell by rho*Cp*(volume transport for that cell). 

T_DENS = HF_mean_DENS ./ (rhoCp.*(tran_mean_ADJ_DENS));
T_seas_DENS = HF_seas_DENS ./ (rhoCp.*(tran_seas_ADJ_DENS));

S_DENS = Sarea_DENS ./ (area_DENS);
S_seas_DENS = Sarea_seas_DENS ./ (area_seas_DENS);

% Quick sanity plots
% Temp
figure; 
subplot(4,1,1); pcolor(v_seas_DENS(:,:,1)); shading flat; caxis([-0.1 0.1]);
subplot(4,1,2); pcolor(v_seas_DENS(:,:,2)); shading flat; caxis([-0.1 0.1]);
subplot(4,1,3); pcolor(v_seas_DENS(:,:,3)); shading flat; caxis([-0.1 0.1]);
subplot(4,1,4); pcolor(v_seas_DENS(:,:,4)); shading flat; caxis([-0.1 0.1]); colorbar;
title('Seasonal velocity binned into density space (v seas DENS)');
figure
subplot(4,1,1); pcolor(T_seas_DENS(:,:,1)); shading flat; caxis([2 14]);
subplot(4,1,2); pcolor(T_seas_DENS(:,:,2)); shading flat; caxis([2 14]);
subplot(4,1,3); pcolor(T_seas_DENS(:,:,3)); shading flat; caxis([2 14]);
subplot(4,1,4); pcolor(T_seas_DENS(:,:,4)); shading flat; caxis([2 14]); colorbar;
title('Seasonal Temperature binned into density space (v seas DENS)');
% Salinity
figure
subplot(4,1,1); pcolor(S_seas_DENS(:,:,1)); shading flat; caxis([34.5 35.6]);
subplot(4,1,2); pcolor(S_seas_DENS(:,:,2)); shading flat; caxis([34.5 35.6]);
subplot(4,1,3); pcolor(S_seas_DENS(:,:,3)); shading flat; caxis([34.5 35.6]);
subplot(4,1,4); pcolor(S_seas_DENS(:,:,4)); shading flat; caxis([34.5 35.6]); colorbar;
title('Seasonal salinity binned into density space (v seas DENS)');





% mean
area_DENS_bar = nansum(area_DENS,2);
Tarea_DENS_bar = nansum(Tarea_DENS,2);
Sarea_DENS_bar = nansum(Sarea_DENS,2);
varea_DENS_bar = nansum(tran_mean_ADJ_DENS,2);
Tbar = repmat((Tarea_DENS_bar./area_DENS_bar),[1,length(dist_grid_vel)]);
Sbar = repmat((Sarea_DENS_bar./area_DENS_bar),[1,length(dist_grid_vel)]);
vbar = repmat((varea_DENS_bar./area_DENS_bar),[1,length(dist_grid_vel)]);

% seasonal
area_seas_DENS_bar = nansum(area_seas_DENS,2);
Tarea_seas_DENS_bar = nansum(Tarea_seas_DENS,2);
Sarea_seas_DENS_bar = nansum(Sarea_seas_DENS,2);
varea_seas_DENS_bar = nansum(tran_seas_ADJ_DENS,2);
Tbar_seas = repmat((Tarea_seas_DENS_bar./area_seas_DENS_bar),[1,length(dist_grid_vel),1]);
Sbar_seas = repmat((Sarea_seas_DENS_bar./area_seas_DENS_bar),[1,length(dist_grid_vel),1]);
vbar_seas = repmat((varea_seas_DENS_bar./area_seas_DENS_bar),[1,length(dist_grid_vel),1]);

% Quick sanity plots
figure; 
subplot(4,1,1); pcolor(vbar_seas(:,:,1)); shading flat; caxis([-0.1 0.1]);
subplot(4,1,2); pcolor(vbar_seas(:,:,2)); shading flat; caxis([-0.1 0.1]);
subplot(4,1,3); pcolor(vbar_seas(:,:,3)); shading flat; caxis([-0.1 0.1]);
subplot(4,1,4); pcolor(vbar_seas(:,:,4)); shading flat; caxis([-0.1 0.1]); colorbar;
title('vbar seas');

figure
subplot(4,1,1); pcolor(Tbar_seas(:,:,1)); shading flat; caxis([2 14]);
subplot(4,1,2); pcolor(Tbar_seas(:,:,2)); shading flat; caxis([2 14]);
subplot(4,1,3); pcolor(Tbar_seas(:,:,3)); shading flat; caxis([2 14]);
subplot(4,1,4); pcolor(Tbar_seas(:,:,4)); shading flat; caxis([2 14]); colorbar;
title('Tbar seas');

figure
subplot(4,1,1); pcolor(Sbar_seas(:,:,1)); shading flat; caxis([34.5 35.6]);
subplot(4,1,2); pcolor(Sbar_seas(:,:,2)); shading flat; caxis([34.5 35.6]);
subplot(4,1,3); pcolor(Sbar_seas(:,:,3)); shading flat; caxis([34.5 35.6]);
subplot(4,1,4); pcolor(Sbar_seas(:,:,4)); shading flat; caxis([34.5 35.6]); colorbar;
title('Sbar seas');


% Then subtract them from the 
% original fields and you'll have 
% T = T^bar + T'
% v = v^bar + v'

% mean
Tprime = T_DENS - Tbar;
Sprime = S_DENS - Sbar;
vprime = v_DENS - vbar;
% seasonal
Tprime_seas = T_seas_DENS - Tbar_seas;
Sprime_seas = S_seas_DENS - Sbar_seas;
vprime_seas = v_seas_DENS - vbar_seas;

% Test: a = nanmean(Tprime,2);

% where the horizontal mean of T' and v' is zero by definition. The fields 
% T^bar and v^bar only contain information about overturning in density space,
% since they are "horizontally" uniform. The overturning heat flux is then 
% computed by integrating rho*Cp*v^bar*T^bar over the area, so it's just 
% heat flux from the "bar" parts of the fields.

% mean
OHF = area_DENS .* (rhoCp.*vbar.*Tbar);
GHF = area_DENS .* (rhoCp.*vprime.*Tprime);
OFF = area_DENS .* -1.*(vbar.*(Sbar-Sref) ./ Sref); % Check:  OFF = area_DENS .* -1.*(vbar.*(Sbar-Sref) ./ Sref);  Should be the same.
GFF = area_DENS .* -1.*(vprime.*(Sprime-Sref) ./ Sref); % Ditto.
CROSS1 = area_DENS .* -1.*(vbar.*(Sprime-Sref) ./ Sref); 
CROSS2 = area_DENS .* -1.*(vprime.*(Sbar-Sref) ./ Sref); 



% seasonal
OHF_seas = area_seas_DENS .* (rhoCp.*vbar_seas.*Tbar_seas);
GHF_seas = area_seas_DENS .* (rhoCp.*vprime_seas.*Tprime_seas);
OFF_seas = area_seas_DENS .* -1.*(vbar_seas.*(Sbar_seas-Sref) ./ Sref); % Check:  OFF = area_DENS .* -1.*(vbar.*(Sbar-Sref) ./ Sref);  Should be the same.
GFF_seas = area_seas_DENS .* -1.*(vprime_seas.*(Sprime_seas-Sref) ./ Sref); % Ditto.



% Check - heat flux computed in density space
% HF_mean_DENS2 = v_DENS .* T_DENS .* rhoCp .* area_DENS;


%% More sanity plots (not saved)

% Heat %%%%%%%%%%%%%%%%%%%%
figure; 
subplot(4,1,1); 
pcolor(dist_grid_vel,pres_grid,HF_mean); shading flat
hold on;
[C,h] = contour(dist_grid_vel,pres_grid,dens_mean,'color','k'); 
clabel(C,h)
%xlabel('Distance along contour');
ylabel('Pressure');
set(gca,'ydir','reverse')
% ylim([0 1050])
grid on;
caxis([-5e12 5e12]);
colorbar;
title('Heat flux into volume - pressure space');

subplot(4,1,2); pcolor(dist_grid_vel,dens_grid,HF_mean_DENS); shading flat
%xlabel('Distance along contour');
ylabel('\sigma0');
set(gca,'ydir','reverse')
grid on;
caxis([-2e13 2e13]);
colorbar;
title('Heat flux into volume - Density space');

subplot(4,1,3); pcolor(dist_grid_vel,dens_grid,OHF); shading flat
%xlabel('Distance along contour');
ylabel('\sigma0');
set(gca,'ydir','reverse')
grid on;
caxis([-4e12 4e12]);
colorbar;
title('Overturning heat flux');

subplot(4,1,4); pcolor(dist_grid_vel,dens_grid,GHF); shading flat
xlabel('Distance along contour');
ylabel('\sigma0');
set(gca,'ydir','reverse')
grid on;
caxis([-2e12 2e12]);
colorbar;
title('Along isopycnal heat flux');


% Freshwater %%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; 
subplot(4,1,1); 
pcolor(dist_grid_vel,pres_grid,FF_mean); shading flat
hold on;
[C,h] = contour(dist_grid_vel,pres_grid,dens_mean,'color','k'); 
clabel(C,h)
%xlabel('Distance along contour');
ylabel('Pressure');
set(gca,'ydir','reverse')
% ylim([0 1050])
grid on;
caxis([-3000 3000]);
colorbar;
title('F flux into volume - pressure space');

subplot(4,1,2); pcolor(dist_grid_vel,dens_grid,FF_mean_DENS); shading flat
%xlabel('Distance along contour');
ylabel('\sigma0');
set(gca,'ydir','reverse')
grid on;
caxis([-3000 3000]);
colorbar;
title('F flux into volume - Density space');

subplot(4,1,3); pcolor(dist_grid_vel,dens_grid,OFF); shading flat
%xlabel('Distance along contour');
ylabel('\sigma0');
set(gca,'ydir','reverse')
grid on;
%caxis([-4e12 4e12]);
colorbar;
title('Overturning F flux');

subplot(4,1,4); pcolor(dist_grid_vel,dens_grid,GFF); shading flat
xlabel('Distance along contour');
ylabel('\sigma0');
set(gca,'ydir','reverse')
grid on;
%caxis([-2e12 2e12]);
colorbar;
title('Along isopycnal F flux');





%% MAIN OUTPUT FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Er = load('D:\Work_computer_sync\OSNAP_postdoc\PAPERS_NEW\N_Atlantic_boundary\matlab\V3_050321\P50_Estimate_errors\datas\error_estimates.mat');
Er = load('D:\Work_computer_sync\OSNAP_postdoc\PAPERS_NEW\N_Atlantic_boundary\matlab\V3_050321\P50_Estimate_errors\datas\error_estimates.mat');

% Also load the fluxes resulting from uncorrected velocities
uncorr = load('fluxes_uncorrected_for_plot');


%% HEAT FLUX PLOT

% Error estimates
% Create cumulative lines so we don't have to nansum every time
HF_CUMUL(:,1) = cumsum(squeeze(nansum(HF_seas_DENS(:,:,1),1)));
HF_CUMUL(:,2) = cumsum(squeeze(nansum(HF_seas_DENS(:,:,2),1)));
HF_CUMUL(:,3) = cumsum(squeeze(nansum(HF_seas_DENS(:,:,3),1)));
HF_CUMUL(:,4) = cumsum(squeeze(nansum(HF_seas_DENS(:,:,4),1)));
HF_CUMUL_MEAN = nanmean(HF_CUMUL,2);

FF_CUMUL(:,1) = cumsum(squeeze(nansum(FF_seas_DENS(:,:,1),1)));
FF_CUMUL(:,2) = cumsum(squeeze(nansum(FF_seas_DENS(:,:,2),1)));
FF_CUMUL(:,3) = cumsum(squeeze(nansum(FF_seas_DENS(:,:,3),1)));
FF_CUMUL(:,4) = cumsum(squeeze(nansum(FF_seas_DENS(:,:,4),1)));
FF_CUMUL_MEAN = nanmean(FF_CUMUL,2);

% Seasonal
upper_JFM = (HF_CUMUL(:,1) ./ 1e15)' + (Er.HEAT_FLUX_error(1,:) ./ 1e15); lower_JFM = (HF_CUMUL(:,1) ./ 1e15)' - (Er.HEAT_FLUX_error(1,:) ./ 1e15);
patchY_JFM = [upper_JFM lower_JFM(end:-1:1)]; patchX_JFM = [dist_grid_vel dist_grid_vel(end:-1:1)];

upper_AMJ = (HF_CUMUL(:,2) ./ 1e15)' + (Er.HEAT_FLUX_error(2,:) ./ 1e15); lower_AMJ = (HF_CUMUL(:,2) ./ 1e15)' - (Er.HEAT_FLUX_error(2,:) ./ 1e15);
patchY_AMJ = [upper_AMJ lower_AMJ(end:-1:1)]; patchX_AMJ = [dist_grid_vel dist_grid_vel(end:-1:1)];

upper_JAS = (HF_CUMUL(:,3) ./ 1e15)' + (Er.HEAT_FLUX_error(3,:) ./ 1e15); lower_JAS = (HF_CUMUL(:,3) ./ 1e15)' - (Er.HEAT_FLUX_error(3,:) ./ 1e15);
patchY_JAS = [upper_JAS lower_JAS(end:-1:1)]; patchX_JAS = [dist_grid_vel dist_grid_vel(end:-1:1)];

upper_OND = (HF_CUMUL(:,4) ./ 1e15)' + (Er.HEAT_FLUX_error(4,:) ./ 1e15); lower_OND = (HF_CUMUL(:,4) ./ 1e15)' - (Er.HEAT_FLUX_error(4,:) ./ 1e15);
patchY_OND = [upper_OND lower_OND(end:-1:1)]; patchX_OND = [dist_grid_vel dist_grid_vel(end:-1:1)];

% Mean
Er.HEAT_FLUX_error_MEAN = nanmean(Er.HEAT_FLUX_error,1);
upper_MEAN = (HF_CUMUL_MEAN ./ 1e15)' + (Er.HEAT_FLUX_error_MEAN ./ 1e15); lower_MEAN = (HF_CUMUL_MEAN ./ 1e15)' - (Er.HEAT_FLUX_error_MEAN ./ 1e15);
patchY_MEAN = [upper_MEAN lower_MEAN(end:-1:1)]; patchX_MEAN = [dist_grid_vel dist_grid_vel(end:-1:1)];

figure(10);
Cl = colororder;
subplot(2,1,1);
hold on;

% plot patches
% plot error patches
patch(patchX_JFM,patchY_JFM,Cl(1,:),'LineStyle','none');
alpha(0.2)
patch(patchX_AMJ,patchY_AMJ,Cl(2,:),'LineStyle','none');
alpha(0.2)
patch(patchX_JAS,patchY_JAS,Cl(3,:),'LineStyle','none');
alpha(0.2)
patch(patchX_OND,patchY_OND,Cl(4,:),'LineStyle','none');
alpha(0.2)

patch(patchX_MEAN,patchY_MEAN,'k','LineStyle','none');
alpha(0.2)



B = plot(dist_grid_vel,cumsum(squeeze(nansum(HF_seas_DENS(:,:,1),1))) ./ 1e15,'color',Cl(1,:),'linewidth',2);
C = plot(dist_grid_vel,cumsum(squeeze(nansum(HF_seas_DENS(:,:,2),1))) ./ 1e15,'color',Cl(2,:),'linewidth',2,'linestyle','--');
D = plot(dist_grid_vel,cumsum(squeeze(nansum(HF_seas_DENS(:,:,3),1))) ./ 1e15,'color',Cl(3,:),'linewidth',2,'linestyle',':');
E = plot(dist_grid_vel,cumsum(squeeze(nansum(HF_seas_DENS(:,:,4),1))) ./ 1e15,'color',Cl(4,:),'linewidth',2,'linestyle','-.');

A = plot(dist_grid_vel,HF_CUMUL_MEAN ./ 1e15,'k','linewidth',3);

for bb = 1:length(OSNAP_crossings)
    line([OSNAP_crossings(bb) OSNAP_crossings(bb)],[-0.5 0.5],'linewidth',2,'linestyle','--','color','k')
end

% plot(dist_grid_vel,cumsum(squeeze(nansum(uncorr.HF_mean(:,:),1))) ./ 1e15,'color','k','linewidth',3,'linestyle','--'); % this is the mean flux w uncorrected velocities

% B = plot(dist_grid_vel,cumsum(squeeze(nansum(OHF_seas(:,:,1),1))) ./ 1e15,'b','linewidth',2,'linestyle','--');
% plot(dist_grid_vel,cumsum(squeeze(nansum(OHF_seas(:,:,2),1))) ./ 1e15,'g','linewidth',2,'linestyle','--');
% plot(dist_grid_vel,cumsum(squeeze(nansum(OHF_seas(:,:,3),1))) ./ 1e15,'r','linewidth',2,'linestyle','--');
% plot(dist_grid_vel,cumsum(squeeze(nansum(OHF_seas(:,:,4),1))) ./ 1e15,'y','linewidth',2,'linestyle','--');
% 
% C = plot(dist_grid_vel,cumsum(squeeze(nansum(GHF_seas(:,:,1),1))) ./ 1e15,'b','linewidth',2,'linestyle',':');
% plot(dist_grid_vel,cumsum(squeeze(nansum(GHF_seas(:,:,2),1))) ./ 1e15,'g','linewidth',2,'linestyle',':');
% plot(dist_grid_vel,cumsum(squeeze(nansum(GHF_seas(:,:,3),1))) ./ 1e15,'r','linewidth',2,'linestyle',':');
% plot(dist_grid_vel,cumsum(squeeze(nansum(GHF_seas(:,:,4),1))) ./ 1e15,'y','linewidth',2,'linestyle',':');

for bb = 1:length(distance_ticks)
   line([distance_ticks(bb) distance_ticks(bb)],[0.48 0.5],'linewidth',3,'color','k'); 
end

xlabel('Distance along contour (km)','fontsize',16)
set(gca,'fontsize',16)
ylim([-0.5 0.5]);
xlim([0 max(dist_grid_vel)]);
grid on
set(gca,'layer','top')
set(gca,'linewidth',1.5)
yticks([-0.5:0.25:0.5]);
% legend([A,B,C],'Total heat flux','Overturning component','Along isopycnal component','Location','SouthWest');
xlabel('Distance along contour (km)','fontsize',16);
ylabel('Heat flux (PW)','fontsize',16)
box on


%% FRESHWATER FLUX PLOT

subplot(2,1,2);
hold on;

% Error estimates
% Seasonal
upper_JFM = (FF_CUMUL(:,1) ./ 1e6)' + (Er.FRESHWATER_FLUX_error(1,:) ./ 1e6); lower_JFM = (FF_CUMUL(:,1) ./ 1e6)' - (Er.FRESHWATER_FLUX_error(1,:) ./ 1e6);
patchY_JFM = [upper_JFM lower_JFM(end:-1:1)]; patchX_JFM = [dist_grid_vel dist_grid_vel(end:-1:1)];

upper_AMJ = (FF_CUMUL(:,2) ./ 1e6)' + (Er.FRESHWATER_FLUX_error(2,:) ./ 1e6); lower_AMJ = (FF_CUMUL(:,2) ./ 1e6)' - (Er.FRESHWATER_FLUX_error(2,:) ./ 1e6);
patchY_AMJ = [upper_AMJ lower_AMJ(end:-1:1)]; patchX_AMJ = [dist_grid_vel dist_grid_vel(end:-1:1)];

upper_JAS = (FF_CUMUL(:,3) ./ 1e6)' + (Er.FRESHWATER_FLUX_error(3,:) ./ 1e6); lower_JAS = (FF_CUMUL(:,3) ./ 1e6)' - (Er.FRESHWATER_FLUX_error(3,:) ./ 1e6);
patchY_JAS = [upper_JAS lower_JAS(end:-1:1)]; patchX_JAS = [dist_grid_vel dist_grid_vel(end:-1:1)];

upper_OND = (FF_CUMUL(:,4) ./ 1e6)' + (Er.FRESHWATER_FLUX_error(4,:) ./ 1e6); lower_OND = (FF_CUMUL(:,4) ./ 1e6)' - (Er.FRESHWATER_FLUX_error(4,:) ./ 1e6);
patchY_OND = [upper_OND lower_OND(end:-1:1)]; patchX_OND = [dist_grid_vel dist_grid_vel(end:-1:1)];

% Mean
Er.FRESHWATER_FLUX_error_MEAN = nanmean(Er.FRESHWATER_FLUX_error,1);
upper_MEAN = (FF_CUMUL_MEAN ./ 1e6)' + (Er.FRESHWATER_FLUX_error_MEAN ./ 1e6); lower_MEAN = (FF_CUMUL_MEAN ./ 1e6)' - (Er.FRESHWATER_FLUX_error_MEAN ./ 1e6);
patchY_MEAN = [upper_MEAN lower_MEAN(end:-1:1)]; patchX_MEAN = [dist_grid_vel dist_grid_vel(end:-1:1)];


% plot patches
% plot error patches
patch(patchX_JFM,patchY_JFM,Cl(1,:),'LineStyle','none');
alpha(0.2)
patch(patchX_AMJ,patchY_AMJ,Cl(2,:),'LineStyle','none');
alpha(0.2)
patch(patchX_JAS,patchY_JAS,Cl(3,:),'LineStyle','none');
alpha(0.2)
patch(patchX_OND,patchY_OND,Cl(4,:),'LineStyle','none');
alpha(0.2)

patch(patchX_MEAN,patchY_MEAN,'k','LineStyle','none');
alpha(0.2)

for bb = 1:length(OSNAP_crossings)
    line([OSNAP_crossings(bb) OSNAP_crossings(bb)],[-0.3 0.2],'linewidth',2,'linestyle','--','color','k')
end



B = plot(dist_grid_vel,cumsum(squeeze(nansum(FF_seas(:,:,1),1)) ./ 1e6),'color',Cl(1,:),'linewidth',2);
C = plot(dist_grid_vel,cumsum(squeeze(nansum(FF_seas(:,:,2),1))) ./ 1e6,'color',Cl(2,:),'linewidth',2,'linestyle','--');
D = plot(dist_grid_vel,cumsum(squeeze(nansum(FF_seas(:,:,3),1))) ./ 1e6,'color',Cl(3,:),'linewidth',2,'linestyle',':');
E = plot(dist_grid_vel,cumsum(squeeze(nansum(FF_seas(:,:,4),1))) ./ 1e6,'color',Cl(4,:),'linewidth',2,'linestyle','-.');

A = plot(dist_grid_vel,FF_CUMUL_MEAN ./ 1e6,'k','linewidth',3);

% plot(dist_grid_vel,cumsum(squeeze(nansum(uncorr.FF_mean(:,:),1))) ./ 1e6,'color','k','linewidth',3,'linestyle','--'); % this is the mean flux w uncorrected velocities

% B = plot(dist_grid_vel,cumsum(squeeze(nansum(OFF_seas(:,:,1),1))) ./ 1e6,'b','linewidth',2,'linestyle','--');
% plot(dist_grid_vel,cumsum(squeeze(nansum(OFF_seas(:,:,2),1))) ./ 1e6,'g','linewidth',2,'linestyle','--');
% plot(dist_grid_vel,cumsum(squeeze(nansum(OFF_seas(:,:,3),1))) ./ 1e6,'r','linewidth',2,'linestyle','--');
% plot(dist_grid_vel,cumsum(squeeze(nansum(OFF_seas(:,:,4),1))) ./ 1e6,'y','linewidth',2,'linestyle','--');
% 
% C = plot(dist_grid_vel,cumsum(squeeze(nansum(GFF_seas(:,:,1),1))) ./ 1e6,'b','linewidth',2,'linestyle',':');
% plot(dist_grid_vel,cumsum(squeeze(nansum(GFF_seas(:,:,2),1))) ./ 1e6,'g','linewidth',2,'linestyle',':');
% plot(dist_grid_vel,cumsum(squeeze(nansum(GFF_seas(:,:,3),1))) ./ 1e6,'r','linewidth',2,'linestyle',':');
% plot(dist_grid_vel,cumsum(squeeze(nansum(GFF_seas(:,:,4),1))) ./ 1e6,'y','linewidth',2,'linestyle',':');

% for bb = 1:length(distance_ticks)
%    line([distance_ticks(bb) distance_ticks(bb)],[0.78 0.8],'linewidth',3,'color','k'); 
% end

xlabel('Distance along contour (km)','fontsize',16)
set(gca,'fontsize',16)
% ylim([-1 0.8]);
xlim([0 max(dist_grid_vel)]);
grid on
set(gca,'layer','top')
set(gca,'linewidth',1.5)
% legend([A,B,C],'Total freshwater flux','Overturning component','Along isopycnal component','Location','SouthWest');
xlabel('Distance along contour (km)','fontsize',16);
ylabel('Freshwater flux (Sv)','fontsize',16)
legend([A,B,C,D,E],'Annual mean','JFM','AMJ','JAS','OND','Location','SouthWest');
box on

%% print figure
width  = 2000;  % frame width
height = 1400;  % frame height
 pngname = ('plots/p17a_cumulative_fluxes.png');
%pngname = ('plots/p17a_cumulative_fluxes_EFFECT_OF_CORRECTION.png');

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












%% Overturning stream function plot
% Er = load('D:\Work_computer_sync\OSNAP_postdoc\PAPERS_NEW\N_Atlantic_boundary\matlab\V3_050321\P50_Estimate_errors\datas\error_estimates.mat');
figure(20);
subplot(1,3,1);
hold on; box on;
%plot(nanmean(OTSF_seas,2),dens_grid,'m','linewidth',1); % mean of seasonal (test)
% seasonal

%% Create coloured patches for error regions
% Mean
upper_MEAN = (OTSF ./ 1e6)' + Er.OTSF_MEAN_error; lower_MEAN = (OTSF ./ 1e6)' - Er.OTSF_MEAN_error;
patchY_MEAN = [upper_MEAN lower_MEAN(end:-1:1)]; patchX_MEAN = [dens_grid dens_grid(end:-1:1)];

% Seasonal
upper_JFM = (OTSF_seas(:,1) ./ 1e6)' + Er.OTSF_JFM_error; lower_JFM = (OTSF_seas(:,1) ./ 1e6)' - Er.OTSF_JFM_error;
patchY_JFM = [upper_JFM lower_JFM(end:-1:1)]; patchX_JFM = [dens_grid dens_grid(end:-1:1)];

upper_AMJ = (OTSF_seas(:,2) ./ 1e6)' + Er.OTSF_AMJ_error; lower_AMJ = (OTSF_seas(:,2) ./ 1e6)' - Er.OTSF_AMJ_error;
patchY_AMJ = [upper_AMJ lower_AMJ(end:-1:1)]; patchX_AMJ = [dens_grid dens_grid(end:-1:1)];

upper_JAS = (OTSF_seas(:,3) ./ 1e6)' + Er.OTSF_JAS_error; lower_JAS = (OTSF_seas(:,3) ./ 1e6)' - Er.OTSF_JAS_error;
patchY_JAS = [upper_JAS lower_JAS(end:-1:1)]; patchX_JAS = [dens_grid dens_grid(end:-1:1)];

upper_OND = (OTSF_seas(:,4) ./ 1e6)' + Er.OTSF_OND_error; lower_OND = (OTSF_seas(:,4) ./ 1e6)' - Er.OTSF_OND_error;
patchY_OND = [upper_OND lower_OND(end:-1:1)]; patchX_OND = [dens_grid dens_grid(end:-1:1)];

% plot error patches
patch(patchY_JFM,patchX_JFM,Cl(1,:),'LineStyle','none');
alpha(0.2)
patch(patchY_AMJ,patchX_AMJ,Cl(2,:),'LineStyle','none');
alpha(0.2)
patch(patchY_JAS,patchX_JAS,Cl(3,:),'LineStyle','none');
alpha(0.2)
patch(patchY_OND,patchX_OND,Cl(4,:),'LineStyle','none');
alpha(0.2)
patch(patchY_MEAN,patchX_MEAN,'k','LineStyle','none');
alpha(0.2)



B = plot(OTSF_seas(:,1) ./ 1e6,dens_grid,'color',Cl(1,:),'linewidth',1.5); 
C = plot(OTSF_seas(:,2) ./ 1e6,dens_grid,'color',Cl(2,:),'linewidth',1.5,'linestyle','--'); 
D = plot(OTSF_seas(:,3) ./ 1e6,dens_grid,'color',Cl(3,:),'linewidth',1.5,'linestyle',':'); 
E = plot(OTSF_seas(:,4) ./ 1e6,dens_grid,'color',Cl(4,:),'linewidth',1.5,'linestyle','-.'); 

A = plot(OTSF ./ 1e6,dens_grid,'k','linewidth',2); % mean

% Plot first maxima
[mx,ind] = max(OTSF_seas(:,1)); plot(mx ./ 1e6,dens_grid(ind),'o','color',Cl(1,:),'linewidth',2);
[mx,ind] = max(OTSF_seas(:,2)); plot(mx ./ 1e6,dens_grid(ind),'o','color',Cl(2,:),'linewidth',2);
[mx,ind] = max(OTSF_seas(:,3)); plot(mx ./ 1e6,dens_grid(ind),'o','color',Cl(3,:),'linewidth',2);
[mx,ind] = max(OTSF_seas(:,4)); plot(mx ./ 1e6,dens_grid(ind),'o','color',Cl(4,:),'linewidth',2);
[mx,ind] = max(OTSF); plot(mx ./ 1e6,dens_grid(ind),'ok','linewidth',3);

% Plot second maxima
% [mx,ind] = max(OTSF_seas(51:end,1)); plot(mx ./ 1e6,dens_grid(50+ind),'ob','linewidth',2);
[mx,ind] = max(OTSF_seas(51:end,2)); plot(mx ./ 1e6,dens_grid(50+ind),'o','color',Cl(2,:),'linewidth',2);
[mx,ind] = max(OTSF_seas(51:end,3)); plot(mx ./ 1e6,dens_grid(50+ind),'o','color',Cl(3,:),'linewidth',2);
[mx,ind] = max(OTSF_seas(51:end,4)); plot(mx ./ 1e6,dens_grid(50+ind),'o','color',Cl(4,:),'linewidth',2);
[mx,ind] = max(OTSF(51:end)); plot(mx ./ 1e6,dens_grid(50+ind),'ok','linewidth',3);



% Plot densities of interest as horizontal lines
% line([-10 15],[26.86 26.86],'linewidth',1.5,'linestyle','--','color',[0.5 0.5 0.5]); 
% text(11,26.86-0.04,'26.86','fontsize',16,'color',[0.5 0.5 0.5]);

line([-10 15],[27.30 27.30],'linewidth',1.5,'linestyle','--','color',[0.5 0.5 0.5]); 
text(-3,27.30-0.04,'27.30','fontsize',16,'color',[0.5 0.5 0.5]);

line([-10 15],[27.42 27.42],'linewidth',1.5,'linestyle','--','color',[0.5 0.5 0.5]); 
text(10,27.42-0.04,'27.42','fontsize',16,'color',[0.5 0.5 0.5]);

line([-10 15],[27.54 27.54],'linewidth',1.5,'linestyle','--','color',[0.5 0.5 0.5]); 
text(10,27.54-0.04,'27.54','fontsize',16,'color',[0.5 0.5 0.5]);

line([-10 15],[27.7 27.7],'linewidth',1.5,'linestyle','--','color',[0.5 0.5 0.5]); 
text(10,27.7-0.04,'27.7','fontsize',16,'color',[0.5 0.5 0.5]);

% Hatched lower region
% Apply Hatch Fill
patchX = [-10 15 15 -10];
patchY = [27.7 27.7 28 28];
h1 = patch(patchX,patchY,[0.5 0.5 0.5]);
set([h1], 'LineStyle', 'none')
hh1 = hatchfill(h1, 'single', -45, 3);
set(hh1,'Color',[0.4 0.4 0.4]);

set(gca,'YDir','reverse');
set(gca,'fontsize',18);
xlabel('\Psi (Sv)','fontsize',18);
ylabel('\sigma_0 (kg m^-^3)','fontsize',18);
grid on
ylim([26.5 28]);
yticks([26.6:0.2:28]);
xlim([-10 15]);
set(gca,'layer','top')
set(gca,'linewidth',1.5)

%legend([A,B,C,D,E],'Annual mean','JFM','AMJ','JAS','OND','Location','NorthWest');






subplot(1,3,2);
hold on; box on;


%% Create coloured patches for error regions
% Mean
upper_BOUND = (OTSF_BOUND ./ 1e6)' + Er.OTSF_BOUND_error; lower_BOUND = (OTSF_BOUND ./ 1e6)' - Er.OTSF_BOUND_error;
patchY_BOUND = [upper_BOUND lower_BOUND(end:-1:1)]; patchX_BOUND = [dens_grid dens_grid(end:-1:1)];

% Seasonal
upper_JFM = (OTSF_BOUND_seas(:,1) ./ 1e6)' + Er.OTSF_BOUND_JFM_error; lower_JFM = (OTSF_BOUND_seas(:,1) ./ 1e6)' - Er.OTSF_BOUND_JFM_error;
patchY_JFM = [upper_JFM lower_JFM(end:-1:1)]; patchX_JFM = [dens_grid dens_grid(end:-1:1)];

upper_AMJ = (OTSF_BOUND_seas(:,2) ./ 1e6)' + Er.OTSF_BOUND_AMJ_error; lower_AMJ = (OTSF_BOUND_seas(:,2) ./ 1e6)' - Er.OTSF_BOUND_AMJ_error;
patchY_AMJ = [upper_AMJ lower_AMJ(end:-1:1)]; patchX_AMJ = [dens_grid dens_grid(end:-1:1)];

upper_JAS = (OTSF_BOUND_seas(:,3) ./ 1e6)' + Er.OTSF_BOUND_JAS_error; lower_JAS = (OTSF_BOUND_seas(:,3) ./ 1e6)' - Er.OTSF_BOUND_JAS_error;
patchY_JAS = [upper_JAS lower_JAS(end:-1:1)]; patchX_JAS = [dens_grid dens_grid(end:-1:1)];

upper_OND = (OTSF_BOUND_seas(:,4) ./ 1e6)' + Er.OTSF_BOUND_OND_error; lower_OND = (OTSF_BOUND_seas(:,4) ./ 1e6)' - Er.OTSF_BOUND_OND_error;
patchY_OND = [upper_OND lower_OND(end:-1:1)]; patchX_OND = [dens_grid dens_grid(end:-1:1)];

% plot error patches
patch(patchY_JFM,patchX_JFM,Cl(1,:),'LineStyle','none');
alpha(0.2)
patch(patchY_AMJ,patchX_AMJ,Cl(2,:),'LineStyle','none');
alpha(0.2)
patch(patchY_JAS,patchX_JAS,Cl(3,:),'LineStyle','none');
alpha(0.2)
patch(patchY_OND,patchX_OND,Cl(4,:),'LineStyle','none');
alpha(0.2)
patch(patchY_BOUND,patchX_BOUND,'k','LineStyle','none');
alpha(0.2)



B = plot(OTSF_BOUND_seas(:,1) ./ 1e6,dens_grid,'color',Cl(1,:),'linewidth',1.5); 
C = plot(OTSF_BOUND_seas(:,2) ./ 1e6,dens_grid,'color',Cl(2,:),'linewidth',1.5,'linestyle','--'); 
D = plot(OTSF_BOUND_seas(:,3) ./ 1e6,dens_grid,'color',Cl(3,:),'linewidth',1.5,'linestyle',':'); 
E = plot(OTSF_BOUND_seas(:,4) ./ 1e6,dens_grid,'color',Cl(4,:),'linewidth',1.5,'linestyle','-.'); 

A = plot(OTSF_BOUND ./ 1e6,dens_grid,'k','linewidth',2); % mean



% Densities of interest
line([-28 0],[27.30 27.30],'linewidth',1.5,'linestyle','--','color',[0.5 0.5 0.5]); 
text(-22,27.30-0.04,'27.30','fontsize',16,'color',[0.5 0.5 0.5]);

line([-28 0],[27.42 27.42],'linewidth',1.5,'linestyle','--','color',[0.5 0.5 0.5]); 
text(-8,27.42-0.04,'27.42','fontsize',16,'color',[0.5 0.5 0.5]);

line([-28 0],[27.54 27.54],'linewidth',1.5,'linestyle','--','color',[0.5 0.5 0.5]); 
text(-8,27.54-0.04,'27.54','fontsize',16,'color',[0.5 0.5 0.5]);

line([-28 0],[27.7 27.7],'linewidth',1.5,'linestyle','--','color',[0.5 0.5 0.5]); 
text(-8,27.7-0.04,'27.7','fontsize',16,'color',[0.5 0.5 0.5]);



set(gca,'YDir','reverse');
set(gca,'fontsize',18);
xlabel('\Psi (Sv)','fontsize',18);
% ylabel('\sigma_0 (kg m^-^3)','fontsize',18);
grid on
ylim([26.5 28]);
yticks([26.6:0.2:28]);
%xlim([-10 15]);
set(gca,'layer','top')
set(gca,'linewidth',1.5)





subplot(1,3,3);
hold on; box on;


%% Create coloured patches for error regions
% Mean
upper_47N = (OTSF_47N ./ 1e6)' + Er.OTSF_47N_error; lower_47N = (OTSF_47N ./ 1e6)' - Er.OTSF_47N_error;
patchY_47N = [upper_47N lower_47N(end:-1:1)]; patchX_47N = [dens_grid dens_grid(end:-1:1)];

% Seasonal
upper_JFM = (OTSF_47N_seas(:,1) ./ 1e6)' + Er.OTSF_47N_JFM_error; lower_JFM = (OTSF_47N_seas(:,1) ./ 1e6)' - Er.OTSF_47N_JFM_error;
patchY_JFM = [upper_JFM lower_JFM(end:-1:1)]; patchX_JFM = [dens_grid dens_grid(end:-1:1)];

upper_AMJ = (OTSF_47N_seas(:,2) ./ 1e6)' + Er.OTSF_47N_AMJ_error; lower_AMJ = (OTSF_47N_seas(:,2) ./ 1e6)' - Er.OTSF_47N_AMJ_error;
patchY_AMJ = [upper_AMJ lower_AMJ(end:-1:1)]; patchX_AMJ = [dens_grid dens_grid(end:-1:1)];

upper_JAS = (OTSF_47N_seas(:,3) ./ 1e6)' + Er.OTSF_47N_JAS_error; lower_JAS = (OTSF_47N_seas(:,3) ./ 1e6)' - Er.OTSF_47N_JAS_error;
patchY_JAS = [upper_JAS lower_JAS(end:-1:1)]; patchX_JAS = [dens_grid dens_grid(end:-1:1)];

upper_OND = (OTSF_47N_seas(:,4) ./ 1e6)' + Er.OTSF_47N_OND_error; lower_OND = (OTSF_47N_seas(:,4) ./ 1e6)' - Er.OTSF_47N_OND_error;
patchY_OND = [upper_OND lower_OND(end:-1:1)]; patchX_OND = [dens_grid dens_grid(end:-1:1)];

% plot error patches
patch(patchY_JFM,patchX_JFM,Cl(1,:),'LineStyle','none');
alpha(0.2)
patch(patchY_AMJ,patchX_AMJ,Cl(2,:),'LineStyle','none');
alpha(0.2)
patch(patchY_JAS,patchX_JAS,Cl(3,:),'LineStyle','none');
alpha(0.2)
patch(patchY_OND,patchX_OND,Cl(4,:),'LineStyle','none');
alpha(0.2)
patch(patchY_47N,patchX_47N,'k','LineStyle','none');
alpha(0.2)



B = plot(OTSF_47N_seas(:,1) ./ 1e6,dens_grid,'color',Cl(1,:),'linewidth',1.5); 
C = plot(OTSF_47N_seas(:,2) ./ 1e6,dens_grid,'color',Cl(2,:),'linewidth',1.5,'linestyle','--'); 
D = plot(OTSF_47N_seas(:,3) ./ 1e6,dens_grid,'color',Cl(3,:),'linewidth',1.5,'linestyle',':'); 
E = plot(OTSF_47N_seas(:,4) ./ 1e6,dens_grid,'color',Cl(4,:),'linewidth',1.5,'linestyle','-.'); 

A = plot(OTSF_47N ./ 1e6,dens_grid,'k','linewidth',2); % mean



% Plot second maxima
[mx,ind] = max(OTSF_47N_seas(1:55,1)); plot(mx ./ 1e6,dens_grid(ind),'o','color',Cl(1,:),'linewidth',2);
[mx,ind] = max(OTSF_47N_seas(1:55,2)); plot(mx ./ 1e6,dens_grid(ind),'o','color',Cl(2,:),'linewidth',2);
[mx,ind] = max(OTSF_47N_seas(1:55,3)); plot(mx ./ 1e6,dens_grid(ind),'o','color',Cl(3,:),'linewidth',2);
[mx,ind] = max(OTSF_47N_seas(1:55,4)); plot(mx ./ 1e6,dens_grid(ind),'o','color',Cl(4,:),'linewidth',2);
[mx,ind] = max(OTSF_47N(1:55)); plot(mx ./ 1e6,dens_grid(ind),'ok','linewidth',3);


% Densities of interest
line([-2 28],[27.30 27.30],'linewidth',1.5,'linestyle','--','color',[0.5 0.5 0.5]); 
text(0,27.30-0.04,'27.30','fontsize',16,'color',[0.5 0.5 0.5]);

line([-2 28],[27.42 27.42],'linewidth',1.5,'linestyle','--','color',[0.5 0.5 0.5]); 
text(0,27.42-0.04,'27.42','fontsize',16,'color',[0.5 0.5 0.5]);

line([-2 28],[27.54 27.54],'linewidth',1.5,'linestyle','--','color',[0.5 0.5 0.5]); 
text(0,27.54-0.04,'27.54','fontsize',16,'color',[0.5 0.5 0.5]);

line([-2 28],[27.7 27.7],'linewidth',1.5,'linestyle','--','color',[0.5 0.5 0.5]); 
text(0,27.7-0.04,'27.7','fontsize',16,'color',[0.5 0.5 0.5]);



% Hatched lower region
% Apply Hatch Fill
patchX = [-2 28 28 -2];
patchY = [27.7 27.7 28 28];
h1 = patch(patchX,patchY,[0.5 0.5 0.5]);
set([h1], 'LineStyle', 'none')
hh1 = hatchfill(h1, 'single', -45, 3);
set(hh1,'Color',[0.4 0.4 0.4]);


set(gca,'YDir','reverse');
set(gca,'fontsize',18);
xlabel('\Psi (Sv)','fontsize',18);
% ylabel('\sigma_0 (kg m^-^3)','fontsize',18);
grid on
ylim([26.5 28]);
yticks([26.6:0.2:28]);
% xlim([-10 15]);
set(gca,'layer','top')
set(gca,'linewidth',1.5)

legend([A,B,C,D,E],'Annual mean','JFM','AMJ','JAS','OND','Location','NorthEast');







%% print figure
width  = 2200;  % frame width
height = 1200;  % frame height
pngname = ('plots/p17b_OTSF_3parts.png');

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



%% Reporting double integral fluxes

disp('***** ANNUAL MEAN HEAT FLUX INTEGRALS ********');
disp(['Total HF: ' num2str(nansum(nansum(HF_mean_DENS)) / 1e15 ) ' PW']);
disp(['Overturning HF: ' num2str(nansum(nansum(OHF))  / 1e15 ) ' PW']);
disp(['Gyre HF: ' num2str(nansum(nansum(GHF))  / 1e15 ) ' PW']);
disp(['Overturning + gyre: ' num2str((nansum(nansum(OHF)) + nansum(nansum(GHF)))  / 1e15 ) ' PW']);
disp(['overturning / total: ' num2str((nansum(nansum(OHF)) ./ nansum(nansum(HF_mean_DENS))))]);

disp('********************************************');

disp('***** SEASONAL HEAT FLUX INTEGRALS ********');
disp(['Total HF: ' num2str(nansum(nansum(HF_seas_DENS)) ./ 1e15 ) ' PW']);
disp(['Overturning HF: ' num2str(nansum(nansum(OHF_seas))  ./ 1e15 ) ' PW']);
disp(['Gyre HF: ' num2str(nansum(nansum(GHF_seas))  ./ 1e15 ) ' PW']);
disp(['Overturning + gyre: ' num2str((nansum(nansum(OHF_seas)) + nansum(nansum(GHF_seas)))  ./ 1e15 ) ' PW']);
disp(['overturning / total: ' num2str((nansum(nansum(OHF_seas)) ./ nansum(nansum(HF_seas_DENS))))]);

disp('********************************************');

disp('***** ANNUAL MEAN FRESHWATER FLUX INTEGRALS ********');
disp(['Total FF: ' num2str(nansum(nansum(FF_mean_DENS)) / 1e6 ) ' Sv']);
disp(['Overturning FF: ' num2str(nansum(nansum(OFF))  / 1e6 ) ' Sv']);
disp(['Gyre FF: ' num2str(nansum(nansum(GFF))  / 1e6 ) ' Sv']);
disp(['Overturning + gyre: ' num2str((nansum(nansum(OFF)) + nansum(nansum(GFF)))  / 1e6 ) ' Sv']);
disp(['overturning / total: ' num2str((nansum(nansum(OFF)) ./ nansum(nansum(FF_mean_DENS))))]);

disp('********************************************');

disp('***** SEASONAL FRESHWATER FLUX INTEGRALS ********');
disp(['Total FF: ' num2str(nansum(nansum(FF_seas_DENS)) ./ 1e6 ) ' Sv']);
disp(['Overturning FF: ' num2str(nansum(nansum(OFF_seas))  ./ 1e6 ) ' Sv']);
disp(['Gyre FF: ' num2str(nansum(nansum(GFF_seas))  ./ 1e6 ) ' Sv']);
disp(['Overturning + gyre: ' num2str((nansum(nansum(OFF_seas)) + nansum(nansum(GFF_seas)))  ./ 1e6 ) ' Sv']);
disp(['overturning / total: ' num2str((nansum(nansum(OFF_seas)) ./ nansum(nansum(FF_seas_DENS))))]);

