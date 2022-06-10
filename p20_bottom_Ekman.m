% p20_bottom_Ekman
% Estimate bottom Ekman from Argo drift speeds
clear; close all;

addpath(genpath('D:\Work_computer_sync\MATLAB_functions'));


% Load 1000m contour data
load('D:\Work_computer_sync\OSNAP_postdoc\PAPERS_NEW\N_Atlantic_boundary\matlab\V3_050321\intermediate_saves/contour_data_1000_EN4inserted.mat');

% Up to end of boundary bit
for aa = 1:63
    lat(aa) = lat_grid(aa);
    dist(aa) = 150;
end

% Greater distance associated with first grid cell
dist(1) = 275;

% Compute coriolis parameter f
% f = 2omega .sind.lat
omega = 7.27e-5;

for aa = 1:length(lat)
f(aa) = 2*omega*sind(lat(aa));
end

%% calculate Ekman drain per unit distance for each grid cell
kb = 0.0025;
for aa = 1:length(lat)
    V(aa) = (kb*0.11^2)/f(aa);
end

%% Integrate over x distance
tran = V.*(dist*1000);
total_tran = sum(tran)
