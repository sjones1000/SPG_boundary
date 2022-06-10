% p1a_load_and_subset_WOD
%
% Initial load and sift of WOD profiles to get the dataset down to a manageable size.
% NOTE contour extraction part removed from this version.
% SJ 310321
%%%%% Takes about 45 mins to run on full data %%%%%%

clear; close all;

addpath(genpath('D:\Work_computer_sync\MATLAB_functions'));
addpath(genpath('D:\Work_computer_sync\OSNAP_postdoc\PAPERS_NEW\N_Atlantic_boundary\datas'));

dist_from_bathy = 75; % km (was 50)
%ocean.pres = (0:10:2000)';
basetime = datenum('1770-01-01'); % <- WOD convention
min_date = datenum('01-01-2000'); % <- the earliest data we are interested in
min_prof_length = 15;




%% Batch file read. Create a list of all .nc files
cd E:\WOD_CTD_N_Atlantic\CTD
%!dir *.nc /o/b >list.dat
% Read list
fid=fopen('list.dat');
[filelist]=textread('list.dat','%s');
fclose(fid);

filename = filelist{1};

lon = double(ncread(filename,'lon'));
lat = double(ncread(filename,'lat'));



%% Load and format bathy
load('GEBCO_world_1D.mat')

minlon = -65;
maxlon = 5;
minlat = 45;
maxlat = 67;

% Subset to region of interest
latind = find(latitude > minlat & latitude < maxlat);
lonind = find(longitude > minlon & longitude < maxlon);
latitude = latitude(latind);
longitude = longitude(lonind);
bathy = bathy(lonind,latind);

% subset further
latitude = latitude(1:2:end);
longitude = longitude(1:2:end);
bathy = bathy(1:2:end,1:2:end);

% Get array size
NX = length(longitude);
NY = length(latitude);

% Prepare bathy for interpolating prof locations onto...
veclen = NX*NY;
[lonm,latm] = meshgrid(longitude,latitude);
Xin = reshape(lonm,veclen,1);
Yin = reshape(latm,veclen,1);
datain = reshape(bathy',veclen,1);

%% PART 1: EXTRACT CONTOUR FROM BATHY % <- removed.

load('D:\Work_computer_sync\OSNAP_postdoc\PAPERS_NEW\N_Atlantic_boundary\matlab\V3_050321\intermediate_saves\contour_data_1000.mat');

% Coarse contour for first pass subset
cont_lon_S = cont_lon(1:50:end);
cont_lat_S = cont_lat(1:50:end);
cont_index = (1:50:length(cont_lat));





lost_due_to_sal_catch = 0;
lost_due_to_pres_catch = 0;
lost_due_to_time_clause = 0;

%% PART 2: LOAD EACH PROFILE. IS IT IN THE VICINITY OF THE CONTOUR? IF NOT, MOVE ON.
dd = 1; nopres = 0;
for aa = 1:length(filelist)
    filename = filelist{aa};
    
    lon = double(ncread(filename,'lon'));
    lat = double(ncread(filename,'lat'));
    
    
    
    
    %% 0.5th PASS FILTER: NO POINT CONTINUING IF DATE IS OUT
    temp_juld = double(ncread(filename,'time')) + basetime;
    
    if temp_juld >= min_date
        
        
        
        
        
        %% Interpolate bathy grid onto profile locations (nearest neighbour)
        lon_delete = abs(lon - longitude); lonminind = find(lon_delete == min(lon_delete)); lonind = lonminind(1); % X
        lat_delete = abs(lat - latitude); latminind = find(lat_delete == min(lat_delete)); latind = latminind(1); % Y
        profdepth = bathy(lonind,latind);
        
        
        
        
        
        
        
        %% FIRST PASS FILTER: IS THE PROFILE IN > (1000) M WATER?
        if profdepth < -1000 % | (lon >= -46.5 & lon <= -43.5 & lat >= 47 & lat <= 47.8) % Note additional clause as we want to include the Flemish cap
            % If so: Coarse comparision between ALL profiles and contour
            %disp('Coarse subset');
            
            lon_in = [repmat(lon,length(cont_lon_S),1) cont_lon_S'];
            lat_in = [repmat(lat,length(cont_lat_S),1) cont_lat_S'];
            
            % distance = gsw_distance(long,lat,{p})
            Cdist = gsw_distance(lon_in,lat_in);
            % locoate where on the contour subset the min distance is.
            min_loc = find(Cdist == min(Cdist));
            
            % now throw out all other distances
            Cdist = min(Cdist);
            % to km
            Cdist = Cdist./1000;
            
            
            
            
            %% SECOND FILTER: IS THE PROFILE APPROXIMATELY NEAR 1800 (1000) M CONTOUR?
            if Cdist < dist_from_bathy*2
                % If so:
                %disp('fine subset');
                cont_vicinity = cont_index(min_loc); % near the right place on the 1800 m contour.
                
                % Need a clause to catch the ends of the contour strip
                if cont_vicinity < 200 % if near the start
                    cont_local_index = 1 : cont_vicinity +200;
                elseif cont_vicinity > (length(cont_lon) - 200) % if near the end
                    cont_local_index = cont_vicinity -200 : length(cont_lon);
                else
                    cont_local_index = cont_vicinity -200 : cont_vicinity +200;
                end
                
                cont_local_lon = cont_lon(cont_local_index);
                cont_local_lat = cont_lat(cont_local_index);
                
                lon_in = [repmat(lon,length(cont_local_lon),1) cont_local_lon'];
                lat_in = [repmat(lat,length(cont_local_lat),1) cont_local_lat'];
                
                % distance = gsw_distance(long,lat,{p})
                Cdist = gsw_distance(lon_in,lat_in);
                % locoate where on the contour subset the min distance is.
                min_loc = find(Cdist == min(Cdist)); min_loc = min_loc(1);
                
                % now throw out all other distances
                Cdist = min(Cdist);
                % to km
                Cdist = Cdist./1000;
                
                
                
                no_salinity = 0;
                %% THIRD FILTER: IS THE PROFILE WITHIN THE TRUE RANGE OF BATHY?
                if Cdist < dist_from_bathy
                    
                    %% 3.5th FILTER: NO POINT IN CONTINUING IF PROFILE IS TOO SHORT OR NO SALINITY
                    
                    try
                        temp_sal = double(ncread(filename,'Salinity'));
                    catch
                        no_salinity = 1;
                    end
                    
                    if no_salinity == 0
                        
                        %% V3: 3.7th FILTER: IS LAT > 47? WANT TO EXCLUDE ALL PROFS ON WRONG SIDE OF 47 TRANSECT.
                        if lat > 47
                            
                            temp_depth = double(ncread(filename,'z'));
                            if length(temp_depth) >= min_prof_length
                                        
                                
                                %% FINALLY >>>>>> If so, populate 'ocean' variable.
                                ocean.lon(1,dd) = lon;
                                ocean.lat(1,dd) = lat;
                                ocean.profdepth(1,dd) = profdepth;
                                
                                % From netcdf
                                ocean.juld(1,dd) = double(ncread(filename,'time')) + basetime;
                                
                                %% These variables need interpolating onto ocean.pres (defined at top of script)
                                raw_temp = double(ncread(filename,'Temperature'));
                                raw_tempF = double(ncread(filename,'Temperature_WODflag'));
                                
                                raw_sal = double(ncread(filename,'Salinity'));
                                raw_salF = double(ncread(filename,'Salinity_WODflag'));
                                
                                raw_depth = double(ncread(filename,'z'));
                                raw_depthF = double(ncread(filename,'z_WODflag'));
                                
                                try
                                    raw_pres = double(ncread(filename,'Pressure'));
                                catch % if no pressure present
                                    raw_pres = gsw_p_from_z(-raw_depth,lat);
                                    disp('***********NO PRES... calculating from depth ***************');
                                end
                                
                                
                                
                                % apply flags
                                raw_temp(raw_tempF ~= 0) = nan;
                                raw_sal(raw_salF ~= 0) = nan;
                                raw_depth(raw_depthF ~= 0) = nan;
                                
                                % also remove random fill values which have
                                % crept in!
                                raw_sal(raw_sal > 40 | raw_sal < 25) = nan;
                                raw_temp(raw_temp > 40 | raw_temp < -2) = nan;
                                
                                % Fill variables.  Note no interpolation in
                                % this version! %%%%%%%%%%%%%%%%%%%%
                                % prefill with nans
                                ocean.temp(:,dd) = nan(10000,1);
                                ocean.sal(:,dd) = nan(10000,1);
                                ocean.pres(:,dd) = nan(10000,1);
                                ocean.depth(:,dd) = nan(10000,1);
                                
                                ocean.temp(1:length(raw_temp),dd) = raw_temp;
                                ocean.sal(1:length(raw_sal),dd) = raw_sal;
                                ocean.depth(1:length(raw_depth),dd) = raw_depth;
                                ocean.pres(1:length(raw_pres),dd) = raw_pres;
                                
                                %% These variables are tacked on (useful info on contour location)
                                ocean.cont_vicinity(1,dd) = cont_local_index(min_loc); % near the right place on the 1800 m contour.
                                ocean.cont_distance(1,dd) = cont_dist_cumul(cont_local_index(min_loc));
                                ocean.contour_lon(1,dd) = cont_lon(cont_local_index(min_loc));
                                ocean.contour_lat(1,dd) = cont_lat(cont_local_index(min_loc));
                                ocean.dist_from_contour(1,dd) = Cdist;
                                
                                dd = dd+1;
                                nopres = 0;
                                
                                
                                
                            else
                                lost_due_to_pres_catch = lost_due_to_pres_catch+1;
                            end
                            
                        end
                        
                    else
                        lost_due_to_sal_catch = lost_due_to_sal_catch+1;
                    end % no salinity else end
                 
                end % third filter end
                
            else
            end % second filter end
        else
        end
        
    else
        lost_due_to_time_clause = lost_due_to_time_clause+1;
    end
    disp(aa)
end

%% Temporary save for sanity purposes.
save('D:\Work_computer_sync\OSNAP_postdoc\PAPERS_NEW\N_Atlantic_boundary\matlab\V3_050321\intermediate_saves\raw_profs_1000m','ocean','-v7.3');
