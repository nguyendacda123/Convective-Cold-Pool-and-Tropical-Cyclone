
clearvars; close all;
CP_general_parameters


%% load raw data
load([matdir,'TC_data/','NAtl_TC_info_1998_2023.mat'],...
    'reg_code','region','TC_dist','loni','lati','timei','Vmax',...
    'Vtrans','storm_code','reg_code','Ylabel','dist_2_land','TCI',...
    'R34','R50','R64','Rmax','TCI_period','TC_name','USA_R_meaning');


% process PIRATA data ----------------------------------------
load([matdir,'PIRATA_data.mat'],'time','sID','sname','lon','lat','Ta','sst','wspd');
load([matdir,'HF/pirata_HF.mat'],'shflx','lhflx');

P_glon = [];
P_glat = [];
for s = 1:max(sID)
    filt = sID == s;
    P_glon = cat(1,P_glon,ones(sum(filt),1)*lon(s));
    P_glat = cat(1,P_glat,ones(sum(filt),1)*lat(s));
end
P_gtime = time;
P_gTa = Ta;
P_gID = sID;
P_gcode = ones(size(P_glon));
P_SST = sst;
P_wspd = wspd;
P_shflx = shflx;
P_lhflx = lhflx;

% process NDBC data -----------------------------------------
datadir = [orig_dir, 'NDBC_buoys/'];
maxID = nanmax(sID);
load([matdir,'NDBC_offshore_buoys_raw_data.mat'],'N_gtime','N_gTa','N_RH',...
    'N_gID','N_glon','N_glat','N_SST','N_WSPD','lons_10','lats_10','sname_10')
load([matdir,'HF/ndbc_HF.mat'],'shflx','lhflx');
N_shflx = shflx;
N_lhflx = lhflx;
N_wspd = N_WSPD;
N_gID = N_gID + maxID;
N_gcode = ones(size(N_glon))*2;


% process saildrone data ---------------------------------------
load([matdir,'Saildrone_rawdata.mat'],'lon','lat','time','date0','sID','Ta','SST','Wspd');
load([matdir,'HF/saildrone_HF.mat'],'shflx','lhflx');
lon(lon < 0) = lon(lon<0) + 360;
filt = lon >= LONR1 & lon<= LONR2 & lat >= LATR1 & lat <= LATR2;
% filt = mod(sID,1e4) >= 2021;
S_gID = sID(filt);
S_glon = lon(filt);
S_glat = lat(filt);
S_gtime = time(filt);
S_gTa = Ta(filt);

S_shflx = shflx(filt);
S_lhflx = lhflx(filt);
S_SST = SST(filt);
S_wspd = Wspd;

u_gID = unique(S_gID);
ID = S_gID;
maxID = nanmax(N_gID);

for k = 1:length(u_gID)
    S_gID(ID == u_gID(k)) = k+maxID;
end
 
S_gcode = ones(size(S_glon))*3;



gtime = cat(1,P_gtime,N_gtime,S_gtime);
glon = cat(1,P_glon,N_glon,S_glon);
glat = cat(1,P_glat,N_glat,S_glat);
gID = cat(1,P_gID,N_gID,S_gID);
gTa = cat(1,P_gTa,N_gTa,S_gTa);
gSST = cat(1,P_SST,N_SST,S_SST);
gWspd = cat(1,P_wspd,N_wspd,S_wspd);
gcode = cat(1,P_gcode,N_gcode,S_gcode);
g_shflx = cat(1,P_shflx,N_shflx,S_shflx);
g_lhflx = cat(1,P_lhflx,N_lhflx,S_lhflx);
glon(glon<=0) = glon(glon<=0) + 360;



% filter missing Ta data
filt = isfinite(gTa);
gtime = gtime(filt);
glon = glon(filt);
glat = glat(filt);
gcode = gcode(filt);
gTa = gTa(filt);
gID = gID(filt);
gSST = gSST(filt);
gWspd = gWspd(filt);
g_shflx = g_shflx(filt);
g_lhflx = g_lhflx(filt);





%% load CP data
for opt = 2:3
    if opt == 1
        tag = '';
    elseif opt == 2
        tag = '_2prctile';
    else
        tag = '_0.5prctile';
    end

    saildrone = load([matdir,'Saildrone_cold_pools_2018_2023_3Tcrits',tag,'.mat'],'CP_period','CP_tstart','CP_tend','CP_sID','Tslope',...
        'CP_lon','CP_lat','dat_sample','dat_units','datname','method');
    S_clon = saildrone.CP_lon{1,1};
    S_clat = saildrone.CP_lat{1,1};
    S_ctstart = saildrone.CP_tstart{1,1};
    S_ctend = saildrone.CP_tend{1,1};
    dat_sample = saildrone.dat_sample{1,1};
    S_cSST = squeeze(nanmean(dat_sample(1:30,5,:),1));
    S_cID = saildrone.CP_sID{1,1};


    pirata = load([matdir,'PIRATA_buoys_cold_pools',tag,'.mat'],'lon','lat','sname',...
        'CP_tstart','CP_tend','CP_sID','CP_period','dat_sample','crit','datname');

    P_clon = nan(size(pirata.CP_tstart));
    P_clat = nan(size(pirata.CP_tstart));
    for k = 1:length(pirata.lon)
        filt = pirata.CP_sID == k;
        P_clon(filt) = pirata.lon(k);
        P_clat(filt) = pirata.lat(k);
    end
    P_ctstart = pirata.CP_tstart;
    P_ctend = pirata.CP_tend;
    P_cSST = squeeze(nanmean(pirata.dat_sample(1:30,4,:)));
    P_cID = pirata.CP_sID;

    ndbc = load([matdir,'NDBC_buoys_cold_pools_offshore',tag,'.mat'],'lons_10','lats_10','sname_10',...
        'CP_tstart','CP_tend','CP_sID','CP_period','dat_sample','datname','wlevel_obs_payload','carlcoos_payload','scoop_payload');

    N_clon = nan(size(ndbc.CP_tstart));
    N_clat = nan(size(ndbc.CP_tstart));
    for k = 1:length(ndbc.lons_10)
        filt = ndbc.CP_sID == k;
        N_clon(filt) = ndbc.lons_10(k);
        N_clat(filt) = ndbc.lats_10(k);
    end
    N_ctstart = ndbc.CP_tstart;
    N_ctend = ndbc.CP_tend;
    N_cSST = squeeze(nanmean(ndbc.dat_sample(1:30,4,:)));
    N_cID = ndbc.CP_sID;

    % join data ---------------------------------
    CP_sst = cat(1,P_cSST,N_cSST,S_cSST);
    CP_lon = cat(1,P_clon,N_clon,S_clon);
    CP_lat = cat(1,P_clat,N_clat,S_clat);

    CP_lon(CP_lon < 0) = CP_lon(CP_lon < 0) + 360;

    % CP_tstart = cat(1,P_ctstart,N_ctstart,S_ctstart);
    dt1 = 120/86400;
    dt2 = 600/86400;
    CP_tstart = cat(1,P_ctstart-dt1,N_ctstart-dt2,S_ctstart);


    CP_tend = cat(1,P_ctend,N_ctend,S_ctend);
    CP_period = (CP_tend - CP_tstart)*24*60; % in minutes
    CP_code = cat(1,ones(size(P_clon)),ones(size(N_clon))*2,ones(size(S_clon))*3);



    maxID = max(P_gID);

    uID = unique(N_cID);
    ID = N_cID;
    for k = 1:length(uID)
        filt = ID == uID(k);
        N_cID(filt) = uID(k) + maxID;
    end


    maxID = max(N_gID);
    uID = unique(S_cID);
    ID = S_cID;
    for k = 1:length(uID)
        filt = ID == uID(k);
        S_cID(filt) = k + maxID;
    end

    CP_ID = cat(1,P_cID,N_cID,S_cID);




    %% Interpolate TC locations to CP time
    % at time t, and location (loni,lati) are there CPs within 500km
    dt = 600/86400;
    nt = length(timei);
    N_gauge = max(gID); % number of observing system
    lon_center = [];
    lat_center = [];
    time_center = [];
    Vmax_center = [];
    Vtrans_center = [];
    TCI_center = [];
    R34_center = [];
    R50_center = [];
    R64_center = [];
    Rmax_center = [];
    d2land_center=[];

    lon_gauge = [];
    lat_gauge = [];
    dist_gauge = [];
    angle_gauge = [];
    quadrant_gauge = [];
    sst_gauge_pre = [];
    sst_gauge_post = [];

    wspd_gauge_pre = [];
    wspd_gauge_post = [];

    Ta_gauge_pre = [];
    Ta_gauge_post = [];

    shflx_gauge_pre = []; % pre-CP
    shflx_gauge_post = []; % post-CP

    lhflx_gauge_pre = [];
    lhflx_gauge_post = [];
    code_gauge = []; % 1 = Pirata, 2 = NDBC, 3 = saildrone
    phase_gauge = []; % -1 = pre-TC, 1 = post-TC


    CP_gauge = []; % 1 = exist, 0 = not exist
    CP_presst = []; % 30min averaged SST before CP
    period = [];
    TC_code = [];



    count = 0;
    %%
    for s = 1:max(storm_code)
        s
        tc_filt = storm_code == s;
        if sum(tc_filt) > 1
            tc_lon = loni(tc_filt);
            tc_lat = lati(tc_filt);
            tc_time = timei(tc_filt);
            tc_intense = Vmax(tc_filt);
            tc_trans = Vtrans(tc_filt);
            tc_R34 = R34(:,tc_filt);
            tc_R50 = R50(:,tc_filt);
            tc_R64 = R64(:,tc_filt);
            tc_Rmax = Rmax(tc_filt);
            tc_d2land=dist_2_land(tc_filt);
            dVmax_dt = TCI(tc_filt);

            % interpolate TC locations to every 10min
            tvec = [nanmin(tc_time):dt:nanmax(tc_time)]';
            lonc = interp1(tc_time,tc_lon,tvec);
            latc = interp1(tc_time,tc_lat,tvec);

            Vfilt = isfinite(tc_intense);
            if sum(Vfilt) > 1
                Vmaxc = interp1(tc_time(Vfilt),tc_intense(Vfilt),tvec);
            else
                Vmaxc = nan(size(tvec));
            end

            Vtransc = interp1(tc_time,tc_trans,tvec);
            d2landc = interp1(tc_time,tc_d2land,tvec);
            R34c = interp1(tc_time,tc_R34',tvec);
            R50c = interp1(tc_time,tc_R50',tvec);
            R64c = interp1(tc_time,tc_R64',tvec);
            Rmaxc = interp1(tc_time,tc_Rmax,tvec);
            dVdt = interp1(tc_time,dVmax_dt,tvec);
            nt = length(tvec);


            R = nan(nt,N_gauge);
            phi = nan(nt,N_gauge);
            quad = nan(nt,N_gauge);

            Lon_g = nan(nt,N_gauge);
            Lat_g = nan(nt,N_gauge);
            code_g = nan(nt,N_gauge);
            phase_g = ones(nt,N_gauge);
            sst_g_pre = nan(nt,N_gauge);
            sst_g_post = nan(nt,N_gauge);
            wspd_g_pre = nan(nt,N_gauge);
            wspd_g_post = nan(nt,N_gauge);
            Ta_g_pre = nan(nt,N_gauge);
            Ta_g_post = nan(nt,N_gauge);

            CP_g = zeros(nt,N_gauge);
            sst_cp = nan(nt,N_gauge);
            dist_500 = zeros(nt,N_gauge);
            dist_200 = zeros(nt,N_gauge);

            shflx_g_pre = nan(nt,N_gauge);
            shflx_g_post = nan(nt,N_gauge);
            lhflx_g_pre = nan(nt,N_gauge);
            lhflx_g_post = nan(nt,N_gauge);

            % find active gauges
            filt2 = gtime >= tvec(1) - dt & gtime <= tvec(end) + dt;

            if sum(filt2 > 0)
                gtime_s = gtime(filt2);
                glon_s = glon(filt2);
                glat_s = glat(filt2);
                gID_s = gID(filt2);
                gcode_s = gcode(filt2);
                gsst_s = gSST(filt2);
                gshflx_s = g_shflx(filt2);
                glhflx_s = g_lhflx(filt2);
                gWspd_s = gWspd(filt2);
                gTa_s = gTa(filt2);
                for i = 1:nt-1
                    % find all active gauges within 20-min
                    g_filt = gtime_s >= tvec(i) - 1.01*dt & gtime_s <= tvec(i) + 1.01*dt;               
                    if sum(g_filt) > 0
                        g_filt2 = gtime_s >= tvec(i) - 3.01*dt & gtime_s < tvec(i);
                        g_filt3 = gtime_s >= tvec(i) & gtime_s <= tvec(i) + 3.01*dt; % for heat flux
                        ID_g = gID_s(g_filt);
                        lon = glon_s(g_filt);
                        lat = glat_s(g_filt);
                        code = gcode_s(g_filt); 

                        id_pre = gID_s(g_filt2);
                        id_post = gID_s(g_filt3);

                        SST_pre = gsst_s(g_filt2);
                        SST_post = gsst_s(g_filt3);

                        wspd_pre = gWspd_s(g_filt2);
                        wspd_post = gWspd_s(g_filt3);

                        Ta_pre = gTa_s(g_filt2);
                        Ta_post = gTa_s(g_filt3);

                        sflx_pre = gshflx_s(g_filt2);
                        lflx_pre = glhflx_s(g_filt2);

                        sflx_post = gshflx_s(g_filt3);
                        lflx_post = glhflx_s(g_filt3); 


                        [uID_g,ind] = unique(ID_g);
                        lon = lon(ind);
                        lat = lat(ind);
                        code = code(ind);
                        center1 = [lonc(i),latc(i)]; 

                        % to avoid stationary TC (center1 = center2)
                        flag = 0;
                        for j = i+1:nt
                            if lonc(j) ~= lonc(i) || latc(j) ~= latc(i) % center 2 differs from center1
                                center2 = [lonc(j),latc(j)];
                                flag = 1;
                                break
                            end
                        end


                        if flag == 1
                            for j = 1:length(uID_g)                       
                                [R(i,uID_g(j)), phi(i,uID_g(j)), quad(i,uID_g(j))] = HR_function_TC_coordinates(center1,center2,[lon(j),lat(j)]);                            
                            end
                        end
                        Lon_g(i,uID_g) = lon;
                        Lat_g(i,uID_g) = lat;
                        code_g(i,uID_g) = code;

                        % compute SST, wspd, Ta, fluxes
                        for j = 1:length(uID_g)
                            filt = id_pre == uID_g(j); % sst value associated with the station
                            sst_g_pre(i,uID_g(j)) = nanmean(SST_pre(filt));                      
                            wspd_g_pre(i,uID_g(j)) = nanmean(wspd_pre(filt)); 
                            Ta_g_pre(i,uID_g(j)) = nanmean(Ta_pre(filt));                         
                            shflx_g_pre(i,uID_g(j)) = nanmean(sflx_pre(filt));
                            lhflx_g_pre(i,uID_g(j)) = nanmean(lflx_pre(filt));

                            filt = id_post == uID_g(j);
                            sst_g_post(i,uID_g(j)) = nanmean(SST_post(filt));                      
                            wspd_g_post(i,uID_g(j)) = nanmean(wspd_post(filt)); 
                            Ta_g_post(i,uID_g(j)) = nanmean(Ta_post(filt)); 
                            shflx_g_post(i,uID_g(j)) = nanmean(sflx_post(filt));
                            lhflx_g_post(i,uID_g(j)) = nanmean(lflx_post(filt)); 
                        end

                    end

                    % find active cold pools
                    cp_filt = CP_tstart <= tvec(i) & CP_tend >= tvec(i);
                    if sum(cp_filt) > 0
                        [ID_c, ind] = unique(CP_ID(cp_filt));
                        sst = CP_sst(cp_filt);
                        sst2 = sst(ind);
                        if length(intersect(uID_g,ID_c)) == length(ID_c)
                            CP_g(i,ID_c) = 1;
                            sst_cp(i,ID_c) = sst2;

                        else 
                            ind = intersect(uID_g,ID_c);                       
                            if ~isempty(ind)
                                CP_g(i,ind) = 1;
                            end 
                            count = count + length(setdiff(ID_c,uID_g));
                        end
                    end                

                end
            end



            lon_center = cat(1,lon_center,lonc);
            lat_center = cat(1,lat_center,latc);
            time_center = cat(1,time_center,tvec);        
            Vmax_center = cat(1,Vmax_center,Vmaxc);
            Vtrans_center = cat(1,Vtrans_center,Vtransc);
            R34_center = cat(1,R34_center,R34c);
            R50_center = cat(1,R50_center,R50c);
            R64_center = cat(1,R64_center,R64c); 
            Rmax_center = cat(1,Rmax_center,Rmaxc);
            d2land_center = cat(1,d2land_center,d2landc);

            TCI_center = cat(1,TCI_center,dVdt);

            lon_gauge = cat(1,lon_gauge,Lon_g);
            lat_gauge = cat(1,lat_gauge,Lat_g);
            dist_gauge = cat(1,dist_gauge,R);
            angle_gauge = cat(1,angle_gauge,phi);
            quadrant_gauge = cat(1,quadrant_gauge,quad);
            code_gauge = cat(1,code_gauge,code_g);


            sst_gauge_pre = cat(1,sst_gauge_pre,sst_g_pre);
            wspd_gauge_pre = cat(1,wspd_gauge_pre,wspd_g_pre);
            Ta_gauge_pre = cat(1,Ta_gauge_pre,Ta_g_pre);        
            lhflx_gauge_pre = cat(1,lhflx_gauge_pre,lhflx_g_pre);
            shflx_gauge_pre = cat(1,shflx_gauge_pre,shflx_g_pre);

            sst_gauge_post = cat(1,sst_gauge_post,sst_g_post);
            wspd_gauge_post = cat(1,wspd_gauge_post,wspd_g_post);
            Ta_gauge_post = cat(1,Ta_gauge_post,Ta_g_post);        
            lhflx_gauge_post = cat(1,lhflx_gauge_post,lhflx_g_post);
            shflx_gauge_post = cat(1,shflx_gauge_post,shflx_g_post);  

            phase_gauge = cat(1,phase_gauge,phase_g);
            CP_gauge = cat(1,CP_gauge,CP_g);
            CP_presst = cat(1,CP_presst,sst_cp);
            TC_code = cat(1,TC_code,ones(length(tvec),1)*s);



        end
    end


    save([matdir,'Active_gauges_and_CPs_in_TC_coordinates',tag,'.mat'],'lon_center','lat_center',...
        'time_center','Vmax_center','Vtrans_center','TCI_center','lon_gauge','lat_gauge','dist_gauge',...
        'code_gauge','angle_gauge','quadrant_gauge','TC_code','phase_gauge','CP_gauge','CP_presst',...
        'shflx_gauge_pre','lhflx_gauge_pre',...
        'shflx_gauge_post','lhflx_gauge_post','R34_center','R50_center','R64_center','Rmax_center','USA_R_meaning',...
        'wspd_gauge_pre','wspd_gauge_post','sst_gauge_pre','sst_gauge_post',...
        'Ta_gauge_pre','Ta_gauge_post','shflx_gauge_pre','lhflx_gauge_pre',...
        'shflx_gauge_post','lhflx_gauge_post','d2land_center')
end
% save([matdir,'Active_gauges_and_CPs_in_TC_coordinates.mat'],'Rmax_center','-append')  








































