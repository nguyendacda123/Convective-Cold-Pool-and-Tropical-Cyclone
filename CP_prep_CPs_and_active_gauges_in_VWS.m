
clearvars; close all;
CP_general_parameters


%% load raw data

for choice = 1:2
    if choice == 1
        opt = 'VWS';
    else
        opt = 'motion';
    end
    load([matdir,'TC_data/SHIPS_VWS.mat'],'tc_ID','timei','lati','loni','Vmax',...
        'SHRD','SHDC','SDDC','R34','Rmax','R50','R64','d2land');


    % process PIRATA data ----------------------------------------
    load([matdir,'PIRATA_data.mat'],'time','sID','sname','lon','lat','Ta','sst');
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
    N_gID = N_gID + maxID;
    N_gcode = ones(size(N_glon))*2;


    % process saildrone data ---------------------------------------
    load([matdir,'Saildrone_rawdata.mat'],'lon','lat','time','date0','sID','Ta','SST');
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
    g_shflx = g_shflx(filt);
    g_lhflx = g_lhflx(filt);





    %% load CP data
    tag = '';
    % tag = '_2prctile';
    % tag = '_0.5prctile';
    saildrone = load([matdir,'Saildrone_cold_pools_2018_2023_3Tcrits',tag,'.mat'],'CP_period','CP_tstart','CP_tend','CP_sID','Tslope',...
        'CP_lon','CP_lat','dat_sample','dat_units','datname','method');
    S_clon = saildrone.CP_lon{1,1};
    S_clat = saildrone.CP_lat{1,1};
    S_ctstart = saildrone.CP_tstart{1,1};
    S_ctend = saildrone.CP_tend{1,1};
    SST = saildrone.dat_sample{1,1};
    S_cSST = squeeze(nanmean(SST(1:30,6,:),1));
    S_cID = saildrone.CP_sID{1,1};


    pirata = load([matdir,'PIRATA_buoys_cold_pools',tag,'.mat'],'lon','lat','sname',...
        'CP_tstart','CP_tend','CP_sID','CP_period','dat_sample','crit');

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
        'CP_tstart','CP_tend','CP_sID','CP_period','dat_sample','wlevel_obs_payload','carlcoos_payload','scoop_payload');

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
    SHRD_center = [];
    SHDC_center = [];
    SDDC_center = [];
    Rmax_center = [];
    R34_center = [];
    R50_center = [];
    R64_center = [];
    d2land_center=[];

    lon_gauge = [];
    lat_gauge = [];
    dist_gauge = [];
    angle_gauge = [];
    quadrant_gauge = [];
    sst_gauge = [];
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
    for s = 1:max(tc_ID)
        s
        tc_filt = tc_ID == s;

        if sum(tc_filt) > 1 
            tc_lon = loni(tc_filt);
            tc_lat = lati(tc_filt);
            tc_time = timei(tc_filt);
            tc_intense = Vmax(tc_filt); 
            tc_SHRD = SHRD(tc_filt);
            tc_SHDC = SHDC(tc_filt);
            tc_SDDC = SDDC(tc_filt);

            tc_Rmax = Rmax(tc_filt);
            tc_R34 = R34(tc_filt);
            tc_R50 = R50(tc_filt);
            tc_R64 = R64(tc_filt);        
            tc_d2land=d2land(tc_filt);

            % interpolate TC locations to every 10min
            tvec = [nanmin(tc_time):dt:nanmax(tc_time)]';
            lonc = interp1(tc_time,tc_lon,tvec);
            latc = interp1(tc_time,tc_lat,tvec);
            Vmaxc = interp1(tc_time,tc_intense,tvec);
            SHRDc = interp1(tc_time,tc_SHRD,tvec);
            SHDCc = interp1(tc_time,tc_SHDC,tvec);
            Rmaxc = interp1(tc_time,tc_Rmax,tvec);
            R34c = interp1(tc_time,tc_R34,tvec);
            R50c = interp1(tc_time,tc_R50,tvec);
            R64c = interp1(tc_time,tc_R64,tvec);
            d2landc = interp1(tc_time,tc_d2land,tvec);
            nt = length(tvec);
            % exclude time difference > 6h
            tdiff = diff(tc_time);
            for k = 1:length(tc_time)-1
                if tdiff(k) > 0.25
                    filt = tvec > tc_time(k) & tvec < tc_time(k+1);
                    lonc(filt) = nan;
                    latc(filt) = nan;
                    Vmaxc(filt) = nan;
                    SHRDc(filt) = nan;
                    SHDCc(filt) = nan;
                    R34c(filt) = nan;
                    Rmaxc(filt) = nan;
                    R50c(filt) = nan;
                    R64c(filt) = nan;
                    d2lanc(filt) = nan;
                end
            end

            % special treatment for SDDC
            SDDCc = nan(nt,1);
            for k = 1:length(tc_time)-1
                filt = tvec > tc_time(k) & tvec < tc_time(k+1);
                if tdiff(k) <= 0.25
                    % between the downshear quadrants
                    if tc_SDDC(k) > 270 && tc_SDDC(k+1) <= 90
                        SDDCc(filt) = interp1([tc_time(k),tc_time(k+1)],...
                            [tc_SDDC(k),tc_SDDC(k+1)+360],tvec(filt));
                    elseif tc_SDDC(k) <= 90 && tc_SDDC(k+1) > 270
                        SDDCc(filt) = interp1([tc_time(k),tc_time(k+1)],...
                            [tc_SDDC(k),tc_SDDC(k+1)+360],tvec(filt));
                    % between two consecutive quadrant
                    elseif abs(tc_SDDC(k+1) - tc_SDDC(k)) <= 90
                        SDDCc(filt) = interp1([tc_time(k),tc_time(k+1)],...
                            [tc_SDDC(k),tc_SDDC(k+1)],tvec(filt));
                    end
                end
            end


            R = nan(nt,N_gauge);
            phi = nan(nt,N_gauge);
            quad = nan(nt,N_gauge);

            Lon_g = nan(nt,N_gauge);
            Lat_g = nan(nt,N_gauge);
            code_g = nan(nt,N_gauge);
            phase_g = ones(nt,N_gauge);
            sst_g = nan(nt,N_gauge);
            CP_g = zeros(nt,N_gauge);
            sst_cp = nan(nt,N_gauge);


            shflx_g_pre = nan(nt,N_gauge);
            lhflx_g_pre = nan(nt,N_gauge);
            shflx_g_post = nan(nt,N_gauge);
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

                        SST = gsst_s(g_filt2);
                        SST_id = gID_s(g_filt2);

                        flx_id_pre = gID_s(g_filt2);
                        sflx_pre = gshflx_s(g_filt2);
                        lflx_pre = glhflx_s(g_filt2);

                        flx_id_post = gID_s(g_filt3);
                        sflx_post = gshflx_s(g_filt3);
                        lflx_post = glhflx_s(g_filt3);                    

                        [uID_g,ind] = unique(ID_g);
                        lon = lon(ind);
                        lat = lat(ind);
                        code = code(ind);
                        center1 = [lonc(i),latc(i)]; 

                        % determining quadrants based on VWS or motion

                        if strcmp(opt,'VWS')
                            lon2 = lonc(i) + sin(SDDCc(i)/180*pi);
                            lat2 = latc(i) + cos(SDDCc(i)/180*pi);
                            center2 = [lon2,lat2];

                            for j = 1:length(uID_g)                       
                                [R(i,uID_g(j)), phi(i,uID_g(j)), quad(i,uID_g(j))] = HR_function_TC_coordinates(center1,center2,[lon(j),lat(j)]);                            
                            end                        
                        elseif strcmp(opt,'motion')
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
                        end

                        Lon_g(i,uID_g) = lon;
                        Lat_g(i,uID_g) = lat;
                        code_g(i,uID_g) = code;

                        % compute pre-30min SST
                        for j = 1:length(uID_g)
                            filt = SST_id == uID_g(j); % sst value associated with the station
                            sst_g(i,uID_g(j)) = nanmean(SST(filt));
                        end

                        % compute mean heat flux

                        for j = 1:length(uID_g)
                            filt = flx_id_pre == uID_g(j); % flux value associated with the station
                            shflx_g_pre(i,uID_g(j)) = nanmean(sflx_pre(filt));
                            lhflx_g_pre(i,uID_g(j)) = nanmean(lflx_pre(filt));

                            filt = flx_id_post == uID_g(j);
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

            Rmax_center = cat(1,Rmax_center,Rmaxc);
            R34_center = cat(1,R34_center,R34c);
            R50_center = cat(1,R50_center,R50c);
            R64_center = cat(1,R64_center,R64c);
            d2land_center = cat(1,d2land_center,d2landc);
            
            SHRD_center = cat(1,SHRD_center,SHRDc);
            SHDC_center = cat(1,SHDC_center,SHDCc);
            SDDC_center = cat(1,SDDC_center,SDDCc);

            lon_gauge = cat(1,lon_gauge,Lon_g);
            lat_gauge = cat(1,lat_gauge,Lat_g);
            dist_gauge = cat(1,dist_gauge,R);
            angle_gauge = cat(1,angle_gauge,phi);
            quadrant_gauge = cat(1,quadrant_gauge,quad);
            code_gauge = cat(1,code_gauge,code_g);
            sst_gauge = cat(1,sst_gauge,sst_g);
            lhflx_gauge_pre = cat(1,lhflx_gauge_pre,lhflx_g_pre);
            shflx_gauge_pre = cat(1,shflx_gauge_pre,shflx_g_pre);
            lhflx_gauge_post = cat(1,lhflx_gauge_post,lhflx_g_post);
            shflx_gauge_post = cat(1,shflx_gauge_post,shflx_g_post);        

            phase_gauge = cat(1,phase_gauge,phase_g);
            CP_gauge = cat(1,CP_gauge,CP_g);
            CP_presst = cat(1,CP_presst,sst_cp);
            TC_code = cat(1,TC_code,ones(length(tvec),1)*s);


        end
    end


    save([matdir,'Active_gauges_and_CPs_in_SHIPS_',opt,'.mat'],'lon_center','lat_center',...
        'time_center','Vmax_center','Rmax_center','R34_center','R50_center','R64_center',...
        'SHRD_center','SHDC_center','lon_gauge','lat_gauge','dist_gauge',...
        'code_gauge','angle_gauge','quadrant_gauge','TC_code','phase_gauge','CP_gauge','CP_presst',...
        'sst_gauge','shflx_gauge_pre','lhflx_gauge_pre',...
        'shflx_gauge_post','lhflx_gauge_post','d2land_center')   
end
% save([matdir,'Active_gauges_and_CPs_in_TC_coordinates.mat'],'Rmax_center','-append')  








































