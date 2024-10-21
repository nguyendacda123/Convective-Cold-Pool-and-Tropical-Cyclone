
clearvars; close all;

CP_general_parameters
for opt = 1:3  % sensitivity to cold pool criteria
    if opt == 1
        Tcrits = [-0.70, -0.57, -0.65]; % 1 percentile
        tag = '_1prctile';
    elseif opt == 2

        Tcrits = [-0.41, -0.33, -0.38]; % 2 percentile
        tag = '_2prctile';
    else
        Tcrits = [-1.09, -0.85, -0.98]; % 0.5 percentile 
        tag = '_0.5prctile';
    end


    load([matdir,'Saildrone_rawdata.mat'],'lon','lat','time','date0','Y_label','SD_id','SD_Y',...
                                          'Ta','Rh','Pa','PAR','SWR','Wspd','Wdir','Uwind','Vwind','Gust',...
                                          'SST','SST_skin','SSS','Wav_T','Hs','C_O2','S_O2', ...
                                          'Cspd','Cdir','Chla','CDom','sID')
    tref = [-60:1:120]/(24*60); % time window centered at coldpool event
    dt = 60/86400; % 1 min -> day
    Q = function_humidity_conversion(Rh,Ta);
    data0 = cat(2,Ta,Wspd,Q,Rh,SST);
    datname0 = {'Ta','Wspd','Q','Rh','SST'};
    dat_units = {'o^C','m.s^-^1','g.m^-^3','%','^oC'};


    load([matdir,'Cold_pools_grad_vs_diff_2018_2023.mat'],'TW_slope','time_slope','method','datname','slope_sID')


    Tslope = squeeze(TW_slope(:,1,:));


    %% compute Ta slopes
    SD = unique(slope_sID);
    nm = length(method);
    Radt = [5, 10, 6]; % time searching radius


    ns = 9; % shifting step
    CP_period = cell(nm,ns);
    CP_tstart = cell(nm,ns);
    CP_lon = cell(nm,ns);
    CP_lat = cell(nm,ns);

    CP_tend = cell(nm,ns);
    CP_sID = cell(nm,ns);
    dat_sample = cell(nm,ns);

    % Tcrits = [-0.70, -0.57, -0.65]; % remember to update this in eParams
    diff_tstamp = {};


    for s = 1:length(SD)
        disp(num2str(SD(s)))

        % filter raw data
        filt1 = sID == SD(s);   
        t_s = time(filt1);
        lon_s = lon(filt1);
        lat_s = lat(filt1);    
        data_s = data0(filt1,:);    


        if length(unique(t_s)) ~= length(t_s)
            [t_s,ind] = unique(t_s);
            data_s = data_s(ind,:);
            lon_s = lon_s(ind);
            lat_s = lat_s(ind);
        end    

        % interpolate for missing data
        for k = 1:length(datname)-1
            indnan = find(isnan(data_s(:,k)));
            ind = find(isfinite(data_s(:,k)));
            if length(ind) > 0 && length(indnan) > 0
                data_s(indnan,k) = interp1(t_s(ind),data_s(ind,k),t_s(indnan));
            end
        end    


        % filter slope data
        filt2 = slope_sID == SD(s);
        tvec = time_slope(filt2);
        slope = Tslope(filt2,:);
        nt = length(tvec);


        for m = 1%:nm

            for j = 1%:ns

                flag = zeros(nt,1);
                slope2 = nan(size(slope));
                slope2(:,1) = slope(:,1); % no change for saildrone

                idx = j:10:length(tvec);
                slope2(idx,m) = slope(idx,m); % only change if m = 2

                % flag slopes < Tcrit
                for k = 1:length(tvec)            
                    if slope2(k,m) <= Tcrits(m)         
                        ind = find(t_s >= tvec(k) - dt*Radt(m) & t_s < tvec(k) + dt*Radt(m));           
                        flag(ind) = 1;
                    end
                end

                % compute number of cold pool events
                i = 0;
                count = 0;
                N = 0;
                ts = [];
                te = [];
                lon1 = [];
                lat1 = [];
                period = [];

                for i = 1:nt
                    if flag(i) == 1
                        if i < nt
                            count = count + 1;
                        else % end of series
                            count = count + 1;
                            N = N +1;
                            te(N) = t_s(nt);
                            ts(N) = t_s(nt-count+1); 
                            period(N) = (te(N)-ts(N))*24*60 + 1;
                            lon1(N) = nanmean(lon_s(nt-count+1:nt));
                            lat1(N) = nanmean(lat_s(nt-count+1:nt)); 

                            t3h = ts(N) + tref;
                            dat_sample{m,j} = cat(3,dat_sample{m,j},interp1(t_s,data_s,t3h));
    %                         if period(N) < 6; pause; end

                        end
                    else
                        if count > 0          
                            N = N +1;
                            te(N) = t_s(i-1);
                            ts(N) = t_s(i-count);
                            period(N) = (te(N)-ts(N))*24*60 + 1;

    %                         if period(N) < 6; pause; end                    
                            lon1(N) = nanmean(lon_s(i-count:i-1));
                            lat1(N) = nanmean(lat_s(i-count:i-1));
                            t3h = ts(N) + tref;
                            dat_sample{m,j} = cat(3,dat_sample{m,j},interp1(t_s,data_s,t3h));                                
                            count = 0;
                        end

                    end
                end
                CP_period{m,j} = cat(1,CP_period{m,j},period');
                CP_tstart{m,j} = cat(1,CP_tstart{m,j},ts');
                CP_tend{m,j} = cat(1,CP_tend{m,j},te');
                CP_sID{m,j} = cat(1,CP_sID{m,j},ones(N,1)*SD(s));
                CP_lon{m,j} = cat(1,CP_lon{m,j},lon1');
                CP_lat{m,j} = cat(1,CP_lat{m,j},lat1'); 
            end
        end



    end

    datname = datname0;
    % save([matdir,'Cold_pools_grad_vs_diff_2021_2023.mat'],'TW_slope','time_slope','method','datname','slope_sID')
    save([matdir,'Saildrone_cold_pools_2018_2023_3Tcrits',tag,'.mat'],'CP_period','CP_tstart','CP_tend','CP_sID','Tslope',...
        'CP_lon','CP_lat','dat_sample','dat_units','datname','method')
end
    
  
    


        
    


    






































