
clearvars; close all;
CP_general_parameters

figout = [figdir,'Cold_pool4/definition/']; 
if ~exist(figout,'dir'); mkdir(figout); end

load([matdir,'Saildrone_rawdata.mat'],'lon','lat','time','date0','sID','Ta','Wspd','SST','Rh')
lon(lon<0) = lon(lon<0) + 360;
time = time+ date0;

Q = function_humidity_conversion(Rh,Ta);
data0 = cat(2,Ta,Wspd,Q);
datname = {'Ta','Wspd','Q'};
method = {'dTa_S_D','dTa_N_D_B_C','dTa_P_I_R_A_T_A'};
nd = length(datname);


%% check data period
SD = unique(sID);
reg_filt = lon >= LONR1 & lon<= LONR2 & lat >= LATR1 & lat <= LATR2;
for s = 1:length(SD)
    
    filt = sID == SD(s) & reg_filt;
    if sum(filt) > 0
        disp(num2str(SD(s)))
        datestr(nanmin(time(sID == SD(s))))
        datestr(nanmax(time(sID == SD(s))))
        pause
    end
end

%% check data resolution
SD = unique(sID);
reg_filt = lon >= LONR1 & lon <= LONR2 & lat >= LATR1 & lat <= LATR2;
res = nan(length(SD),size(dat,2));
Nday = 0;
for s = 1:length(SD)
    disp(num2str(SD(s)))
    filt = sID == SD(s) & reg_filt;
    if sum(filt) == 0
        disp('No Data')
    else
        T = Ta(filt);
        S = SST(filt);
        R = Rh(filt);
        W = Wspd(filt);
        ti = time(filt);
        
        dat = cat(2,ti,T,S,W,R);
        
        for k = 1:size(dat,2)
            if k == 1
                dt = dat(2:end,k) - dat(1:end-1,k);
                res(s,k) = mode(dt)*86400/60;
                Nday = Nday + sum(isfinite(dat(:,k)))*res(s,k)*60/86400;
            else
                t = ti(isfinite(dat(:,k)));
                dt = t(2:end) - t(1:end-1);
                res(s,k) = mode(dt)*86400/60;
            end
        end

    end
end
disp('tres - tlen - Ta - SST - Wspd - Rh')
res
Nday

%% compute Ta slopes
SD = unique(sID);
reg_filt = lon >= LONR1 & lon<= LONR2 & lat >= LATR1 & lat <= LATR2;

TW_slope = [];
time_slope = [];
slope_sID = [];

dt = 60/86400; % time step
tshift = (4*60+30)/86400; % 4min 30 sec
mincount = [5,5,1]; % minimum data to be considered
for s = 1:length(SD)
    disp(num2str(SD(s)))
    filt = sID == SD(s) & reg_filt;
    if sum(filt) == 0
        disp('No Data')
    else
        t_s = time(filt);

        data_s = data0(filt,:);
        if length(unique(t_s)) ~= length(t_s)
            [t_s,ind] = unique(t_s);
            data_s = data_s(ind,:);
        end    

        % interpolate for missing data
        for k = 1:nd
            indnan = find(isnan(data_s(:,k)));
            ind = find(isfinite(data_s(:,k)));
            if ~isempty(ind) && ~isempty(indnan)
                data_s(indnan,k) = interp1(t_s(ind),data_s(ind,k),t_s(indnan));
            end
        end

        tvec = nanmin(t_s)+tshift : dt : nanmax(t_s);
        nt = length(tvec);

        slope = nan(nt,nd,3);

        % compute grad-10
        parfor i = 1:nt           
            filt = t_s >= tvec(i) - 5*dt & t_s < tvec(i) + 5*dt;
            if sum(filt) >= mincount(1)
                for k = 1:nd
                    p = polyfit(t_s(filt),data_s(filt,k),1);
                    slope(i,k,1) = p(1)/(60*24)*10;
                end

            end
        end
        % compute diff-mean10
        parfor i = 1:nt
    %         if mod(i,10) == 1

                filt1 = t_s >= tvec(i) - 10*dt & t_s < tvec(i);
                filt2 = t_s >= tvec(i) & t_s < tvec(i) + 10*dt;

                if sum(filt1) >= mincount(2) && sum(filt2) >= mincount(2)
                    for k = 1:nd
                        slope(i,k,2) = nanmean(data_s(filt2,k)) - nanmean(data_s(filt1,k));
                    end

                end
    %         end
        end
        % compute diff-mean2
        parfor i = 1:nt
    %         if mod(i,10) == 1
                filt1 = t_s >= tvec(i) - 7*dt & t_s < tvec(i) - 5*dt;
                filt2 = t_s >= tvec(i) + 3*dt & t_s < tvec(i) + 5*dt;
                if sum(filt1) >= mincount(3) && sum(filt2) >= mincount(3)
                    for k = 1:nd
                        slope(i,k,3) = nanmean(data_s(filt2,k)) - nanmean(data_s(filt1,k));  
                    end
                end
    %         end        

        end
        TW_slope = cat(1,TW_slope,slope);   
        time_slope = cat(2,time_slope,tvec);
        slope_sID = cat(1,slope_sID,ones(size(slope,1),1)*SD(s));
    end
 

end


save([matdir,'Cold_pools_grad_vs_diff_2018_2023.mat'],'TW_slope','time_slope','method','datname','slope_sID')


%% check   
load([matdir,'Cold_pools_grad_vs_diff_2018_2023.mat'],'TW_slope','time_slope','method','datname','slope_sID')  




   


    


        
    


    







































