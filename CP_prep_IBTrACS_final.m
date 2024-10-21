
clearvars; close all;
CP_general_parameters

datadir = '/media/da/DA5T/matdir/storm_tracks/IBTrACKS_2024/';
outdir = [matdir,'TC_data/'];
if ~exist(outdir,'dir'); mkdir(outdir); end

%% load TC tracks locations
Ymin = 1998;
Ymax = 2023;
region = {'North_Atlantic';  'Western_Pacific'; 'East_Pacific'; 'North_Indian'; 
             'South_Indian'; 'South_Pacific'};   

lon2 = [];
lat2 = [];
t2 = [];
TC_dist = []; 
TC_name = [];
Vmax2 = [];
Vtrans2 = [];
storm_code2 = [];
R34 = [];
R50 = [];
R64 = [];
Rmax = [];


reg_code = [];
Ylabel = [];
land_dist = [];

TCI_period = 24/24; % day
TCI = [];


for s = 1%:length(region)
    
    for Y = Ymin:Ymax
        Y    
        load([datadir,region{s},'_IBTRACKS_Y',num2str(Y),'.mat'],'numstorm','storm_name','loni','lati','timei',...
        'Vmaxi','Vtransi','Pmini','dist_2_land','land_fall','storm_code','USA_Rmax','USA_R34','USA_R50','USA_R64','USA_R_meaning');
    

        loni(loni<0) = loni(loni<0) + 360;
        if isempty(storm_code2)
            prev_ind = 0;
        else
            prev_ind = nanmax(storm_code2);
        end

        for k = 1:numstorm
            filt = storm_code == k;                
            if sum(filt) >= 2
                lon1 = loni(filt);
              
                lat1 = lati(filt);
                Vmax1 = Vmaxi(filt);
                Vtrans1 = Vtransi(filt);
                R34_1 = USA_R34(:,filt);
                R50_1 = USA_R50(:,filt);
                R64_1 = USA_R64(:,filt);
                Rmax_1 = USA_Rmax(filt);
                t1 = timei(filt);
                
                dist1 = nan(length(lon1)-1,1);
                for i = 1:length(lon1)-1
                    dist1(i) = lldistkm([lat1(i),lon1(i)],[lat1(i+1),lon1(i+1)])*1000;                                        
                end
                dist2 = cat(1,dist1(1)/2,0.5*(dist1(1:end-1,1)+dist1(2:end,1)),dist1(end)/2);
                TC_dist = cat(1,TC_dist,dist2);
                
%                 if max(dist2) > 2e6
%                     plot(lon1,lat1,'-ro'); pause;
%                 end
                lon2 = cat(1,lon2,lon1);
                lat2 = cat(1,lat2,lat1);                
                Vmax2 = cat(1,Vmax2,Vmax1);
                Vtrans2 = cat(1,Vtrans2,Vtrans1);
                R34 = cat(2,R34,R34_1);
                R50 = cat(2,R50,R50_1);
                R64 = cat(2,R64,R64_1);
                Rmax = cat(1,Rmax,Rmax_1);
                
                land_dist = cat(1,land_dist,dist_2_land(filt));
                t2 = cat(1,t2,t1);
                
                
                code = ones(size(lon1))*s;
                if s == 2
                    code(lon1 >= 180) = s+1;
                elseif s == 3
                    code(lon1 <= 180) = s-1;
                elseif code == 5
                    code(lon1 >= 140) = s+1;
                elseif code == 6
                    code(lon1 <= 140) = s-1;
                end
                
                reg_code = cat(1,reg_code,code);
                

                
                Ylabel = cat(1,Ylabel,ones(size(lon1))*Y);
                
                
                storm_code2 = cat(1,storm_code2,storm_code(filt)+prev_ind);
                TC_name = cat(2,TC_name,storm_name(:,k));
              
                % compute TC intensification rate

                rate = nan(size(t1));
                for i = 1:length(t1)
                    filt = t1 < t1(i) + TCI_period + 1/24 & t1 >= t1(i) & isfinite(Vmax1);
                    if sum(filt) > 2
                        p = polyfit(t1(filt),Vmax1(filt),1); % unit = knot/day
                        rate(i) = p(1);
                  
                    end

                end

                TCI = cat(1,TCI,rate);                  
                
            
            end
        end
        
    end 
end

Vmax = Vmax2;
Vtrans = Vtrans2;
loni = lon2;
lati = lat2;
timei = t2;
storm_code = storm_code2;
dist_2_land = land_dist;




save([outdir,'NAtl_TC_info_1998_2023.mat'],...
    'reg_code','region','TC_dist','loni','lati','timei','Vmax',...
    'Vtrans','storm_code','reg_code','Ylabel','dist_2_land','TCI',...
    'TCI_period','TC_name','USA_R_meaning','R34','R50','R64','Rmax');















