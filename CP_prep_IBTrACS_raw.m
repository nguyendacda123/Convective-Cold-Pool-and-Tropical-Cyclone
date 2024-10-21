
clearvars;
close all;

CP_general_parameters

indir = [orig_dir,'storm_tracks/IBTrACS_2024/'];
outdir = [matdir,'storm_tracks/IBTrACKS_2024/']; 
if ~exist(outdir,'dir'); mkdir(outdir); end
% 
fname = [indir,'IBTrACS.ALL.v04r00.nc']; 
% ncdisp(fname); pause

ref_date = datenum(1858,11,17);


numobs = ncread(fname,'numobs');


time = ncread(fname,'time');
time = time + ref_date;

lon = ncread(fname,'lon');
lat = ncread(fname,'lat');

basin = ncread(fname,'basin');
basin_code1 = {'EP','NA','NI','SA','SI','SP','WP'};
basin_name = {'East_Pacific'; 'North_Atlantic'; 'North_Indian'; 
    'South_Atlantic'; 'South_Indian'; 'South_Pacific'; 'Western_Pacific'};

% subbasin = ncread(fname,'subbasin');
% subbasin_code1 = {'AS','BB','CP','CS','GM','NA','EA','WA','MM'};
% subbasin_name = {'Arabian_Sea Bay_of_Bengal'; 'Central_Pacific'; 
%     'Caribbean_Sea'; 'Gulf_of_Mexico'; 'North_Atlantic'; 'Eastern_Australia';
%     'Western_Australia'; 'No_subbasin_for_this_position'};

name = ncread(fname,'name');
Vmax = ncread(fname,'wmo_wind'); 
Vtrans = ncread(fname,'storm_speed');
Pmin = ncread(fname,'wmo_pres');
dist2land = ncread(fname,'dist2land'); % unit = km
landfall = ncread(fname,'landfall');

usa_r34 = ncread(fname,'usa_r34');
usa_r50 = ncread(fname,'usa_r50');
usa_r64 = ncread(fname,'usa_r64');
usa_rmax = ncread(fname,'usa_rmw');


R1 = usa_r34*1.852; % nautical miles to km
R2 = usa_r50*1.852;
R3 = usa_r64*1.852;
R4 = usa_rmax*1.852;
%%
verif = 0;
Ymin = 1998;
Ymax = 2023;
reg = squeeze(basin(:,1,:))';
for Y = Ymin:Ymax
    Y
    t1 = datenum(Y,1,1);
    t2 = datenum(Y,12,31);
    filt1 = time(1,:) >= t1 & time(1,:) <= t2;    
    
    for r = 1:length(basin_name)
        
        filt2 = sum(reg == basin_code1{r},2) == 2;
        
        filt = filt1&filt2'; 
        if verif == 1
            name(:,filt)'
            pause
        end
        
        N0 = numobs(filt');

        numstorm = length(N0);
        storm_name = name(:,filt);
        
        time0 = time(:,filt);
        lon0 = lon(:,filt);
        lat0 = lat(:,filt);        

        Vmax0 = Vmax(:,filt);
        Vtrans0 = Vtrans(:,filt);
        Pmin0 = Pmin(:,filt);
        dist2land0 = dist2land(:,filt);
        landfall0 = landfall(:,filt);
        TC_size1 = R1(:,:,filt);
        TC_size2 = R2(:,:,filt);
        TC_size3 = R3(:,:,filt);
        TC_size4 = R4(:,filt);
        
        % joining data to 1D vector
        
        timei = [];
        lati = [];
        loni = [];

        storm_code = [];

   
        Vmaxi = [];
        Vtransi = [];
        Pmini = [];
        dist_2_land = [];
        land_fall = [];
        
        USA_R34 = [];
        USA_R50 = [];
        USA_R64 = [];
        USA_Rmax = [];

        for k = 1:length(N0)
            N = N0(k);
            timei = cat(1,timei,time0(1:N,k));
            loni = cat(1,loni,lon0(1:N,k));
            lati = cat(1,lati,lat0(1:N,k));
            

 
            storm_code = cat(1,storm_code,ones(N,1)*k);
            Vmaxi = cat(1,Vmaxi,Vmax0(1:N,k));
            
            Vtrans2 = nan(N,1);
            if N > 1
                lon2 = cat(1,lon0(1,k),0.5*(lon0(2:N,k) + lon0(1:N-1,k)),lon0(N,k));
                lat2 = cat(1,lat0(1,k),0.5*(lat0(2:N,k) + lat0(1:N-1,k)),lat0(N,k));
                t2 = cat(1,time0(1,k),0.5*(time0(2:N,k) + time0(1:N-1,k)),time0(N,k));
             

                dist = [];
                for i = 1:length(lon2)-1
                    dist(i,1) = lldistkm([lat2(i),lon2(i)],[lat2(i+1),lon2(i+1)])*1000; % unit = m
                end


                dt = (t2(2:end) - t2(1:end-1))*86400; % unit = second

                Vtrans2 = dist./dt;
                Vtrans2(Vtrans2 > 70) = nan;
            end

            Vtransi = cat(1,Vtransi,Vtrans2);            
            USA_R34 = cat(2,USA_R34,TC_size1(:,1:N,k));
            USA_R50 = cat(2,USA_R50,TC_size2(:,1:N,k));
            USA_R64 = cat(2,USA_R64,TC_size3(:,1:N,k));
            USA_Rmax = cat(1,USA_Rmax,TC_size4(1:N,k));
            
%             Vtransi = cat(1,Vtransi,Vtrans0(1:N,k)*0.514444); % knot to m/s;
            Pmini = cat(1,Pmini,Pmin0(1:N,k));
            dist_2_land = cat(1,dist_2_land,dist2land0(1:N,k));
            land_fall = cat(1,land_fall,landfall0(1:N,k));
            
            
        end  
        USA_R_meaning = 'from 1 to 4 NE clockwise to NW';
        if length(Vmaxi) ~= length(USA_R34); pause; end
        fname = [basin_name{r},'_IBTRACKS_Y',num2str(Y),'.mat'];
        save([outdir,fname],'numstorm','storm_name','loni','lati','timei',...
            'Vmaxi','Vtransi','Pmini','dist_2_land','land_fall','USA_R34',...
            'USA_R50','USA_R64','USA_Rmax','storm_code','USA_R_meaning')
        
        
    end
end


%% verification


% r = 3; 
% figure; hold on;
% for Y = 2019:2021
%     fname = [basin_name{r},'_IBTRACKS_Y',num2str(Y),'.mat'];
%     load([outdir,fname],'numstorm','storm_name','loni','lati','timei',...
%         'Vmaxi','Vtransi','Pmini','dist_2_land','land_fall','storm_code')
%     [sum(isfinite(Vmaxi)), sum(isfinite(loni))]
%     
%     plot_projection_ROI(min(loni),max(loni),min(lati),max(lati));
%     filt = isfinite(Vmaxi);
%     m_scatter(loni(filt),lati(filt),'bo');
%     m_scatter(loni(~filt),lati(~filt),'rs');
%     title(num2str(Y));
%     
%     
%     pause
%     cla
% end




