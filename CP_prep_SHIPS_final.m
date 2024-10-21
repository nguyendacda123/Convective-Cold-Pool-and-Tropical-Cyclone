
clearvars; close all;
CP_general_parameters

datadir = [orig_dir,'SHIPS/'];
outdir = [matdir,'TC_data/'];
if ~exist(outdir,'dir'); mkdir(outdir); end

%% load TC tracks locations
Ymin = 1998;
Ymax = 2023;
  
fname = [datadir,'SHIPS_extracted.csv'];
t = readtable(fname);


%% 

yymmdd = t.YYMMDD;
hh = t.UTC;
ymd = yymmdd;
yymmdd(ymd < 250000) = yymmdd(ymd < 250000) + 2e7;
yymmdd(ymd > 250000) = yymmdd(ymd > 250000) + 19e6;

timei = nan(length(hh),1);
for i = 1:length(hh)  
    timei(i) = datenum([num2str(yymmdd(i)),num2str(hh(i),'%.2d')],'yyyymmddhh');
end


filt = timei >= datenum(1998,1,1,0,0,0) & timei <= datenum(2023,12,31,23,59,59);
timei = timei(filt);
tc_ID = t.ID(filt);
% tc_ID = tc_ID - nanmin(tc_ID) + 1;
tc_name = t.Name(filt);
Vmax = t.Vmax(filt);
lati = t.Lat(filt);
loni = 360 - t.Lon(filt); % convert from west to east degree
MSLP = t.MSLP(filt);
SHRD = t.SHRD(filt);
SHDC = t.SHDC(filt);
SDDC = t.SDDC(filt);

% re-arrange data according to tc_ID



save([outdir,'SHIPS_VWS.mat'],'tc_name','tc_ID','timei','lati','loni','Vmax','SHRD','SHDC','SDDC','MSLP');


%% adding radius of maximum winds from IBTRACS

load([matdir,'TC_data/','NAtl_TC_info_1998_2023.mat'],'loni','lati','timei','R34','R50','R64','Rmax','dist_2_land');
lon1 = loni;
lat1 = lati;
time1 = timei;



load([matdir,'TC_data/','SHIPS_VWS.mat'],'timei','lati','loni');
time2 = timei;
lon2 = loni;
lat2 = lati;



Rmax_2 = nan(size(time2));
R34_2 = nan(size(time2));
R50_2 = nan(size(time2));
R64_2 = nan(size(time2));
d2land = nan(size(time2));

count = 0;
err = nan(size(time2));
for i = 1:length(lon2)
    dlon = abs(lon1 - lon2(i));
    dlat = abs(lat1 - lat2(i));
    dt = abs(time1 - time2(i));
    ind = find(dlon < 0.001 & dlat < 0.001 & dt < 1/(24*60),1,'first');
    
    
    if ~isempty(ind)
        count = count+1;
        Rmax_2(i) = Rmax(ind);
        R34_2(i) = R34(ind);
        R50_2(i) = R50(ind);
        R64_2(i) = R64(ind);
        d2land(i) = dist_2_land(ind);
        err(i) = dlon(ind) + dlat(ind) + dt(ind);
    end
    
    
end
Rmax = Rmax_2;
R34 = R34_2;
R50 = R50_2;
R64 = R64_2;
save([matdir,'TC_data/','SHIPS_VWS.mat'],'Rmax','R34','R50','R64','d2land','-append')


%% verification on shear vector
load([matdir,'TC_data/','SHIPS_VWS.mat'],'tc_name','tc_ID','timei','lati','loni','Vmax','SHRD','SHDC','SDDC','MSLP');

Q = ones(length(SDDC),1);
Q(SDDC > 90 & SDDC <= 180) = 2;
Q(SDDC > 180 & SDDC <= 270) = 3;
Q(SDDC > 270) = 4;
IDs = unique(tc_ID)

tdiff = [];
SDDC_diff = [];
Qdiff = [];
for k = 1:length(IDs)
    filt = tc_ID == IDs(k);
    Qdiff = diff(Q(filt));    
    tdiff = diff(timei(tc_ID == IDs(k)));
    SDDC_k = SDDC(filt);
    
    Adiff = nan(length(Qdiff),1);
    for i = 1:length(tdiff)
        if tdiff(i) == 0.25
            if abs(Qdiff(i)) <=1
                Adiff(i) = SDDC_k(i+1) - SDDC_k(i);
            elseif abs(Qdiff(i)) == 3
                vals = SDDC_k(i:i+1);
                vals(vals < 90) = vals(vals < 90) + 360;
                Adiff(i) = diff(vals);
            end
        end
    end
            
    
    SDDC_diff = cat(1,SDDC_diff,Adiff);
    
end

figure;
hist(SDDC_diff)
prctile(SDDC_diff,98)




















































