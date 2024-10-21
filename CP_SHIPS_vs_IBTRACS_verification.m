
clearvars; close all;
saildrone_eparams

datadir = '/media/da/DA5T/Datasets/storm_tracks/SHIPS_predictors/';
outdir = [matdir,'TC_data/'];
if ~exist(outdir,'dir'); mkdir(outdir); end

%% load TC tracks locations
Ymin = 1998;
Ymax = 2023;
  
fname = [datadir,'VWS_extracted.csv'];
ships = readtable(fname);

load([matdir,'TC_data/','NAtl_TC_info_1998_2023.mat'],...
    'reg_code','region','TC_dist','loni','lati','timei','Vmax','TC_name','storm_code','R34','R50','R64','Rmax');
lon1 = loni;
lat1 = lati;
time1 = timei;
Vmax1 = Vmax;
code1 = storm_code;
name1 = {};
for k = 1:length(TC_name)
    name1{k} = deblank(TC_name(:,k)');
end


load([matdir,'TC_data/','SHIPS_VWS.mat'],'tc_name','tc_ID','timei','lati','loni','Vmax','SHRD','SHDC','SDDC','MSLP','d2land');
time2 = timei;
lon2 = loni;
lat2 = lati;
Vmax2 = Vmax;
name2 = tc_name;
code2 = tc_ID;

%% plot SHIPS track
lonmin = 250;
lonmax = 350;
latmin = -10;
latmax = 50;

lonR1 = LONR1;
lonR2 = LONR2;
latR1 = LATR1;
latR2 = LATR2;




figure('Position',[0 0 800 600]); hold on;
plot_projection_ROI(lonmin,lonmax,latmin,latmax);

clist2 = unique(code2);
h = [];
% plot TC tracks
for k = 1:length(clist2)
    filt = code2 == clist2(k);
    x = lon2(filt);
    y = lat2(filt);
    m_plot(x,y,'-','Color',my_color('light olive drab'),'LineWidth',1); 

end
plot_rectangle(lonR1,lonR2,latR1,latR2,'r');
% export_fig([figout,'SHIPS_tracks_map.png'],'-dpng','-r200')

%% plot distribution of Vmax
wscale = [0,33, 63, 82, 95, 112, 136, 500];
scale = {'TD','TS','Cat1','Cat2','Cat3','Cat4','Cat5'};
tc_life1 = zeros(length(scale),1);

filt = time1 < datenum(2023,1,1);
Vmax1 = Vmax1(filt);
code1 = code1(filt);
time1 = time1(filt);
dist_2_land = dist_2_land(filt);
clist1 = unique(code1);


for k = 1:length(clist1)
    filt = code1 == clist1(k) & isfinite(Vmax1) & dist_2_land > 0;
    t = time1(filt);
    dt = diff(t);
    dt(dt > 0.25) = 0;
    Vmax = Vmax1(filt);
    Vmax = 0.5*(Vmax(2:end) + Vmax(1:end-1));
    for i = 1:length(scale)       
        filt2 = Vmax > wscale(i) & Vmax <= wscale(i+1);
        if sum(filt2) > 0
            tc_life1(i) = tc_life1(i) + sum(dt(filt2));
        end
    end

end

clist2 = unique(code2);
tc_life2 = zeros(length(scale),1);
for k = 1:length(clist2)
    filt = code2 == clist2(k) & d2land > 0;
    t = time2(filt);
    dt = diff(t);
    dt(dt > 0.25) = 0;
    Vmax = Vmax2(filt);
    Vmax = 0.5*(Vmax(2:end) + Vmax(1:end-1));
    for i = 1:length(scale)       
        filt2 = Vmax > wscale(i) & Vmax <= wscale(i+1);
        if sum(filt2) > 0
            tc_life2(i) = tc_life2(i) + sum(dt(filt2));
        end
    end  
 
end    

figure; hold on;

yyaxis left
h(1) = bar((1:length(scale))-0.125,tc_life1,0.2,'FaceColor',my_color('dark green'))
h(2) = bar((1:length(scale))+0.125,tc_life2,0.2,'FaceColor',my_color('light blue'))
ylabel('tracked time(day)'); ylim([0,1500])

yyaxis right
h(3) = plot((tc_life2./tc_life1)*100,'-o','LineWidth',1.5);
legend(h(:),{'IBTrACS','SHIPS','SHIPS/IBTrACS'})
set(gca,'Xtick',1:length(scale),'XtickLabel',scale);
set(gca,'FontSize',12);
set(gca,'Ygrid','on');
ylabel('time ratio (%)'); ylim([40,100])

xlabel('Saffir-Simpson wind scale')

title('SHIPS vs. IBTrACS tracked time (1998-2022)')
set(gcf,'Color','w')
figout = [figdir,'Cold_pool_R1/SHIPS/']; 
export_fig([figout,'SHIPS_vs_IBTrACS_tracked_time_ocean.png'],'-dpng','-r200')

%% co-location check
Rmax_2 = nan(size(time2));
R34_2 = nan(size(time2));
R50_2 = nan(size(time2));
R64_2 = nan(size(time2));



count = 0;
err = nan(size(time2));
for i = 1:length(lon2)
    dlon = abs(lon1 - lon2(i));
    dlat = abs(lat1 - lat2(i));
    dt = abs(time1 - time2(i));
    dVmax = abs(Vmax1-Vmax2(i));
    ind = find(dlon < 0.001 & dlat < 0.001 & dt < 1/(24*60) & dVmax < 1,1,'first');
    
    
    if ~isempty(ind)
        count = count+1;
        Rmax_2(i) = Rmax(ind);
        R34_2(i) = R34(ind);
        R50_2(i) = R50(ind);
        R64_2(i) = R64(ind);
        err(i) = dlon(ind) + dlat(ind) + dt(ind);
    end
    
    
end
% save([matdir,'TC_data/','SHIPS_VWS.mat'],'Rmax_2','R34_2','R50_2','R64_2','-append')


%% compare TC name
figure; hold on;

tc_list = {};
for Y = 1998:2022
    filt1 = time1 >= datenum(Y,1,1,0,0,0) & time1 <= datenum(Y,12,31,23,59,59);
    filt2 = time2 >= datenum(Y,1,1,0,0,0) & time2 <= datenum(Y,12,31,23,59,59);
    ID1 = code1(filt1);
    ID2 = code2(filt2);
    
    list1 = name1(unique(ID1))';
    for k = 1:length(list1)
        if length(list1{k}) > 4
            list1{k} = list1{k}(1:4);
        end
    end
    list2 = unique(name2(filt2));
    

    lonY1 = lon1(filt1);
    latY1 = lat1(filt1);
    timeY1 = time1(filt1);
    VmaxY1 = Vmax1(filt1);
    

    lonY2 = lon2(filt2);
    latY2 = lat2(filt2);
    timeY2 = time2(filt2);
    VmaxY2 = Vmax2(filt2); 
    
    u_ID1 = unique(ID1);
    u_ID2 = unique(ID2);
    [common_name,ind1,ind2] = intersect(list1,list2);
    c_ID1 = u_ID1(ind1);
    c_ID2 = u_ID2(ind2);
    
    
    for i = 1:length(common_name)
        filt1 = ID1 == c_ID1(i);
        filt2 = ID2 == c_ID2(i);
        
        plot_projection_NWAtlantic
        m_scatter(lonY1(filt1),latY1(filt1),'bo')
        m_scatter(lonY2(filt2),latY2(filt2),'r^')
        title(common_name{i})
        pause
        cla
    end
        
end

%%
dat1 = cat(2,lon1,lat1,time1);
dat2 = cat(2,lon2,lat2,time2);

[lon,ind1,ind2] = intersect(lon1,lon2);
[time,ind1,ind2] = intersect(time1,time2);    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

