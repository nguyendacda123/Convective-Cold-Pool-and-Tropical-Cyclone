
clearvars; close all;
CP_general_parameters

figout = [figdir,'Cold_pool_R1/']; 
if ~exist(figout,'dir'); mkdir(figout); end



%% plot TC tracks, saildrone tracks and buoys locations
load([matdir,'Saildrone_rawdata.mat'],'lon','lat','sID');
lon(lon<0) = lon(lon<0) + 360;
filt = lon >= 0 & lon < 50;
lon(filt) = lon(filt) + 360;
load([matdir,'TC_data/','NAtl_TC_info_1998_2023.mat'],'loni','lati','storm_code','dist_2_land');
filt = loni >= 0 & loni < 50;
loni(filt) = loni(filt) + 360;
% load([matdir,'TC_data/','SHIPS_VWS.mat'],'loni','lati','tc_ID');
% storm_code = tc_ID;
% loni = 360 - loni;

pirata = load([matdir,'PIRATA_buoys_cold_pools_crit063.mat'],'lon','lat');
ndbc = load([matdir,'NDBC_offshore_buoys_metadata.mat'],'lons_10','lats_10');

load([matdir,'Active_gauges_and_CPs_in_TC_coordinates.mat'],'lon_center','lat_center',...
    'time_center','lon_gauge','lat_gauge','dist_gauge',...
    'code_gauge','TC_code','phase_gauge','CP_gauge','d2land_center')


%%

lonmin = 250;
lonmax = 350;
latmin = -10;
latmax = 50;

lonR1 = LONR1;
lonR2 = LONR2;
latR1 = LATR1;
latR2 = LATR2;


g_dist = reshape(dist_gauge,[],1);
code = reshape(repmat(TC_code,1,size(lon_gauge,2)),[],1);

loni(dist_2_land <= 0) = nan;

figure('Position',[0 0 800 600]); hold on;
plot_projection_ROI(lonmin,lonmax,latmin,latmax);


h = [];
% plot TC tracks
count = 0;
for k = 1:max(storm_code)
    filt = code == k;
    filt2 = storm_code == k;
    if nanmin(g_dist(filt)) < 1000
        x = loni(filt2);
        y = lati(filt2);
        h(1) = m_plot(x,y,'-','Color',my_color('light olive drab'),'LineWidth',1); 
        count = count + 1;
    end

end

% plot saildrone track
ID = unique(sID);
for k = 1:length(ID)
    if mod(ID(k),1e4) < 2023
        filt = sID == ID(k);
        x = lon(filt);
        y = lat(filt);
        h(2) = m_plot(x,y,'m-','LineWidth',1.5);
    end
end

plot_rectangle(lonR1,lonR2,latR1,latR2,'r');
% plot PIRATA location
x = pirata.lon; x(x<=0) = x(x<=0) + 360;
y = pirata.lat;
h(3) = m_plot(x,y,'rs','MarkerFaceColor','r');


% plot NDBC location
x = ndbc.lons_10; x(x<0) = x(x<0) + 360;
y = ndbc.lats_10;
h(4) = m_plot(x,y,'bo','MarkerFaceColor','b');

% legend(h(:),{'TC','Saildrone','PIRATA','NDBC'});
legend(h(:),{'TC','Saildrone','PIRATA','NDBC'});
set(gca,'FontSize',14);
set(gcf,'Color','w'); pause

export_fig([figout,'IBTrACS_tracks_map.png'],'-dpng','-r200')















































