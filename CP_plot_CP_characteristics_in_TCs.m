
clearvars; close all;
saildrone_eparams

figout = [figdir,'Cold_pool_R1/grad_vs_diff/']; 
if ~exist(figout,'dir'); mkdir(figout); end


% reference for wind conversion: http://www.systemsengineeringaustralia.com.au/download/WMO_TC_Wind_Averaging_27_Aug_2010.pdf
load([matdir,'Saildrone_cold_pools_2018_2023_3Tcrits.mat'],'CP_period','CP_tstart','CP_tend','CP_sID','Tslope',...
    'CP_lon','CP_lat','dat_sample','dat_units','datname','method','diff2_mis','diff2_far','diff10_mis','diff10_far')
data = dat_sample{1}(:,1:3,:);

load([matdir,'NDBC_buoys_cold_pools_offshore.mat'],'lons_10','lats_10','sname_10',...
    'CP_tstart','CP_tend','CP_sID','CP_period','dat_sample','datname','wlevel_obs_payload','carlcoos_payload','scoop_payload','crit')
dat_sample(:,2,:) = dat_sample(:,2,:)*1.11;
data = [];
data = cat(3,data,dat_sample(:,1:3,:));


% load([matdir,'PIRATA_buoys_cold_pools.mat'],'lon','lat','sname',...
%     'CP_tstart','CP_tend','CP_sID','CP_period','dat_sample','datname','crit')
% dat_sample(:,2,:) = dat_sample(:,2,:)*1.03;
% data = cat(3,data,dat_sample(:,1:3,:));


%% plot mean composite

dname = {'Air Temperature','Wind Speed','Absolute Humidity'};
unit = dat_units(1:3);
tref = [-60:1:120]/(24*60);
figure('Position',[0 0 600, 900]);
nrow = 3; ncol = 1;
cl = {'b','r',my_color('dark green')};
Msize = 100;
sub = 0;


pre_max_wspd = nanmax(squeeze(data(1:60,2,:)),[],1);
ind = find(pre_max_wspd >=15);



TC_flag = 1;
for k = 1:length(dname) 
    subplot(nrow,ncol,k); hold on;    


    dat = squeeze(data(:,k,:));
    
    if TC_flag == 1
        dat = dat(:,ind);
    end
    dat_anom = dat - nanmean(dat(40:60,:),1);
    [datmed,CI] = function_median_and_CI_bootstrap_2D(dat_anom,Msize);

    X = tref*24*60;

    plot(X,datmed,'-','Color',cl{k},'LineWidth',1.5);
    f = fill([X';flipud(X')],[CI(:,1);flipud(CI(:,2))],cl{k});
    alpha(f,0.5);    

    text(-30,0,fignum{k},'FontSize',16,'FontWeight','bold')
    title([dname{k}])
    grid on;
    xlabel('time (minute)'); xlim([-20, 60]);
    ylabel([unit{k}])
    set(gca,'Xtick',-20:5:60,'FontSize',12)    
    

end
set(gcf,'Color','w'); pause
if TC_flag == 1
    export_fig([figout,'coldpool_composite_saildrone_TCwind_saildrone.png'],'-dpng','-r200')
else
    export_fig([figout,'coldpool_composite_saildrone.png'],'-dpng','-r200')
end
















































