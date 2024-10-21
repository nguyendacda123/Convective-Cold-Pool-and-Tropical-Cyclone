
clearvars; close all;
CP_general_parameters


figout = [figdir,'Cold_pool_R1/definition/']; 
if ~exist(figout,'dir'); mkdir(figout); end
load([matdir,'Saildrone_rawdata.mat'],'lon','lat','sID')

load([matdir,'Saildrone_cold_pools_2018_2023_3Tcrits.mat'],'CP_period','CP_tstart','CP_tend','CP_sID','Tslope',...
    'CP_lon','CP_lat','dat_sample','dat_units','datname','method','diff2_mis','diff2_far','diff10_mis','diff10_far')


%% plot histogram
[nm,nrep] = size(CP_period);
tbin = 5:10:75;
nbin = length(tbin);

mean_count = nan(nm,nbin);

CI = nan(nm,nbin,2);

tot_count = nan(nm,1);
CI_tot = nan(nm,2);

GF_time = nan(nm,1); % gust front time
CI_gust = nan(nm,2);

tab1 = [];
for d = 1:nm
    counts = nan(nbin,nrep);
    GF = nan(nrep,1);
    for k = 1:nrep
        dat = round(CP_period{d,k});
        counts(:,k) = histc(dat,tbin);
        GF(k) = nansum(dat)/60/24; % day
    end
    [mean_count(d,:), CI(d,:,:)] = function_median_and_CI_bootstrap_2D(counts,100);
    [tot_count(d), CI_tot(d,:)] = function_median_and_CI_bootstrap_2D(sum(counts,1),100);
    [GF_time(d), CI_gust(d,:)] = function_median_and_CI_bootstrap(GF,100);
    if d > 1
        row = cat(2,tot_count(d),CI_tot(d,:));
        tab1 = cat(2,tab1,row);
    end
end
mean_count(mean_count == 0) = nan;


% compute missing, far rates
Misdat = {diff10_mis,diff2_mis};
Fardat = {diff10_far,diff2_far};
tab2 = [];
for d = 1:2
    col = [];
    for s = 1:4
        if s < 4
            dat = Misdat{d};
            val = s;
        else
            dat = Fardat{d};
            val = 0;
        end
        N = size(dat,1);
        [meand ,cid] = function_median_and_CI_bootstrap_2D(sum(dat == val,1)/N*100,100);
        row = cat(2,meand,cid);
        col = cat(1,col,row);
    end
    tab2 = cat(2,tab2,col);
end



figure('Position',[0, 0, 600,900]); hold on;
nrow = 3; ncol = 1;
xaxis = 10:10:80;
for d = 1:nm
      
    subplot(nrow,ncol,d); hold on;
    h(1) = bar(xaxis,mean_count(d,:)); 
    if d > 1
        errorbar(xaxis,mean_count(d,:),squeeze(CI(d,:,1)-mean_count(d,:)),squeeze(CI(d,:,2)-mean_count(d,:)),'k.');
    end
    Nevents = num2str(tot_count(d),'%.0f');
    if d == 1
        text(20,50,{['Number of CPs: ',Nevents];...
                    ['CP time (day): ',num2str(GF_time(d),'%.1f')]},'Color','b','FontSize',12)
    else
        lower = num2str(round(CI_tot(d,1)));
        upper = num2str(round(CI_tot(d,2)));
        lower2 = num2str(CI_gust(d,1),'%.1f');
        upper2 = num2str(CI_gust(d,2),'%.1f');
                
        
        text(20,50,{['Number of CPs: ',Nevents,' (',lower,', ',upper,')'];...
            ['CP time (day): ',num2str(GF_time(d),'%.1f'),' (',lower2,', ',upper2,')']},'Color','b','FontSize',12)
    end
    
    xlabel('CP period (minutes)'); ylabel('count')
    xlim([5, 85]); text(7,1500,fignum{d},'FontSize',16,'FontWeight','bold')
    title(method{d})        
    set(gca,'FontSize',12);

end
set(gcf,'Color','w'); pause
export_fig([figout,'Saildrone_CP_periods_hist.png'],'-dpng','-r300') 




%% plot maps of locations
method = {'Saildrone','Simulated NDBC','Simulated PIRATA'};
figure('Position',[0 0 800 1200]);
nrow = 3; ncol = 1;

for d = 1:3
    subplot(nrow,ncol,d); hold on;
    plot_projection_ROI(260,315,12,36);

    SD = unique(CP_sID{d,1});
    for s = 1:length(SD)
        filt = sID == SD(s);
        lon_s = lon(filt);
        lat_s = lat(filt);

        lon_s(lon_s < 0) = lon_s(lon_s < 0) + 360;

        ID = num2str(SD(s));
    %     if ID(5:end) == '2023'
        m_plot(lon_s,lat_s,'m-')
    %     m_text(nanmax(lon_s),nanmin(lat_s),[ID(1:4),'-',ID(7:8)],'FontSize',12)
    %     end



    end

    CPlon = CP_lon{d,1};
    CPlat = CP_lat{d,1};
    CPperiod = CP_period{d,1};

    CPlon(CPlon < 0) = CPlon(CPlon < 0) + 360;
    set(gca,'FontSize',12);
    m_scatter(CPlon,CPlat,15,CPperiod,'filled'); colormap(jet); cb = colorbar; caxis([10, 40])
    m_text(CPlon(1),CPlat(1),fignum{d},'FontSize',16,'FontWeight','bold')
    xlabel(cb,'Cold-pool period (min)','FontSize',12)
    title([method{d},' coldpools'],'FontSize',16)
end
set(gcf,'Color','w'); pause
export_fig([figout,'Colpool_locations_saildrone.png'],'-dpng','-r200')


%% plot mean composite

ind = [1,2,3]; % 
% dname = datname(ind);
cname = {'All conditions','TC conditions'};
dname = {'Ta (^oC)','Wspd (m.s^-^1)','AH (g.m^-^3)'};
unit = dat_units(ind);
tref = [-60:1:120]/(24*60);
figure('Position',[0 0 1000, 900]);
nrow = 3; ncol = 2;
cl = {'b','r',my_color('dark green')};
Msize = 1000;
sub = 0;
order = [1,3,5,2,4,6];
method = {'Saildrone','NDBC','PIRATA'};

y_lim = {[-1.8, 0.5],[-0.5,2.5],[-0.7,0.3]};

for TC_flag = 0:1
    for k = 1:length(dname) 
        sub=sub+1; subplot(nrow,ncol,order(sub)); hold on;    
        h = [];
        for m = 1:length(method)
            dat = squeeze(dat_sample{m}(:,ind(k),:));  
            if TC_flag == 1
                max_wspd = nanmax(dat_sample{m}(:,2,:),[],1);
                index = find(max_wspd >=15);
                dat = dat(:,index);
                length(index)
            end
            dat_anom = dat - nanmean(dat(1:30,:),1);
            [datmed,CI] = function_median_and_CI_bootstrap_2D(dat_anom,Msize);

            X = tref*24*60;

            h(m) = plot(X,datmed,'-','Color',cl{m},'LineWidth',1.5);
            f = fill([X';flipud(X')],[CI(:,1);flipud(CI(:,2))],cl{m});
            alpha(f,0.5);    
        end
        text(-30,0,fignum{order(sub)},'FontSize',16,'FontWeight','bold')
        if k == 1
            title([cname{TC_flag+1}])
        end
        grid on;
        xlabel('time (minute)'); xlim([-30, 60]); ylim(y_lim{k})
        if TC_flag == 0
            ylabel([dname{k}])
        end
        set(gca,'Xtick',-30:10:60,'FontSize',12)    

        if k == 1
            legend(h(:),method,'Location','best')
        end
    end
end
set(gcf,'Color','w'); pause
if TC_flag == 1
    export_fig([figout,'coldpool_composite_saildrone_TCwind.png'],'-dpng','-r200')
else
    export_fig([figout,'coldpool_composite_saildrone.png'],'-dpng','-r200')
end


%% saildrone data days
load([matdir,'Saildrone_rawdata.mat'],'sID','time','date0');
time = time + date0;


SD_list = unique(sID);
Ylabel = mod(SD_list,1e4);
filt = Ylabel >= 2021;
SD_list = SD_list(filt);

Ymin = 2021;
Ymax = 2023;
ny = Ymax-Ymin+1;
Ndays = nan(ny,12);
Ns = nan(ny,12);

for i = 1:ny
    Y = Ymin + i -1;
    for M = 1:12
        filt = time >= datenum(Y,M,1) & time <= datenum(Y,M,nday(Y,M),23,59,59);
        t = time(filt);
        dt = abs(t(2:end) - t(1:end-1));
        dt(dt > 600/86400) = nan;
        
        Ndays(i,M) = nansum(dt);
        Ns(i,M) = length(unique(sID(filt)));
    end
end
Ndays_M = nansum(Ndays,1)';
Ndays_Y = nansum(Ndays,2);



tlabel = {};
for Y = Ymin:Ymax
    tlabel = [tlabel,num2str(Y),'Jul'];
end

Ndays(Ndays == 0) = nan;
Ns(Ns == 0) = nan;

xaxis = 1:ny*12;
figure('Position',[0 0 800 400]); hold on;
yyaxis left
plot(xaxis,reshape(Ndays',[],1),'-o')
ylabel('days of data per month')

yyaxis right
plot(xaxis,reshape(Ns',[],1),'-s')
ylabel('number of active stations')
set(gca,'Xtick',xaxis(1:6:end),'XtickLabel',tlabel,'FontSize',12);
xlabel('time')
xlim([min(xaxis),max(xaxis)])
legend('data days','active stations');
title('Saildrone data coverage')
set(gcf,'Color','w');

export_fig([figout,'Saildrone_data_days_and_No_stations.png'],'-dpng','-r300') 

%%

Mbin = 1:12;
Ybin = 2021:2023;

nbin1 = length(Mbin);
nbin2 = length(Ybin);

M_count = nan(nm,nbin1);
Y_count = nan(nm,nbin2);
M_CI = nan(nm,nbin1,2);
Y_CI = nan(nm,nbin2,2);
for d = 1:nm
    count_M = nan(nbin1,nrep);
    count_Y = nan(nbin2,nrep);
    for k = 1:nrep
        ID = CP_sID{d,k};
        tstart = CP_tstart{d,k};


        N = length(tstart);
        Yc = nan(N,1);
        Mc = nan(N,1);
        Dc = nan(N,1);
        for i = 1:length(tstart)
            [Yc(i),Mc(i),Dc(i)] = datevec(tstart(i));
        end
        for M = 1:12
            count_M(M,k) = nansum(Mc == M);
        end
        for Y = 2021:2023
            count_Y(Y-Ymin+1,k) = nansum(Yc == Y);
        end
        
        count_M(:,k) = count_M(:,k)./Ndays_M;
        count_Y(:,k) = count_Y(:,k)./Ndays_Y;
    end

    [M_count(d,:), M_CI(d,:,:)] = function_median_and_CI_bootstrap_2D(count_M,100);
    [Y_count(d,:), Y_CI(d,:,:)] = function_median_and_CI_bootstrap_2D(count_Y,100);       
        
end
M_count(~isfinite(M_count)) = nan;
Y_count(~isfinite(Y_count)) = nan;
%%
figure('Position',[0 0 1000 800]);
nrow = 3; ncol = 2;
sub = 0;
for d = 1:length(method)
    
    % seasonal --------------------
    sub = sub+1; subplot(nrow,ncol,sub); hold on;
    h(1) = bar(Mbin,M_count(d,:)); 
    if d > 1
        errorbar(Mbin,M_count(d,:),squeeze(M_CI(d,:,1)-M_count(d,:)),squeeze(M_CI(d,:,2)-M_count(d,:)),'k.');
    end
    
    xlabel('month'); ylabel('CP.day^-^1')
    xlim([0.5 12.5]);
    title(method{d}) 
    set(gca,'Xtick',Mbin,'FontSize',12);   

    
    
    % interann ------------------
    sub = sub+1; subplot(nrow,ncol,sub); hold on;
    h(1) = bar(Ybin,Y_count(d,:)); 
    if d > 1
        errorbar(Ybin,Y_count(d,:),squeeze(Y_CI(d,:,1)-Y_count(d,:)),squeeze(Y_CI(d,:,2)-Y_count(d,:)),'k.');
    end

    xlabel('Year'); ylabel('CP.day^-^1')
    xlim([2020 2024]); 
    title(method{d})        
    set(gca,'Xtick',Ybin,'FontSize',12);    
    
end
axes( 'Position', [0, 0.96, 1, 0.04] ) ;
set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
text( 0.5, 0, ['Cold pool frequencies (scaled by number of saildrones)'], 'FontSize', 16', 'FontWeight', 'Bold', ...
  'HorizontalAlignment', 'Center', 'VerticalAlignment', 'bottom','Color','k' ) ; 

set(gcf,'Color','w'); pause
export_fig([figout,'Coldpools_seasonal_interann_saildrone',tag,'.png'],'-dpng','-r200') 













































