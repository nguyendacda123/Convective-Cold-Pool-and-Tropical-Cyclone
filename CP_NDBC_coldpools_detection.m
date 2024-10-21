
% clearvars; close all;
% CP_general_parameters

datadir = [orig_dir, 'NDBC_buoys/'];
figout = [figdir,'Cold_pool4/NDBC/'];
if ~exist(figout,'dir'); mkdir(figout); end



% load([matdir,'NDBC_offshore_buoys_metadata.mat'],'lons_10','lats_10','sname_10')
load([matdir,'NDBC_offshore_buoys_raw_data.mat'],'N_gtime','N_gTa','N_RH',...
    'N_gID','N_glon','N_glat','N_SST','N_WSPD','lons_10','lats_10','sname_10')

tref = [-60:1:120]/(24*60); % time window centered at coldpool event


CP_period = [];
CP_tstart = [];
CP_lon = [];
CP_lat = [];
lon_end = [];
lat_start = [];
lat_end = [];

CP_tend = [];
CP_sID = [];

dat_sample = [];

plt = 0;
crit = Tcrits(2); % from eparams
dt = 600/86400; 
ref_date = datenum(1970,1,1);

for s = 1:length(sname_10)
    filt = N_gID == s;
       
    fname = [datadir, sname_10{s},'h9999.nc']

    Ta = N_gTa(filt);

    rh = N_RH(filt);
    Q = function_humidity_conversion(rh,Ta);
    
    sst = N_SST(filt);
    wspd = N_WSPD(filt);

    time = N_gtime(filt);
    tvec = nanmin(time):dt:nanmax(time);

    nt = length(tvec);

    data0 = cat(2,Ta,wspd,Q,sst);
    datname = {'Ta','Wspd','Q','SST'};    
    
    data = interp1(time,data0,tvec);
    
    % removing nan
    Ta = data(:,1);
    indnan = find(isnan(Ta));
    ind = find(isfinite(Ta));

    if length(ind) > 0 && length(indnan) > 0
        Ta(indnan) = interp1(tvec(ind),Ta(ind),tvec(indnan));
    end    

    
    slope = nan(size(tvec));
    slope(2:end) = Ta(2:end) - Ta(1:end-1);
    

    flag = zeros(size(slope));
    flag(slope <= crit) = 1;
    
    % compute number of cold pool events
    i = 0;
    count = 0;
    N = 0;
    ts = [];
    te = [];
    period = [];

    for i = 1:nt
        if flag(i) == 1
            if i < nt
                count = count + 1;
            else % end of series
                count = count + 1;
                N = N +1;
                te(N) = tvec(nt);
                ts(N) = tvec(nt-count)-600/86400; % to account for averaging period 
                
                t3h = ts(N) + tref;
                dat_sample = cat(3,dat_sample,interp1(tvec,data,t3h));     
            end
        else
            if count > 0          
                N = N +1;
                te(N) = tvec(i-1);
                ts(N) = tvec(i-count-1)-600/86400;  % to account for averaging period 
                period(N) = (count+1)*10;

                t3h = ts(N) + tref;
                dat_sample = cat(3,dat_sample,interp1(tvec,data,t3h));                                
                count = 0;
            end
                
        end
    end
    CP_period = cat(1,CP_period,period');
    CP_tstart = cat(1,CP_tstart,ts');
    CP_tend = cat(1,CP_tend,te');
    CP_sID = cat(1,CP_sID,ones(N,1)*s);
    
    
    
    [length(CP_tstart), length(CP_sID)]
    %% plot
    if plt == 1 
        figure('Position',[0 0 1200, 1200]);
        nrow = 6; ncol = 1;      
        xtick = floor(tvec(1)):5:ceil(tvec(end));
        xtick_minor = floor(tvec(1)):1:ceil(tvec(end));
        tick_label = datestr(xtick);
        tick_label = tick_label(:,1:6);
        tlimit = linspace(floor(tvec(1)),ceil(tvec(end)),nrow+1)      
        

        for k = 1:nrow-3
            b = subplot(nrow-3,ncol,k); hold on;
            h = [];
            colororder(b,cat(1,my_color('dull light blue'),my_color('black')))

            yyaxis left
            h(1) = plot(tvec,data(:,3));
            plot(tvec,17*ones(size(tvec)),'--','Color',my_color('dull light blue'))
            ylabel('Wspd (m.s^-^1)')  
            ylim([0, 45])
            yyaxis right
            h(2) = plot(tvec,Ta,'LineWidth',1);
            h(3) = plot(tvec,data(:,5),'m-')
            Ta1 = interp1(tvec,Ta,ts);
            Ta2 = interp1(tvec,Ta,te);
            
            if ~isempty(Ta1)
                h(4) = plot(ts,Ta1,'go','MarkerFaceColor','g');
                h(5) = plot(te,Ta2,'rs','MarkerFaceColor','r');
            end
            ylabel('T-air (^oC)')

            if k == 1
                legend(h(:),{'Wspd','T-air','sst','start','end'},'Location','best','FontSize',12)
            end
            set(gca,'FontSize',12);
            xlim([tlimit(k),tlimit(k+1)]);
            set(gca,'Xtick',xtick,'XtickLabel',tick_label,'TickDir','out');

        end
        ID = sname_10{s};
        axes( 'Position', [0, 0.96, 1, 0.04] ) ;
        set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
        text( 0.5, 0, ['NDBC Buoy ID: ',ID], 'FontSize', 14', 'FontWeight', 'Bold', ...
          'HorizontalAlignment', 'Center', 'VerticalAlignment', 'bottom','Color','m' ) ;    

        set(gcf,'Color','w'); pause
%         close all;
%         export_fig([figout,'SD',ID,'_individual_plot.png'],'-dpng','-r300')
%         pause(0.1); close all  
    end   
    
    
end

wlevel_obs_payload = {'bara9','chsv3','frcb6','kywf1','lamv3','ltbv3','misp4','vcaf1'};
carlcoos_payload = {'41052'};
scoop_payload = setdiff(sname_10,[wlevel_obs_payload,carlcoos_payload]);



save([matdir,'NDBC_buoys_cold_pools_offshore',tag,'.mat'],'lons_10','lats_10','sname_10',...
    'CP_tstart','CP_tend','CP_sID','CP_period','dat_sample',...
    'datname','wlevel_obs_payload','carlcoos_payload','scoop_payload','crit')













