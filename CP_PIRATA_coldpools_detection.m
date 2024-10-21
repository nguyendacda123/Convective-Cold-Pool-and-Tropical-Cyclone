
% clearvars; close all;
% CP_general_parameters


figout = [figdir,'Cold_pool4/PIRATA/']; 
if ~exist(figout,'dir'); mkdir(figout); end

load([matdir,'PIRATA_data.mat'],'time','Ta','wspd','wdir','sst','rain','sID','sname','lon','lat','rh')
Q = function_humidity_conversion(rh,Ta);
data0 = cat(2,Ta,wspd,Q,sst);
datname = {'Ta','Wspd','Q','SST'};

tref = [-60:1:120]/(24*60); % time window centered at coldpool event

ind = find(lon >= LONR1 & lon <= LONR2 & lat >= LATR1 & lat <= LATR2);

CP_period = [];
CP_tstart = [];
CP_lon = [];
CP_lat = [];
lon_end = [];
lat_start = [];
lat_end = [];
Tslope = [];
CP_tend = [];
CP_sID = [];

dat_sample = [];

plt = 1;
crit = Tcrits(3); % from eparams

dt = 600/86400; % 10min 


for s = ind
    disp(sname{s})
    filt = sID == s;
    tvec = time(filt);

    if sum(round(diff(tvec)*24*60) == 10) == length(tvec) -1
        data = data0(filt,:);
    else
        tvec = nanmin(tvec):dt:nanmax(tvec);
        data = interp1(time(filt),data0(filt,:),tvec);
    end
    
    nt = length(tvec);
    
    % removing nan
    Ta = data(:,1);
    indnan = find(isnan(Ta));
    ind = find(isfinite(Ta));

    if length(ind) > 0 && length(indnan) > 0
        Ta(indnan) = interp1(tvec(ind),Ta(ind),tvec(indnan));
    end    

    
    slope = nan(size(tvec));
    slope(2:end) = Ta(2:end) - Ta(1:end-1);
    Tslope = cat(1,Tslope,slope);

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
                ts(N) = tvec(nt-count) - 120/86400; % account for averaging period  
                
                t3h = ts(N) + tref;
                dat_sample = cat(3,dat_sample,interp1(tvec,data,t3h));
                
                
            end
        else
            if count > 0          
                N = N +1;
                te(N) = tvec(i-1);
                ts(N) = tvec(i-count-1) - 120/86400;
                period(N) = count*10+2;

                t3h = ts(N) + tref;
                dat_sample = cat(3,dat_sample,interp1(tvec,data,t3h));                                
                count = 0;
            end
                
        end
    end
    CP_period = cat(1,CP_period,period');
    CP_tstart = cat(1,CP_tstart,ts');
    CP_tend = cat(1,CP_tend,te');
    CP_sID = cat(1,CP_sID,ones(length(period),1)*s);
 
    
%     %% plot
%     if plt == 1 
%         figure('Position',[0 0 1200, 1200]);
%         nrow = 6; ncol = 1; 
%         if nt > 90*24*6
%             t1 = tvec(end-90*24*6);
%             t2 = tvec(end);
%         else
%             t1 = tvec(1);
%             t2 = tvec(end);
%         end
%         
%         xtick = floor(t1):2:ceil(t2);
%         xtick_minor = floor(t1):1:ceil(t2);
%         tick_label = datestr(xtick);
%         tick_label = tick_label(:,1:6);
%         tlimit = linspace(floor(t1),ceil(t2),nrow+1)      
%         
% 
%         for k = 1:nrow
%             b = subplot(nrow,ncol,k); hold on;
%             h = [];
%             colororder(b,cat(1,my_color('dull light blue'),my_color('black')))
% 
%             yyaxis left
%             h(1) = plot(tvec,data(:,2));
%             h(2) = plot(tvec,data(:,5),'b-');
%             plot(tvec,17*ones(size(tvec)),'--','Color',my_color('dull light blue')) 
%             
%             ylabel('Wspd/rain')  
%             
%             
%             ylim([0, 45])
%             yyaxis right
%             h(3) = plot(tvec,Ta,'LineWidth',1);
%             h(4) = plot(tvec,data(:,4),'m-')
%             Ta1 = interp1(tvec,Ta,ts);
%             Ta2 = interp1(tvec,Ta,te);
%             
%             if ~isempty(Ta1)
%                 h(5) = plot(ts,Ta1,'go','MarkerFaceColor','g');
%                 h(6) = plot(te,Ta2,'rs','MarkerFaceColor','r');
%             end
%             ylabel('T-air (^oC)')
% 
%             if k == 1
%                 legend(h(:),{'Wspd','rain','T-air','sst','start','end'},'Location','best','FontSize',12)
%             end
%             set(gca,'FontSize',12);
%             xlim([tlimit(k),tlimit(k+1)]);
%             set(gca,'Xtick',xtick,'XtickLabel',tick_label,'XMinorTick','on','TickDir','out');
% 
%         end
%         ID = sname{s}; [Y,M] = datevec(tvec(end));
%         axes( 'Position', [0, 0.96, 1, 0.04] ) ;
%         set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
%         text( 0.5, 0, ['PIRATA Buoy ID: ',ID,' (',num2str(Y),')'], 'FontSize', 14', 'FontWeight', 'Bold', ...
%           'HorizontalAlignment', 'Center', 'VerticalAlignment', 'bottom','Color','m' ) ;    
% 
%         set(gcf,'Color','w'); pause
% 
%         export_fig([figout,'PIRATA_',ID,'_individual_plot.png'],'-dpng','-r300')
%         pause(0.1); close all  
%     end   
    
    
end

save([matdir,'PIRATA_buoys_cold_pools.mat'],'lon','lat','sname',...
    'CP_tstart','CP_tend','CP_sID','CP_period','dat_sample','crit','datname')

















