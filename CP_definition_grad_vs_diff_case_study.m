
clearvars; close all;
CP_general_parameters

figout = [figdir,'Cold_pool4/grad_vs_diff/']; 
if ~exist(figout,'dir'); mkdir(figout); end

load([matdir,'Saildrone_rawdata.mat'],'lon','lat','time','date0','Y_label','SD_id','SD_Y',...
                                      'Ta','Rh','Pa','PAR','SWR','Wspd','Wdir','Uwind','Vwind','Gust',...
                                      'SST','SST_skin','SSS','Wav_T','Hs','C_O2','S_O2', ...
                                      'Cspd','Cdir','Chla','CDom','sID')
time = time+ date0;
tref = [-60:1:120]/(24*60); % time window centered at coldpool event

data = cat(2,Ta,Rh,Pa,SWR,Wspd,SST,SSS,Cspd,Chla,CDom);
datname = {'Ta','Rh','Pa','SWR','Wspd','SST','SSS','Cspd','Chla','CDom'};
unit = {'o^C','%','HPa','W.m^-^2','m.s^-^1','^oC','psu','m.s^-^1','kg.m^-^3',' '};


%% Plot example of different methods


SD = unique(sID);
plt  = 0;

filt = mod(SD,1e4) >= 2021;
SD = SD(filt);

datname = {'grad','diff10','diff8'};
nd = length(datname);
dt = 60/86400; % time step

bias = [0, 0.7, 0.79];

tshift = (4*60+30)/86400; % 4min 30 sec
mincount = [5,5,1]; % minimum data to be considered


figure('Position',[0 0 800 400]);
nrow = 1; ncol = 1; sub = 0
casename = {'Hits:','Miss:','False:'};
gr = my_color('dark green');
cl = {gr,'r',my_color('orange')}

version = 0;
Tcrit = -0.68;
for s = 1:length(SD)
    disp(num2str(SD(s)))
    filt = sID == SD(s);
    t_s = time(filt);
    lon_s = lon(filt);
    lat_s = lat(filt);

    data_s = data(filt,:);        
    Ta_s = Ta(filt);

    if length(unique(t_s)) ~= length(t_s)
        [t_s,ind] = unique(t_s);
        Ta_s = Ta_s(ind);
        data_s = data_s(ind,:);
        lon_s = lon_s(ind);
        lat_s = lat_s(ind);
    end    

    % interpolate for missing data

    indnan = find(isnan(Ta_s));
    ind = find(isfinite(Ta_s));

    if length(ind) > 0 && length(indnan) > 0
        Ta_s(indnan) = interp1(t_s(ind),Ta_s(ind),t_s(indnan));
    end

    for k = 1:length(datname)-1
        indnan = find(isnan(data_s(:,k)));
        ind = find(isfinite(data_s(:,k)));
        if length(ind) > 0 && length(indnan) > 0
            data_s(indnan,k) = interp1(t_s(ind),data_s(ind,k),t_s(indnan));
        end
    end        


    tvec = nanmin(t_s)+tshift : dt : nanmax(t_s)-tshift;
    nt = length(tvec);

    % compute Ta slope

    grad10_slope = 0;
    diff10_slope = 0;
    diff2_slope = 0;
    for i = 1:nt          
        if mod(i,10) == 1        
            filt = t_s >= tvec(i) - 5*dt & t_s < tvec(i) + 5*dt;
            if sum(filt) >= mincount(1)
                p = polyfit(t_s(filt),Ta_s(filt),1);
                grad10_slope = p(1)/(60*24)*10;
               
            end

            % 10-min averaging
            filt1 = t_s >= tvec(i) - 10*dt & t_s < tvec(i);
            filt2 = t_s >= tvec(i) & t_s < tvec(i) + 10*dt;

            if sum(filt1) >= mincount(2) && sum(filt2) >= mincount(2)
                diff10_slope = nanmean(Ta_s(filt2)) - nanmean(Ta_s(filt1));             
            end
            
            
            % 2-min averaging
            filt3 = t_s >= tvec(i) - 7*dt & t_s < tvec(i) - 5*dt;
            filt4 = t_s >= tvec(i) + 3*dt & t_s < tvec(i) + 5*dt;
            if sum(filt3) >= mincount(3) && sum(filt4) >= mincount(3)
                diff2_slope = nanmean(Ta_s(filt4)) - nanmean(Ta_s(filt3));             
            end            
            
            
            if sub <= nrow*ncol              
                if (grad10_slope <= Tcrit && diff10_slope <= Tcrit && sub == 0) ||... % hits event
                    (grad10_slope <= Tcrit && diff10_slope > Tcrit  && sub == 1) ||... % miss event
                    (grad10_slope > Tcrit && diff10_slope <= Tcrit  && sub == 2) % false alarm
                
                    filt0 = t_s >= tvec(i) - 15*dt & t_s <= tvec(i) + 15*dt;
                    Ta_case = Ta_s(filt0); % case study
                    t_case = t_s(filt0);
                    
                    Ta_point = interp1(t_case,Ta_case,tvec(i));
                    
                    t_grad10 = [tvec(i)-5*dt, tvec(i)+ 5*dt];
                    Ta_grad10 = p(1)*t_grad10 + p(2); 
                    
                    t_diff10 = [tvec(i), tvec(i)+10*dt];
                    Ta_diff10 = [nanmean(Ta_s(filt1)),nanmean(Ta_s(filt2))]
                    
                    t_diff2 = [tvec(i)-5*dt, tvec(i)+5*dt];
                    Ta_diff2 = [nanmean(Ta_s(filt3)),nanmean(Ta_s(filt4))]                    
                    
                    sub = sub+1; subplot(nrow,ncol,sub); hold on;

                    
                    plot(t_case,Ta_case,'k-o','LineWidth',1.5); 
                    h(1) = plot(t_grad10,Ta_grad10,'b--o','MarkerFaceColor','b','LineWidth',1);
                    
                    plot([t_grad10(1)-5*dt,t_grad10(1)+5*dt],ones(2,1)*Ta_diff10(1),'-r','LineWidth',1.5);
                    plot([t_grad10(2)-5*dt,t_grad10(2)+5*dt],ones(2,1)*Ta_diff10(2),'-r','LineWidth',1.5);
                    h(2) = plot(t_diff10,Ta_diff10, 'r--o','MarkerFaceColor','r','LineWidth',1.5)
                    
                    
                    
                    
                    
                    plot([t_diff2(1),t_diff2(1)-2*dt],ones(2,1)*Ta_diff2(1),'-','LineWidth',1.5,'Color',gr);
                    plot([t_diff2(2),t_diff2(2)-2*dt],ones(2,1)*Ta_diff2(2),'-','LineWidth',1.5,'Color',gr);
                    h(3) = plot(t_diff2,Ta_diff2, '--o','Color',gr,'MarkerFaceColor',gr,'LineWidth',1.5)                    
                    
                    plot(tvec(i),Ta_point,'m^','MarkerFaceColor','m','MarKerSize',5);
                    plot([tvec(i),tvec(i)],[Ta_point-0.3, Ta_point + 0.3],'m-');
                    text(tvec(i),Ta_point+0.3,'t0','Color','m','FontSize',14,'FontWeight','bold'); pause
                    
                    % plot(ones(2,1)*tvec(i),Ta_point+[-1,1],'m-','LineWidth',1.5)
                    
                    xlab = datestr(t_case(1:5:end));
                    ID = num2str(SD(s));                    

                    set(gca,'Xtick',t_case(1:5:end),'XtickLabel',xlab(:,13:17))
                    xlim([min(t_case),max(t_case)])
                    xlabel('time'); ylabel('Ta (^oC)'); grid on;
                    if sub == 1
                        legend(h(:),{'dTa_S_D','dTa_N_D_B_C','dTa_P_I_R_A_T_A'},'Location','best','FontSize',12)
                    end
                    title(['SD',ID(1:4),' (',xlab(1,1:11),')'])                              
                    set(gca,'FontSize',12,'Color',0.95*[1,1,1]);                    

                    set(gcf,'Color','w'); 
                   
                    export_fig([figout,'CP_magnitude_diff_vs_grad_v',num2str(version),'.png'],'-dpng','-r300') 
                    version = version + 1; sub = 0;
                    pause(0.1)   
                    clf
                        
                end 

            end            
        
        end
 
    end

end


%% Plot hits, missings, false alarms cases


SD = unique(sID);
plt  = 0;


filt = mod(SD,1e4) >= 2021;
SD = SD(filt);




CP_period = cell(3,1);
CP_tstart = cell(3,1);
CP_lon = cell(3,1);
CP_lat = cell(3,1);
CP_tend = cell(3,1);
CP_sID = cell(3,1);
Tslope = [];
time_slope = [];


datname = {'grad','diff10','diff8'};
nd = length(datname);
dt = 60/86400; % time step

bias = [0, 0.7, 0.79];

tshift = (4*60+30)/86400; % 4min 30 sec
mincount = [5,5,1]; % minimum data to be considered


figure('Position',[0 0 800 1200]);
nrow = 3; ncol = 1; sub = 0
casename = {'Hits:','Miss:','False:'};
gr = my_color('dark green');
cl = {gr,'r',my_color('orange')}

version = 0;
Tcrit = -0.68;
for s = 1:length(SD)
    disp(num2str(SD(s)))
    filt = sID == SD(s);
    t_s = time(filt);
    lon_s = lon(filt);
    lat_s = lat(filt);

    data_s = data(filt,:);        
    Ta_s = Ta(filt);

    if length(unique(t_s)) ~= length(t_s)
        [t_s,ind] = unique(t_s);
        Ta_s = Ta_s(ind);
        data_s = data_s(ind,:);
        lon_s = lon_s(ind);
        lat_s = lat_s(ind);
    end    

    % interpolate for missing data

    indnan = find(isnan(Ta_s));
    ind = find(isfinite(Ta_s));

    if length(ind) > 0 && length(indnan) > 0
        Ta_s(indnan) = interp1(t_s(ind),Ta_s(ind),t_s(indnan));
    end

    for k = 1:length(datname)-1
        indnan = find(isnan(data_s(:,k)));
        ind = find(isfinite(data_s(:,k)));
        if length(ind) > 0 && length(indnan) > 0
            data_s(indnan,k) = interp1(t_s(ind),data_s(ind,k),t_s(indnan));
        end
    end        


    tvec = nanmin(t_s)+tshift : dt : nanmax(t_s)-tshift;
    nt = length(tvec);

    % compute Ta slope

    grad10_slope = 0;
    diff10_slope = 0;
    diff2_slope = 0;
    for i = 1:nt          
        if mod(i,10) == 1        
            filt = t_s >= tvec(i) - 5*dt & t_s < tvec(i) + 5*dt;
            if sum(filt) >= mincount(1)
                p = polyfit(t_s(filt),Ta_s(filt),1);
                grad10_slope = p(1)/(60*24)*10;
               
            end

            % 10-min averaging
            filt1 = t_s >= tvec(i) - 10*dt & t_s < tvec(i);
            filt2 = t_s >= tvec(i) & t_s < tvec(i) + 10*dt;

            if sum(filt1) >= mincount(2) && sum(filt2) >= mincount(2)
                diff10_slope = nanmean(Ta_s(filt2)) - nanmean(Ta_s(filt1));             
            end
            
            
            % 2-min averaging
            filt3 = t_s >= tvec(i) - 6*dt & t_s < tvec(i) - 4*dt;
            filt4 = t_s >= tvec(i) + 4*dt & t_s < tvec(i) + 6*dt;
            if sum(filt3) >= mincount(3) && sum(filt4) >= mincount(3)
                diff2_slope = nanmean(Ta_s(filt4)) - nanmean(Ta_s(filt3));             
            end            
            
            
            if sub <= nrow*ncol              
                if (grad10_slope <= Tcrit && diff10_slope <= Tcrit && sub == 0) ||... % hits event
                    (grad10_slope <= Tcrit && diff10_slope > Tcrit  && sub == 1) ||... % miss event
                    (grad10_slope > Tcrit && diff10_slope <= Tcrit  && sub == 2) % false alarm
                
                    filt0 = t_s >= tvec(i) - 15*dt & t_s <= tvec(i) + 15*dt;
                    Ta_case = Ta_s(filt0); % case study
                    t_case = t_s(filt0);
                    
                    Ta_point = interp1(t_case,Ta_case,tvec(i));
                    
                    t_grad10 = [tvec(i)-5*dt, tvec(i)+ 5*dt];
                    Ta_grad10 = p(1)*t_grad10 + p(2); 
                    
                    t_diff10 = [tvec(i)-5*dt, tvec(i)+5*dt];
                    Ta_diff10 = [nanmean(Ta_s(filt1)),nanmean(Ta_s(filt2))]
                    
                    t_diff2 = [tvec(i)-5*dt, tvec(i)+5*dt];
                    Ta_diff2 = [nanmean(Ta_s(filt3)),nanmean(Ta_s(filt4))]                    
                    
                    sub = sub+1; subplot(nrow,ncol,sub); hold on;

                    
                    plot(t_case,Ta_case,'k-o','LineWidth',1.5); 
                    h(1) = plot(t_grad10,Ta_grad10,'b--o','MarkerFaceColor','b','LineWidth',1);
                    
                    plot(t_s(filt1),ones(sum(filt1),1)*Ta_diff10(1),'-r','LineWidth',1.5);
                    plot(t_s(filt2),ones(sum(filt2),1)*Ta_diff10(2),'-r','LineWidth',1.5);
                    h(2) = plot(t_diff10,Ta_diff10, 'r--o','MarkerFaceColor','r','LineWidth',1.5)

                    plot(t_s(filt3),ones(sum(filt3),1)*Ta_diff2(1),'-','LineWidth',1.5,'Color',gr);
                    plot(t_s(filt4),ones(sum(filt4),1)*Ta_diff2(2),'-','LineWidth',1.5,'Color',gr);
                    h(3) = plot(t_diff2,Ta_diff2, '--o','Color',gr,'MarkerFaceColor',gr,'LineWidth',1.5)                    
                    
                    plot(tvec(i),Ta_point,'m^','MarkerFaceColor','m','MarKerSize',10);
                    % plot(ones(2,1)*tvec(i),Ta_point+[-1,1],'m-','LineWidth',1.5)
                    
                    xlab = datestr(t_case(1:5:end));
                    ID = num2str(SD(s));                    

                    set(gca,'Xtick',t_case(1:5:end),'XtickLabel',xlab(:,13:17))
                    xlim([min(t_case),max(t_case)])
                    xlabel('time'); ylabel('Ta (^oC)'); grid on;
                    if sub == 1
                        legend(h(:),{'grad-10','diff-mean10','diff-mean2'},'Location','best')
                    end
                    title(['SD',ID(1:4),' (',xlab(1,1:11),')'])
                    text(t_case(2),nanmin(Ta_case)+0.5,{[casename{sub}];
                        ['grad-10 = ',num2str(grad10_slope,'%.2f')];
                        ['diff-mean10 = ',num2str(diff10_slope,'%.2f')];
                        ['diff-mean2 = ',num2str(diff2_slope,'%.2f')]},'FontSize',12,'Color',cl{sub});              
                    
                    set(gca,'FontSize',12,'Color',0.95*[1,1,1]);                    
                end
                if sub == 3
                    set(gcf,'Color','w'); pause
                   
                    export_fig([figout,'CP_magnitude_diff_vs_grad_v',num2str(version),'.png'],'-dpng','-r200') 
                    version = version + 1; sub = 0;
                    version
                    pause(0.1)   
                    clf
                        
                end 

            end            
        
        end
 
    end

end


save([matdir,'Cold_pools_grad_vs_diff_2021_2023_biascorr.mat'],'Tslope','time_slope',...
    'CP_period','CP_tstart','CP_tend','CP_sID','CP_lon','CP_lat')




%% plot hits missing false alarm


load([matdir,'Cold_pools_grad_vs_diff_2021_2023_biascorr.mat'],'Tslope','time_slope',...
    'CP_period','CP_tstart','CP_tend','CP_sID','CP_lon','CP_lat')

filt = sID == 10312021;
t_s = time(filt);
lon_s = lon(filt);
lat_s = lat(filt);
data_s = data(filt,:);        
Ta_s = Ta(filt);

if length(unique(t_s)) ~= length(t_s)
    [t_s,ind] = unique(t_s);
    Ta_s = Ta_s(ind);
    data_s = data_s(ind,:);
    lon_s = lon_s(ind);
    lat_s = lat_s(ind);
end 

% interpolate for missing data

indnan = find(isnan(Ta_s));
ind = find(isfinite(Ta_s));

if length(ind) > 0 && length(indnan) > 0
    Ta_s(indnan) = interp1(t_s(ind),Ta_s(ind),t_s(indnan));
end

for k = 1:length(datname)-1
    indnan = find(isnan(data_s(:,k)));
    ind = find(isfinite(data_s(:,k)));
    if length(ind) > 0 && length(indnan) > 0
        data_s(indnan,k) = interp1(t_s(ind),data_s(ind,k),t_s(indnan));
    end
end 


figure('Position',[0 0 1200, 1200]);
nrow = 6; ncol = 1;

xtick = floor(t_s(1)):1/48:ceil(t_s(end));
tick_label = datestr(xtick);
tick_label = tick_label(:,13:17);
tlimit = linspace(floor(t_s(1)+19),floor(t_s(1))+21,nrow+1);      

for k = 1:nrow
    b = subplot(nrow,ncol,k); hold on;
    h = [];
    colororder(b,cat(1,my_color('dull light blue'),my_color('black')))

    yyaxis left
    h(1) = plot(t_s,data_s(:,5));
    plot(t_s,17*ones(size(t_s)),'--','Color',my_color('dull light blue'))
    ylabel('Wspd (m.s^-^1)')  
    ylim([0, 45])
    yyaxis right
    h(2) = plot(t_s,Ta_s,'LineWidth',1);
    msize = [6,6,6]; ylim([25.5, 29])
    cl = {'g','r'}
 
    for i = 1:2
        filt2 = CP_sID{i} == 10312021;
        ts = CP_tstart{i}(filt2);
        te = CP_tend{i}(filt2);
        Ta1 = interp1(t_s,Ta_s,ts);
        Ta2 = interp1(t_s,Ta_s,te);
        if ~isempty(Ta1)
            h(3) = plot(ts,Ta1,'o','MarkerFaceColor',cl{i},'MarkerSize',msize(i));
            h(4) = plot(te,Ta2,'s','MarkerFaceColor',cl{i},'MarkerSize',msize(i));
        end
    end
    ylabel('T-air (^oC)')

%     if k == 1
%         legend(h(:),{'Wspd','T-air','start','end'},'Location','best','FontSize',12)
%     end
    set(gca,'FontSize',12);
    xlim([tlimit(k),tlimit(k+1)]);
    set(gca,'Xtick',xtick,'XtickLabel',tick_label,'TickDir','out');

end
ID = num2str(SD(s));
axes( 'Position', [0, 0.96, 1, 0.04] ) ;
set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
text( 0.5, 0, ['SD',ID(1:4),' - ',ID(5:end)], 'FontSize', 14', 'FontWeight', 'Bold', ...
  'HorizontalAlignment', 'Center', 'VerticalAlignment', 'bottom','Color','m' ) ;    

% set(gcf,'Color','w'); pause
% export_fig([figout,'SD',ID,'_individual_plot.png'],'-dpng','-r300')
% pause(0.1); close all  
 
    
         
%%

        
    


    







































