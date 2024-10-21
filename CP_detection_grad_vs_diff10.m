
clearvars; close all;
CP_general_parameters

figout = [figdir,'Cold_pool4/definition/']; 
if ~exist(figout,'dir'); mkdir(figout); end


load([matdir,'Saildrone_rawdata.mat'],'time','date0','sID','Ta','Wspd')
time = time + date0;
load([matdir,'Saildrone_cold_pools_2021_2023_3Tcrits.mat'],'CP_period','CP_tstart','CP_tend','CP_sID','Tslope',...
    'CP_lon','CP_lat','dat_sample','dat_units','datname','method')

load([matdir,'Cold_pools_grad_vs_diff_2018_2023.mat'],'TW_slope','time_slope','method','datname','slope_sID')
Tslope = squeeze(TW_slope(:,1,:));




%% plot missing far events
mincount = [5,5,1]; % minimum data to be considered
dt = 60/86400;

gr = my_color('dark green');
figure('Position',[0 0 800, 500]); hold on;
version = 0;

plt = 1;

t11 = CP_tstart{1,1};
t12 = CP_tend{1,1};
ID1 = CP_sID{1,1};

diff10_mis = zeros(length(t11),9);

for j = 1:9
    j
    t21 = CP_tstart{2,j};
    t22 = CP_tend{2,j};
    ID2 = CP_sID{2,j};


    for k = 1:length(t11)
        


        filt_s = sID == ID1(k);
        t_s = time(filt_s);
        wspd_s = Wspd(filt_s);
        Ta_s = Ta(filt_s);

        if length(unique(t_s)) ~= length(t_s)
            [t_s,ind] = unique(t_s);
            wspd_s = wspd_s(ind,:);
            Ta_s = Ta_s(ind);
        end          


        % filter slope data
        filt_sl = slope_sID == ID1(k);
        t_grad = time_slope(filt_sl);
        slope_grad = Tslope(filt_sl,1);


        t_diff1 = min(t_grad):dt*0.1:max(t_grad); % highres diff slope
        slope_diff1 = interp1(time_slope(filt_sl),Tslope(filt_sl,2),t_diff1); 
        
        t_diff2 = t_grad(j:10:end); % actual res diff slope           
        temp = Tslope(filt_sl,2);
        slope_diff2 = temp(1:10:end);           

        filt_case = t_s >= t11(k) - 15*dt & t_s <= t12(k) + 15*dt;       
        Ta_case = Ta_s(filt_case); % case study
        t_case = t_s(filt_case);

        filt0 = t_grad >= t11(k) - 15*dt & t_grad <= t12(k) + 15*dt; 
        filt1 = t_diff1 >= t11(k) - 15*dt & t_diff1 <= t12(k) + 15*dt;   
        filt2 = t_diff2 >= t11(k) - 15*dt & t_diff2 <= t12(k) + 15*dt; 
             
        
        filt = ((t21 >=t11(k) & t21 <= t12(k)) | ... % 1-2-1
                (t22 >= t11(k) & t22 <= t12(k)) |... % 1-2-1
                (t21 <= t11(k) & t22 >= t12(k))) ... % 2-1-2
                & ID2 == ID1(k);
            
        if sum(filt) > 0
            diff10_mis(k,j) = 1;
        else
            dslope1 = slope_diff1(filt1);

            if sum(dslope1 <= Tcrits(2)) == 0
                diff10_mis(k,j) = 2;
            else
                diff10_mis(k,j) = 3;
            end
        end



        if diff10_mis(k,1) == 2 && plt == 1

    %         sub = sub+1; subplot(nrow,ncol,sub); hold on;
            h = [];
            xlab = datestr(t_case(1:5:end));
            yyaxis left;

            h(1) = plot(t_case,Ta_case,'-o','LineWidth',1.5); 
            Tmin = floor(nanmin(Ta_case));
            Tmax = ceil(nanmax(Ta_case));
            set(gca,'Xtick',t_case(1:5:end),'XtickLabel',xlab(:,13:17),'Ytick',Tmin:0.5:Tmax)
            xlim([min(t_case),max(t_case)])
            xlabel('time'); ylabel('Ta (^oC)'); grid on;
            ylim([Tmin,Tmax])

            yyaxis right;
            h(2) = plot(t_grad(filt0),slope_grad(filt0),'m:*','LineWidth',1);
            filt = slope_grad <= Tcrits(1);
            plot(t_grad(filt' & filt0),slope_grad(filt' & filt0),'m-','LineWidth',1.5);


            h(3) = plot(t_diff1(filt1),slope_diff1(filt1),'r:','LineWidth',1);
            filt = slope_diff1 <= Tcrits(2);

            plot(t_diff1(filt & filt1),slope_diff1(filt & filt1),'r-','LineWidth',1.5);
            plot(t_diff2(filt2),slope_diff2(filt2),'rd','LineWidth',1.5,'MarkerFaceColor','r');

            slopemin = floor(nanmin(slope_grad(filt0)));
            slopemax = ceil(nanmax(slope_grad(filt0)));
            if slopemax - slopemin < Tmax - Tmin
                slopemax = slopemin + Tmax - Tmin;           
            end
            set(gca,'Ytick',slopemin:0.5:slopemax)
            ylim([slopemin, slopemax]);
            ylabel('dTa (^oC.10min^-^1)')

            
            ID = num2str(ID1(k));                    


            legend(h(:),{'Ta','dTa_S_D','dTa_N_D_B_C'},'Location','best')

            title(['SD',ID(1:4),' (',xlab(1,1:11),')'])             
            set(gca,'FontSize',12);        
            set(gcf,'Color','w');

            choice = input('Save this fig ? y or n ?')
            if strcmp(choice,'y')
                version = version + 1;
                export_fig([figout,'missing_explanation_v',num2str(version),'_type2.png'],'-dpng','-r200')
            end
            clf; hold on;

        end



    end
    
end



%% plot false alarm CPs
figure('Position',[0 0 800, 500]); hold on;   

t11 = CP_tstart{1,1};
t12 = CP_tend{1,1};
ID1 = CP_sID{1,1};

diff10_far = zeros(size(t21,9));
version = 0;
plt = 1;
for j = 1:9
    j
    t21 = CP_tstart{2,j};
    t22 = CP_tend{2,j};
    ID2 = CP_sID{2,j};

    for k = 1:length(t21)
      
        filt_s = sID == ID2(k);
        t_s = time(filt_s);
        wspd_s = Wspd(filt_s);
        Ta_s = Ta(filt_s);

        if length(unique(t_s)) ~= length(t_s)
            [t_s,ind] = unique(t_s);
            wspd_s = wspd_s(ind,:);
            Ta_s = Ta_s(ind);
        end          
        % interpolate for missing data

        indnan = find(isnan(Ta_s));
        ind = find(isfinite(Ta_s));
        if length(ind) > 0 && length(indnan) > 0
            Ta_s(indnan) = interp1(t_s(ind),Ta_s(ind),t_s(indnan));
        end

        indnan = find(isnan(wspd_s));
        ind = find(isfinite(wspd_s));
        if length(ind) > 0 && length(indnan) > 0
            wspd_s(indnan) = interp1(t_s(ind),wspd_s(ind),t_s(indnan));
        end


        % filter slope data
        filt_sl = slope_sID == ID2(k);
        t_grad = time_slope(filt_sl);
        slope_grad = Tslope(filt_sl,1);


        t_diff1 = min(t_grad):dt*0.1:max(t_grad); % highres diff slope
        slope_diff1 = interp1(time_slope(filt_sl),Tslope(filt_sl,2),t_diff1); 

        t_diff2 = t_grad(1:10:end); % actual res diff slope
        temp = Tslope(filt_sl,2);
        slope_diff2 = temp(1:10:end);   


        filt_case = t_s >= t21(k) - 15*dt & t_s <= t22(k) + 15*dt;       
        Ta_case = Ta_s(filt_case); % case study
        t_case = t_s(filt_case);

        filt0 = t_grad >= t21(k) - 15*dt & t_grad <= t22(k) + 15*dt; 
        filt1 = t_diff1 >= t21(k) - 15*dt & t_diff1 <= t22(k) + 15*dt;
        filt2 = t_diff2 >= t21(k) - 15*dt & t_diff2 <= t22(k) + 15*dt;       


        filt = ((t11 >=t21(k) & t11 <= t22(k)) | (t12 >= t21(k) & t12 <= t22(k)) |...
                (t11 <= t21(k) & t12 >= t22(k))   ) & ID1 == ID2(k);
            
        if sum(filt) > 0
            diff10_far(k,j) = 1;
        elseif plt == 1

            h = [];
            yyaxis left;

            h(1) = plot(t_case,Ta_case,'-o','LineWidth',1.5); 
            Tmin = floor(nanmin(Ta_case));
            Tmax = ceil(nanmax(Ta_case));
            set(gca,'Xtick',t_case(1:5:end),'XtickLabel',xlab(:,13:17),'Ytick',Tmin:0.5:Tmax)
            xlim([min(t_case),max(t_case)])
            xlabel('time'); ylabel('Ta (^oC)'); grid on;
            ylim([Tmin,Tmax])

            yyaxis right;
            h(2) = plot(t_grad(filt0),slope_grad(filt0),'m:*','LineWidth',1);
            filt = slope_grad <= Tcrits(1);
            plot(t_grad(filt' & filt0),slope_grad(filt' & filt0),'m-','LineWidth',1.5);


            h(3) = plot(t_diff1(filt1),slope_diff1(filt1),'r:','LineWidth',1);
            filt = slope_diff1 <= Tcrits(2);

            plot(t_diff1(filt & filt1),slope_diff1(filt & filt1),'r-','LineWidth',1.5);
            plot(t_diff2(filt2),slope_diff2(filt2),'rd','LineWidth',1.5,'MarkerFaceColor','r');

            slopemin = floor(nanmin(slope_grad(filt0)));
            slopemax = ceil(nanmax(slope_grad(filt0)));
            if slopemax - slopemin < Tmax - Tmin
                slopemax = slopemin + Tmax - Tmin;           
            end
            set(gca,'Ytick',slopemin:0.5:slopemax)
            ylim([slopemin, slopemax]);
            ylabel('Ta slope (^oC.10min^-^1)')

            xlab = datestr(t_case(1:5:end));
            ID = num2str(ID1(k));                    


            legend(h(:),{'Ta','dTa_S_D','dTa_N_D_B_C'},'Location','best')

            title(['SD',ID(1:4),' (',xlab(1,1:11),')'])             
            set(gca,'FontSize',12);        
            set(gcf,'Color','w');

            choice = input('Save this fig ? y or n ?')
            if strcmp(choice,'y')
                version = version + 1;
                export_fig([figout,'FAR_explanation_v',num2str(version),'.png'],'-dpng','-r200')
            end
            clf; hold on;

        end
    end
end


diff10_mis_meaning = {'1 = CPs in grad is also in diff; 2 = not in diff due to smoothing; 3 = not in diff due to resolution'};
diff10_far_meaning = {'1 = CPs in diff is also in grad; 0 = not in grad due to slope value'};
% save([matdir,'Saildrone_cold_pools_2018_2023_3Tcrits.mat'],'diff10_mis','diff10_far',...
%     'diff10_mis_meaning','diff10_far_meaning','-append')
    
  
    


        
    


    






































