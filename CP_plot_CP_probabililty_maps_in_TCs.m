
clearvars; close all;
CP_general_parameters

figout = [figdir,'Cold_pool_R1/Maps_in_TCs/Frequency/']; 
if ~exist(figout,'dir'); mkdir(figout); end


% load data              
% load([matdir,'Active_gauges_and_CPs_in_TC_coordinates.mat'],'lon_center','lat_center',...
%     'time_center','Vmax_center','Vtrans_center','TCI_center','lon_gauge','lat_gauge','dist_gauge',...
%     'code_gauge','angle_gauge','quadrant_gauge','TC_code','phase_gauge','CP_gauge',...
%     'mdist_500','mdist_200','Dist_min','tDist_min','R34_center','R50_center','R64_center','Rmax_center')
figure('Position',[0 0 1000 1100]); 
nrow = 3; ncol = 2; 


load([matdir,'NDBC_offshore_buoys_raw_data.mat'],'land_IDs','ocean_IDs');
land_IDs = land_IDs+13; % code of NDBC starting from 14
ocean_IDs = ocean_IDs + 13;

for opt = 1%3:4
    sub = 0;

    for d = 1:3
        if d == 1
            load([matdir,'Active_gauges_and_CPs_in_TC_coordinates.mat'],'lon_gauge','lat_gauge','dist_gauge','time_center','d2land_center',...
                'quadrant_gauge','CP_gauge','R34_center','R50_center','R64_center','Rmax_center','Vmax_center')
        datname = 'IBTrACS motion';
        elseif d == 2
            load([matdir,'Active_gauges_and_CPs_in_SHIPS_motion.mat'],'lon_gauge','lat_gauge','dist_gauge','d2land_center',...
                'quadrant_gauge','CP_gauge','R34_center','R50_center','R64_center','Rmax_center','Vmax_center','time_center')        
        datname = 'SHIPS motion';
        else
            load([matdir,'Active_gauges_and_CPs_in_SHIPS_VWS.mat'],'lon_gauge','lat_gauge','dist_gauge','d2land_center',...
                'quadrant_gauge','CP_gauge','R34_center','R50_center','R64_center','Rmax_center','Vmax_center','time_center')  
        datname = 'SHIPS shear';
        end
        

        if opt == 2
            ind = setdiff(1:size(lon_gauge,2),land_IDs);
            lon_gauge = lon_gauge(:,ind);
            lat_gauge = lat_gauge(:,ind);
            dist_gauge = dist_gauge(:,ind);
            quadrant_gauge = quadrant_gauge(:,ind);
            CP_gauge = CP_gauge(:,ind);
            tag = '_offshore_gauge';
            d2land_min = 0;
        elseif opt == 3
            d2land_min = 250;
            tag = ['_IBTRACS_offshore_TCs',num2str(d2land_min)];
        elseif opt == 4
            d2land_min = 500;
            tag = ['_IBTRACS_offshore_TCs',num2str(d2land_min)];  
        else 
            tag = '';
            d2land_min = 0;
        end
            

        loncp = reshape(lon_gauge,[],1);
        latcp = reshape(lat_gauge,[],1);
        Vmax = reshape(repmat(Vmax_center,1,size(lon_gauge,2)),[],1);
        t = reshape(repmat(time_center,1,size(lon_gauge,2)),[],1);
        d2land = reshape(repmat(d2land_center,1,size(lon_gauge,2)),[],1);
        R = abs(reshape(dist_gauge,[],1));

        cp_flag = reshape(CP_gauge,[],1);
        quad = reshape(quadrant_gauge,[],1);
        RMW = reshape(repmat(Rmax_center,1,size(lon_gauge,2)),[],1);
        rRMW = R./RMW;

        filt0 = loncp >= LONR1 & loncp <= LONR2 & latcp >= LATR1 & latcp <= LATR2 & ...
            d2land > d2land_min;



        % compute probability
        Rbin = 0:50:500;
        nbin = length(Rbin)-1;
        freq_g = nan(nbin,4);
        freq_cp = nan(nbin,4);
        prob = nan(361,Rbin(end)+1); 
        for q = 1:4
            for k = 1:nbin
                filt1 = filt0 & quad == q & R >= Rbin(k) & R < Rbin(k+1);
                filt2 = filt1 & cp_flag == 1;
                freq_g(k,q) = sum(filt1);
                freq_cp(k,q) = sum(filt2);

                prob((q-1)*90+1:q*90+1,Rbin(k)+1:Rbin(k+1)+1) = freq_cp(k,q)/freq_g(k,q)*100;
            end

        end

        % compute area and mean probability
        S = nan(360,Rbin(end));
        for k = 2:Rbin(end)
            S(:,k) = (pi*k^2 - pi*(k-1)^2)/360;
        end
        meanprob = nansum(nansum(prob(1:360,1:500).*S))/nansum(S(:))
        minprob = min(prob(:))
        maxprob = max(prob(:))
        % compute probability in relative coordinate
        Rbin2 = 0:0.5:5;
        nbin = length(Rbin2)-1;
        freq_g = nan(nbin,4);
        freq_cp = nan(nbin,4);
        prob2 = nan(361,nbin+1); 
        for q = 1:4
            for k = 1:nbin
                filt1 = filt0 & quad == q & rRMW >= Rbin2(k) & rRMW <= Rbin2(k+1);
                filt2 = filt1 & cp_flag == 1;
                freq_g(k,q) = sum(filt1);
                freq_cp(k,q) = sum(filt2);
                prob2((q-1)*90+1:q*90+1,Rbin(k)+1:Rbin(k+1)+1) = freq_cp(k,q)/freq_g(k,q)*100;
            end

        end




        sub = sub+1; subplot(nrow,ncol,sub); 
        ticklabel = {'0','100', '200', '300', '400', '500 km'};
        polarPcolor(0:500,0:360,prob,'Ncircles',6,'Nspokes',1,'RtickLabel',ticklabel,'FontSize',12)
%         quiver(0,1,0,0.8,'LineWidth',2,'Color','k','MaxHeadSize',1)
%         xlim([-1.1,1.1]);ylim([-1.1,1.1])
        if d < 3
            text([0.8,0.8,-1,-1],[0.7,-0.7,-0.7,0.7],{'FR','RR','RL','FL'},...
                'FontSize',14,'Color','k')
        else
            text([0.8,0.8,-1,-1],[0.7,-0.7,-0.7,0.7],{'DR','UR','UL','DL'},...
                'FontSize',14,'Color','k')
        end
        text(-0.8,0.8,datname,'FontSize',16,'Color','k','FontWeight','bold','Rotation',90)
%         text(-0.8,-0.8,'%','FontSize',14,'Color','k','FontWeight','bold')
        text(-0.95,0.95,fignum{sub},'FontSize',16,'Color','k','FontWeight','bold')
        if d == 1
            text(0,0.8,'R (km)','FontSize',16,'Color','k','FontWeight','bold')
        end
        set(gca,'FontSize',14);
        caxis([0,15])


        sub = sub+1; subplot(nrow,ncol,sub); 
        ticklabel = {'0','1', '2', '3', '4', '5'};
        polarPcolor(0:0.01:5,0:360,prob2,'Ncircles',6,'Nspokes',1,'RtickLabel',ticklabel,'FontSize',12)
%         quiver(0,1,0,0.8,'LineWidth',2,'Color','k','MaxHeadSize',1)
%         xlim([-1.1,1.1]);ylim([-1.1,1.1])
        
        text(-0.95,0.95,fignum{sub},'FontSize',16,'Color','k','FontWeight','bold')
        set(gca,'FontSize',14); set(gcf,'Color','w'); 
        if d == 1
            text(-0.8,-0.8,'%','FontSize',14,'Color','k','FontWeight','bold')
            text(0,0.8,'R/RMW','FontSize',16,'Color','k','FontWeight','bold')
        end
        
        caxis([0,15])
        
%         if d < 3
%             text([0.8,0.8,-0.8,-0.8],[0.7,-0.7,-0.7,0.7],{'FR','RR','RL','FL'},...
%                 'FontSize',14,'Color','k')
%         else
%             text([0.8,0.8,-0.8,-0.8],[0.7,-0.7,-0.7,0.7],{'DR','UR','UL','DL'},...
%                 'FontSize',14,'Color','k')
%         end        
    end
    colormap('bluewhitered'); 
    pause
    export_fig([figout,'CP_probability_quadrant_distance',tag,'.png'],'-dpng','-r200')
    pause(0.1); clf
end


%% using R30/50/64

clearvars; close all;
CP_general_parameters

figout = [figdir,'Cold_pool4/Maps_in_TCs/']; 
if ~exist(figout,'dir'); mkdir(figout); end


% load data              
% load([matdir,'Active_gauges_and_CPs_in_TC_coordinates.mat'],'lon_center','lat_center',...
%     'time_center','Vmax_center','Vtrans_center','TCI_center','lon_gauge','lat_gauge','dist_gauge',...
%     'code_gauge','angle_gauge','quadrant_gauge','TC_code','phase_gauge','CP_gauge',...
%     'mdist_500','mdist_200','Dist_min','tDist_min','R34_center','R50_center','R64_center','Rmax_center')
load([matdir,'Active_gauges_and_CPs_in_TC_coordinates.mat'],'lon_gauge','lat_gauge','dist_gauge',...
    'quadrant_gauge','CP_gauge','mdist_500','R34_center','R50_center','R64_center','Rmax_center')
          
LONR1 = LONR1;
LONR2 = LONR2;
LATR1 = LATR1;
LATR2 = LATR2;


loncp = reshape(lon_gauge,[],1);
latcp = reshape(lat_gauge,[],1);
R = abs(reshape(dist_gauge,[],1));
mind500 = reshape(mdist_500,[],1);
cp_flag = reshape(CP_gauge,[],1);
quad = reshape(quadrant_gauge,[],1);

filt0 = loncp >= LONR1 & loncp <= LONR2 & latcp >= LATR1 & latcp <= LATR2 & mind500 == 1;

Rmax_center = repmat(Rmax_center,1,4);
RMW = cat(3,R34_center,R50_center,R64_center,Rmax_center);
Rname = {'R34','R50','R64','Rmax'};
nR = length(Rname)
figure('Position',[0 0 1200 1000]); nrow = 2; ncol = 2;
for i = 1:nR
    dat = squeeze(RMW(:,:,i));
    Rbin = 0:50:500;
    Rbin2 = 0:0.5:5;
    nbin = length(Rbin2)-1;
    freq_g = nan(nbin,4);
    freq_cp = nan(nbin,4);
    prob2 = nan(361,nbin+1); 
    for q = 1:4
        rRMW = R./reshape(repmat(dat(:,q),1,size(lon_gauge,2)),[],1);
        for k = 1:nbin
            filt1 = filt0 & quad == q & rRMW >= Rbin2(k) & rRMW <= Rbin2(k+1);
            filt2 = filt1 & cp_flag == 1;
            freq_g(k,q) = sum(filt1);
            freq_cp(k,q) = sum(filt2);
            prob2((q-1)*90+1:q*90+1,Rbin(k)+1:Rbin(k+1)+1) = freq_cp(k,q)/freq_g(k,q)*100;
        end

    end    
    
    h2 = subplot(nrow,ncol,i);
    ticklabel = {'0','1', '2', '3', '4', '5'};
    polarPcolor(0:0.01:5,0:360,prob2,'Ncircles',6,'Nspokes',1,'RtickLabel',ticklabel,'FontSize',12)

    text(-0.8,-0.8,'%','FontSize',14,'Color','k','FontWeight','bold')
    text(-0.8,0.8,fignum{i},'FontSize',16,'Color','k','FontWeight','bold')
    if i == 1
        text([0.8,0.8,-0.8,-0.8],[0.7,-0.7,-0.7,0.7],{'FR','RR','RL','FL'},...
        'FontSize',14,'Color','k')
    end
    set(gca,'FontSize',14);
    set(gcf,'Color','w'); 
    text(0,0.8,['R/',Rname{i}],'FontSize',16,'Color','k','FontWeight','bold')
    
    
end
colormap('jet'); set(gcf,'Color','w'); pause
export_fig([figout,'CP_probability_quadrant_distance_relative_multiR.png'],'-dpng','-r200')































