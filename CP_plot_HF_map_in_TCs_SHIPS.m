
clearvars; close all;
CP_general_parameters

figout = [figdir,'Cold_pool_R1/Maps_in_TCs/HeatFlux/']; 
if ~exist(figout,'dir'); mkdir(figout); end


% load data    
% load([matdir,'Active_gauges_and_CPs_in_SHIPS_motion.mat'])
tag2 = 'motion';
for opt = 2
    load([matdir,'Active_gauges_and_CPs_in_SHIPS_',tag2,'.mat'],'lon_gauge','lat_gauge','dist_gauge',...
        'code_gauge','quadrant_gauge','CP_gauge','shflx_gauge_pre',...
        'lhflx_gauge_pre','lhflx_gauge_post','shflx_gauge_post','Vmax_center',...
        'd2land_center','time_center')

          
    lonR1 = 360-98;
    lonR2 = 360-34;
    latR1 = 11;
    latR2 = 37;

    load([matdir,'NDBC_offshore_buoys_raw_data.mat'],'land_IDs','ocean_IDs');
    land_IDs = land_IDs+13; % code of NDBC starting from 14
    ocean_IDs = ocean_IDs + 13;

    if opt == 2
        ind = setdiff(1:size(lon_gauge,2),land_IDs);
        lon_gauge = lon_gauge(:,ind);
        lat_gauge = lat_gauge(:,ind);
        dist_gauge = dist_gauge(:,ind);
        code_gauge = code_gauge(:,ind);
        quadrant_gauge = quadrant_gauge(:,ind);

        CP_gauge = CP_gauge(:,ind);
        tag = ['_SHIPS_',tag2,'_offshore_gauge'];
    elseif opt == 3
        d2land_min = 250;
        tag = ['_SHIPS_',tag2,'_offshore_TCs',num2str(d2land_min)];
    elseif opt == 4
        d2land_min = 500;
        tag = ['_SHIPS_',tag2,'_offshore_TCs',num2str(d2land_min)];  
    else 
        tag = ['_SHIPS_',tag2];
    end   

    %% plot histogram of quadrant data

    smin = 0;
    smax = 45;
    lmin = 100;
    lmax = 240;



    g_lon = reshape(lon_gauge,[],1);
    g_lat = reshape(lat_gauge,[],1);
    g_dist = reshape(dist_gauge,[],1);
    g_quad = reshape(quadrant_gauge,[],1);
    g_code = reshape(code_gauge,[],1);
    cp_flag = reshape(CP_gauge,[],1);

    g_shflx_pre = reshape(shflx_gauge_pre,[],1);
    g_lhflx_pre = reshape(lhflx_gauge_pre,[],1);
    g_shflx_post = reshape(shflx_gauge_post,[],1);
    g_lhflx_post = reshape(lhflx_gauge_post,[],1);



    [nt,nb] = size(lon_gauge);
    Vmax = reshape(repmat(Vmax_center,[1, nb]),[],1);
    t = reshape(repmat(time_center,[1, nb]),[],1);
    d2land = reshape(repmat(d2land_center,[1, nb]),[],1);
    filt0 = g_lon >= lonR1 & g_lon <= lonR2 & g_lat >= latR1 & g_lat <= latR2 & d2land > 0;


    figure('Position',[0 0 1200 900]);
    nrow = 2; ncol = 2;


    bins = 0.5:1:4.5; % for quadrant count
    Msize = 1000;

    cl1 = {my_color('dull light blue'),my_color('deep sky blue'),my_color('orange')};
    cl2 = {my_color('chartreuse'),my_color('dark green'),my_color('violet red')};

    nbin = length(bins);

    sdat2 = [];
    ldat2 = [];

    

    codefilt = isfinite(g_code);  

    
    sub = 0;
    for i = 1:2 % CP vs pre CP
        CP_filt = cp_flag == 1;
        if i == 1           
            g_shflx = g_shflx_post;
            g_lhflx = g_lhflx_post;
        else
            g_shflx = g_shflx_pre;
            g_lhflx = g_lhflx_pre;
        end  
        
        sdat = []; 
        sdat_ci = [];
        
        ldat = []; 
        ldat_ci = [];        
        
        for opt = 1:2 % inside vs. outside

            if opt == 1
                dist_filt = abs(g_dist) <= 500;
            else
                dist_filt = abs(g_dist) > 500 & abs(g_dist) <= 1000 ;
            end    


            
            filt = filt0 & codefilt & dist_filt & CP_filt;
            quad = g_quad(filt);
            sflux = g_shflx(filt);
            lflux = g_lhflx(filt);            

            [Smean, S_CI] = coldpool3_function_histogram_bootstrap(sflux,quad,bins,Msize);
            [Lmean, L_CI] = coldpool3_function_histogram_bootstrap(lflux,quad,bins,Msize);
            
            if opt == 1
                sdat = cat(1,sdat,Smean);
                sdat_ci = cat(1,sdat_ci,S_CI);
                ldat = cat(1, ldat,Lmean);
                ldat_ci = cat(1,ldat_ci,L_CI); 
            else
                sdat = cat(1,sdat,Smean(nbin));
                sdat_ci = cat(1,sdat_ci,S_CI(nbin,:));
                ldat = cat(1,ldat,Lmean(nbin));
                ldat_ci = cat(1,ldat_ci,L_CI(nbin,:));
            end
            
        end
        sdat2 = cat(2,sdat2,sdat);
        ldat2 = cat(2,ldat2,ldat);
        
        sdat_err = sdat_ci - sdat;   
        ldat_err = ldat_ci - ldat;  

        sub = sub+1; subplot(nrow,ncol,sub); hold on; 
        bar(1:4,sdat(1:4),'FaceColor',cl1{1})
        bar(5,sdat(5),'FaceColor',cl1{2})
        bar(6,sdat(6),'FaceColor',cl1{3})
        errorbar(1:6,sdat,sdat_err(:,1),sdat_err(:,2),'k.'); 
        
        ylabel('W.m^-^2'); xlabel('TC quadrant');
        if strcmp(tag2,'VWS')
            set(gca,'Xtick',1:6,'XtickLabel',{'DR','UR','UL','DL','inside','outside'},'FontSize',12)
        else
            set(gca,'Xtick',1:6,'XtickLabel',{'FR','RR','RL','FL','inside','outside'},'FontSize',12)
        end
        text(0.5,smax,fignum{sub},'FontSize',16,'FontWeight','bold')
        ylim([smin, smax]);
        if i == 1
            title('sensible flux (TC-CP)')
        else
            title('sensible flux (TC-preCP)')
        end
        set(gca,'Ygrid','on')
        
        sub = sub+1; subplot(nrow,ncol,sub); hold on; 
        bar(1:4,ldat(1:4),'FaceColor',cl2{1})
        bar(5,ldat(5),'FaceColor',cl2{2})
        bar(6,ldat(6),'FaceColor',cl2{3})
        errorbar(1:6,ldat,ldat_err(:,1),ldat_err(:,2),'k.');         
        
        ylabel('W.m^-^2'); xlabel('TC quadrant');
        if strcmp(tag2,'VWS')
            set(gca,'Xtick',1:6,'XtickLabel',{'DR','UR','UL','DL','inside','outside'},'FontSize',12)
        else
            set(gca,'Xtick',1:6,'XtickLabel',{'FR','RR','RL','FL','inside','outside'},'FontSize',12)
        end
        set(gca,'Ygrid','on')
        ylim([lmin, lmax]);
        if i == 1
            title('latent flux (TC-CP)')
        else
            title('latent flux (TC-preCP)')
        end    
        text(0.5,lmax,fignum{sub},'FontSize',16,'FontWeight','bold')

    end
    
%     axes( 'Position', [0, 0.96, 1, 0.04] ) ;
%     set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
%     text( 0.5, 0, ['Mean heat flux exchange in TC quadrants (',datname{d},')'], 'FontSize', 20, 'FontWeight', 'Bold', ...
%       'HorizontalAlignment', 'Center', 'VerticalAlignment', 'bottom','Color','r' ) ;  
    set(gcf,'Color','w');   pause
    export_fig([figout,'HF_in_TC_quadrants',tag,'.png'],'-dpng','-r200')
    clf;    
    
end       





 








































