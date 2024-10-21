
clearvars; close all;
CP_general_parameters

figout = [figdir,'Cold_pool_R1/Maps_in_TCs/HeatFlux/']; 
if ~exist(figout,'dir'); mkdir(figout); end


% load data    
% load([matdir,'Active_gauges_and_CPs_in_SHIPS_motion.mat'])
for opt = 1%2:4
    load([matdir,'Active_gauges_and_CPs_in_TC_coordinates.mat'],'lon_gauge','lat_gauge','dist_gauge',...
        'code_gauge','quadrant_gauge','CP_gauge','shflx_gauge_pre',...
        'lhflx_gauge_pre','lhflx_gauge_post','shflx_gauge_post','wspd_gauge_pre','Vmax_center',...
        'wspd_gauge_post','sst_gauge_pre','sst_gauge_post','Ta_gauge_pre','Ta_gauge_post','d2land_center','time_center')

          
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
        tag = '_IBTRACS_offshore_gauge';
    elseif opt == 3
        d2land_min = 250;
        tag = ['_IBTRACS_offshore_TCs',num2str(d2land_min)];
    elseif opt == 4
        d2land_min = 500;
        tag = ['_IBTRACS_offshore_TCs',num2str(d2land_min)];  
    else 
        tag = '_IBTRACS';
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
        set(gca,'Xtick',1:6,'XtickLabel',{'FR','RR','RL','FL','inside','outside'},'FontSize',12)
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
        set(gca,'Xtick',1:6,'XtickLabel',{'FR','RR','RL','FL','inside','outside'},'FontSize',12)
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



%% impacts of Vmax and SST
load([matdir,'Active_gauges_and_CPs_in_TC_coordinates.mat'],'lon_gauge','lat_gauge','dist_gauge',...
    'code_gauge','quadrant_gauge','CP_gauge','shflx_gauge_pre',...
    'lhflx_gauge_pre','lhflx_gauge_post','shflx_gauge_post','wspd_gauge_pre','Vmax_center',...
    'wspd_gauge_post','sst_gauge_pre','sst_gauge_post','Ta_gauge_pre','Ta_gauge_post','d2land_center','time_center')


lonR1 = 360-98;
lonR2 = 360-34;
latR1 = 11;
latR2 = 37;


datname = {'PIRATA','NDBC','Saildrone','All'};
minval = [28, 5.5, 26.3];
maxval = [29.5,10,28.8];

distname = {'outside 500 km','within 500 km'};

g_lon = reshape(lon_gauge,[],1);
g_lat = reshape(lat_gauge,[],1);
g_dist = reshape(dist_gauge,[],1);
g_quad = reshape(quadrant_gauge,[],1);
g_code = reshape(code_gauge,[],1);
cp_flag = reshape(CP_gauge,[],1);


g_sst_pre = reshape(sst_gauge_pre,[],1);
g_wspd_pre = reshape(wspd_gauge_pre,[],1);
g_Ta_pre = reshape(Ta_gauge_pre,[],1);

g_sst_post = reshape(sst_gauge_post,[],1);
g_wspd_post = reshape(wspd_gauge_post,[],1);
g_Ta_post = reshape(Ta_gauge_post,[],1);

g_sst = reshape(sst_gauge,[],1);
[nt,nb] = size(lon_gauge);
Vmax = reshape(repmat(Vmax_center,[1, nb]),[],1);
d2land = reshape(repmat(d2land_center,[1, nb]),[],1);
filt0 = g_lon >= lonR1 & g_lon <= lonR2 & g_lat >= latR1 & g_lat <= latR2 & d2land > 0;


figure('Position',[0 0 1200 900]);
nrow = 2; ncol = 3;
subnum = [1 4 7 10 2 5 8 11 3 6 9 12];

bins = 0.5:1:4.5; % for quadrant count
Msize = 1000;

cl1 = {my_color('dull light blue'),my_color('deep sky blue'),my_color('orange')};
cl2 = {my_color('chartreuse'),my_color('dark green'),my_color('violet red')};

nbin = length(bins);

sst2 = [];
wspd2 = [];
Ta2 = [];
for d = 4%1:length(datname)
    
    if d < 4
        codefilt = g_code == d;        
    else
        codefilt = isfinite(g_code);  
    end 
    
    sub = 0;
    for i = 1:2 % CP vs pre CP
        CP_filt = cp_flag == 1;
        if i == 1           
            g_sst = g_sst_post;
            g_wspd = g_wspd_post;
            g_Ta = g_Ta_post;
        else
            g_sst = g_sst_pre;
            g_wspd = g_wspd_pre;
            g_Ta = g_Ta_pre;
        end  
        
        sst = []; 
        sst_ci = [];
        
        wspd = []; 
        wspd_ci = []; 
        
        Ta = [];
        Ta_ci = [];
        
        for opt = 1:2 % inside vs. outside

            if opt == 1
                dist_filt = abs(g_dist) <= 500;
            else
                dist_filt = abs(g_dist) > 500 & abs(g_dist) <= 1000 ;
            end    


            
            filt = filt0 & codefilt & dist_filt & CP_filt;
            quad = g_quad(filt);
            sst_opt = g_sst(filt);
            wspd_opt = g_wspd(filt);  
            Ta_opt = g_Ta(filt);

            [sst_mean, s_CI] = coldpool3_function_histogram_bootstrap(sst_opt,quad,bins,Msize);
            [wspd_mean, w_CI] = coldpool3_function_histogram_bootstrap(wspd_opt,quad,bins,Msize);
            [Ta_mean, T_CI] = coldpool3_function_histogram_bootstrap(Ta_opt,quad,bins,Msize);
            if opt == 1
                sst = cat(1,sst,sst_mean);
                sst_ci = cat(1,sst_ci,s_CI);
                
                wspd = cat(1, wspd,wspd_mean);
                wspd_ci = cat(1,wspd_ci,w_CI); 
                
                Ta = cat(1,Ta,Ta_mean);
                Ta_ci = cat(1,Ta_ci,T_CI);                
                
            else
                sst = cat(1,sst,sst_mean(nbin));
                sst_ci = cat(1,sst_ci,s_CI(nbin,:));
                
                Ta = cat(1,Ta,Ta_mean(nbin));
                Ta_ci = cat(1,Ta_ci,T_CI(nbin,:));                
                
                wspd = cat(1,wspd,wspd_mean(nbin));
                wspd_ci = cat(1,wspd_ci,w_CI(nbin,:));
            end
            
        end
        sst2 = cat(2,sst2,sst);
        wspd2 = cat(2,wspd2,wspd);
        Ta2 = cat(2,Ta2,Ta);
        
        
        sst_err = sst_ci - sst;   
        wspd_err = wspd_ci - wspd;  
        Ta_err = Ta_ci - Ta;

        sub = sub+1; subplot(nrow,ncol,sub); hold on; 
        bar(1:4,sst(1:4),'FaceColor',cl1{1})
        bar(5,sst(5),'FaceColor',cl1{2})
        bar(6,sst(6),'FaceColor',cl1{3})
        errorbar(1:6,sst,sst_err(:,1),sst_err(:,2),'k.'); 
        
        ylabel('^oC'); xlabel('TC quadrant');
        set(gca,'Xtick',1:6,'XtickLabel',{'FR','RR','RL','FL','inside','outside'},'FontSize',12)
        set(gca,'Ygrid','on')
        text(0.5,maxval(1),fignum{sub},'FontSize',16,'FontWeight','bold')
        ylim([minval(1), maxval(1)]);
        if i == 1
            title('SST (TC-CP)')
        else
            title('SST (TC-preCP)')
        end

        sub = sub+1; subplot(nrow,ncol,sub); hold on; 
        bar(1:4,wspd(1:4),'FaceColor',cl2{1})
        bar(5,wspd(5),'FaceColor',cl2{2})
        bar(6,wspd(6),'FaceColor',cl2{3})
        errorbar(1:6,wspd,wspd_err(:,1),wspd_err(:,2),'k.');         
        
        ylabel('m.s^-^1'); xlabel('TC quadrant');
        set(gca,'Xtick',1:6,'XtickLabel',{'FR','RR','RL','FL','inside','outside'},'FontSize',12)
        set(gca,'Ygrid','on')
        ylim([minval(2), maxval(2)]);

        if i == 1
            title('WSPD (TC-CP)')
        else
            title('WSPD (TC-preCP)')
        end    
        text(0.5,maxval(2),fignum{sub},'FontSize',16,'FontWeight','bold')
        
        
        sub = sub+1; subplot(nrow,ncol,sub); hold on; 
        bar(1:4,Ta(1:4),'FaceColor',cl2{1})
        bar(5,Ta(5),'FaceColor',cl2{2})
        bar(6,Ta(6),'FaceColor',cl2{3})
        errorbar(1:6,Ta,Ta_err(:,1),Ta_err(:,2),'k.');         
        
        ylabel('^oC'); xlabel('TC quadrant');
        set(gca,'Xtick',1:6,'XtickLabel',{'FR','RR','RL','FL','inside','outside'},'FontSize',12)
        set(gca,'Ygrid','on')
        
        ylim([minval(3),maxval(3)]);
        if i == 1
            title('Ta (TC-CP)')
        else
            title('Ta (TC-preCP)')
        end    
        text(0.5,maxval(3),fignum{sub},'FontSize',16,'FontWeight','bold')        

    end
    
%     axes( 'Position', [0, 0.96, 1, 0.04] ) ;
%     set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
%     text( 0.5, 0, ['Mean heat flux exchange in TC quadrants (',datname{d},')'], 'FontSize', 20, 'FontWeight', 'Bold', ...
%       'HorizontalAlignment', 'Center', 'VerticalAlignment', 'bottom','Color','r' ) ;  
    set(gcf,'Color','w');   pause
    export_fig([figout,'Ta_Wspd_SST_quadrants_barplot_IBTRACS.png'],'-dpng','-r200')
    clf;    
    
end   


 








































