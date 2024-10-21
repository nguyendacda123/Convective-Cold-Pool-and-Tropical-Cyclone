
clearvars; close all;
saildrone_eparams

figout = [figdir,'Cold_pool4/Role_of_SST/']; 
if ~exist(figout,'dir'); mkdir(figout); end


%% load data
load([matdir,'Active_gauges_and_CPs_in_TC_coordinates.mat'],'lon_center','lat_center',...
    'time_center','Vmax_center','Vtrans_center','TCI_center','lon_gauge','lat_gauge','dist_gauge',...
    'code_gauge','angle_gauge','quadrant_gauge','TC_code','phase_gauge','CP_gauge','CP_presst',...
    'mdist_500','mdist_200','Dist_min','tDist_min','sst_gauge','lhflx_gauge','shflx_gauge')             
           
                
lonR1 = 360-98;
lonR2 = 360-34;
latR1 = 11;
latR2 = 37;
filt0 = lon_center >= lonR1 & lon_center <= lonR2 & lat_center >= latR1 & lat_center <= latR2;
[nt,nb] = size(quadrant_gauge);

Vmax = repmat(Vmax_center,[1, nb]);
Vtrans = repmat(Vtrans_center,[1, nb]);
% compute cold pool vs. distance
tag = {'TD-TS','Cat1-5','slow6.2','fast6.2','All'};



%% plot SST distribution associated with CPs, with TCs and the ratio




datname = {'PIRATA','NDBC','Saildrone','All'};
s1 = [0, 0, 0, 0];
s2 = [20, 30, 140, 45];
l1 = [70, 100, 130, 100];
l2 = [160, 240, 450, 240];

distname = {'outside 500 km','within 500 km'};

g_lon = reshape(lon_gauge,[],1);
g_lat = reshape(lat_gauge,[],1);
g_dist = reshape(dist_gauge,[],1);
g_quad = reshape(quadrant_gauge,[],1);
g_code = reshape(code_gauge,[],1);
cp_flag = reshape(CP_gauge,[],1);
mind500 = reshape(mdist_500,[],1);
g_sst = reshape(sst_gauge,[],1);
g_lhflx = reshape(lhflx_gauge,[],1);

filt0 = g_lon >= lonR1 & g_lon <= lonR2 & g_lat >= latR1 & g_lat <= latR2 & mind500 == 1;


figure('Position',[0 0 1200 400]);
nrow = 1; ncol = 2;
subnum = [1 4 7 10 2 5 8 11 3 6 9 12];

bins = 0.5:1:4.5; % for quadrant count
Msize = 1000;

cl1 = {my_color('dull light blue'),my_color('deep sky blue'),my_color('orange')};
cl2 = {my_color('chartreuse'),my_color('dark green'),my_color('violet red')};

nbin = length(bins);

sdat2 = [];
ldat2 = [];
for d = 4%length(datname)
    
    if d < 4
        codefilt = g_code == d;        
    else
        codefilt = isfinite(g_code);  
    end 
    
    sub = 0;
    for i = 1:2 % CP vs non CP
        if i == 1
            CP_filt = cp_flag == 1;
        else
            CP_filt = cp_flag == 0;
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
            sst = g_sst(filt);
            lflux = g_lhflx(filt);            

            [Smean, S_CI] = coldpool3_function_histogram_bootstrap(sst,quad,bins,Msize);
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
        
        ylabel('^oC'); %xlabel('TC quadrant');
        set(gca,'Xtick',1:6,'XtickLabel',{'FR','RR','RL','FL','inside','outside'},'FontSize',12)
        text(0.8,29.4,fignum{sub},'FontSize',16,'FontWeight','bold')
        ylim([28, 29.5]);
        if i == 1
            title('SST (TC-CPs)')
        else
            title('SST (TC-nonCPs)')
        end


    end
    
%     axes( 'Position', [0, 0.96, 1, 0.04] ) ;
%     set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
%     text( 0.5, 0, ['Mean heat flux exchange in TC quadrants (',datname{d},')'], 'FontSize', 20, 'FontWeight', 'Bold', ...
%       'HorizontalAlignment', 'Center', 'VerticalAlignment', 'bottom','Color','r' ) ; 
    
    set(gcf,'Color','w');   pause
    export_fig([figout,datname{d},'_SST_in_TC_quadrants_barplot.png'],'-dpng','-r200')
    clf;    
    
end














































