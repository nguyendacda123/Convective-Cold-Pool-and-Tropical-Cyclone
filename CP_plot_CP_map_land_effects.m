
clearvars; close all;
CP_general_parameters

     
load([matdir,'NDBC_offshore_buoys_raw_data.mat'],'land_IDs','ocean_IDs');
coastal_IDs = land_IDs+13; % code of NDBC starting from 14
offshore_IDs = ocean_IDs + 13;

figout = [figdir,'Cold_pool_R1/summary/']; 
if ~exist(figout,'dir'); mkdir(figout); end

datname = {' (>50 km)',' (>200 km)'};

figure('Position',[0 0 1100 1000]);
nrow = 2; ncol = 2;
sub_order = [1, 3, 2, 4];
sub = 0;
for d = 1:length(datname)
    
    load([matdir,'Active_gauges_and_CPs_in_TC_coordinates.mat'],'lon_center','lat_center',...
        'time_center','Vmax_center','Vtrans_center','TCI_center','lon_gauge','lat_gauge','dist_gauge',...
        'code_gauge','angle_gauge','quadrant_gauge','sst_gauge','TC_code','phase_gauge','CP_gauge','CP_presst',...
        'd2land_center') 
    
    if d == 1
        ind = cat(2,coastal_IDs, offshore_IDs);
    else

        ind = offshore_IDs;
    end
    lon_gauge = lon_gauge(:,ind);
    lat_gauge = lat_gauge(:,ind);
    dist_gauge = dist_gauge(:,ind);
    code_gauge = code_gauge(:,ind);
    angle_gauge = angle_gauge(:,ind);
    quadrant_gauge = quadrant_gauge(:,ind);
    phase_gauge = phase_gauge(:,ind);
    CP_gauge = CP_gauge(:,ind);
        
    
    


    loncp = reshape(lon_gauge,[],1);
    latcp = reshape(lat_gauge,[],1);
    R = abs(reshape(dist_gauge,[],1));
    phi = reshape(angle_gauge,[],1);
    type = reshape(code_gauge,[],1);

    cp_flag = reshape(CP_gauge,[],1);
    quad = reshape(quadrant_gauge,[],1);
    d2land = reshape(repmat(d2land_center,1,size(lon_gauge,2)),[],1);

    filt0 = loncp >= LONR1 & loncp <= LONR2 & latcp >= LATR1 & latcp <= LATR2 & cp_flag == 1 & d2land >0;

    R = R(filt0);
    phi = phi(filt0);
    type = type(filt0);
    quad = quad(filt0);

    phic = linspace(0,2*pi,721);
    ind1 = find(phic == 0);
    ind2 = find(phic == pi);
    ind3 = find(phic == pi/2);
    ind4 = find(phic == 3*pi/2);

    sub = sub+1; subplot(nrow,ncol,sub_order(sub)); hold on;            
    Rc = [100:200:500];
    filt = R <= 500 & isfinite(type);
    X = R(filt).*sin(phi(filt)/180*pi);
    Y = R(filt).*cos(phi(filt)/180*pi); 
    scatter(X,Y,15,'filled','MarkerFaceColor',my_color('dark green')); 

    for k = 1:length(Rc)
        Xc = Rc(k).*sin(phic);
        Yc = Rc(k).*cos(phic);
        plot(Xc,Yc,'k--','LineWidth',2);
        text(Xc(181),Yc(181)+20,[num2str(Rc(k))],'FontSize',14,'Color','r','FontWeight','bold');

    end        
    plot(Xc([ind1,ind2]),Yc([ind1,ind2]),'k--','LineWidth',2);
    plot(Xc([ind3,ind4]),Yc([ind3,ind4]),'k--','LineWidth',2);        

    xlabel('Distance from TC center (km)'); ylabel('Distance from TC center (km)');

    if d == 1
        text([400,400,-400,-400],[400,-400,-400,400],{'FR','RR','RL','FL'},...
            'Color','k','FontSize',14,'FontWeight','bold')       
    end
    text(-400,400,fignum{sub},'FontSize',16,'FontWeight','bold')
    set(gca,'FontSize',14);
           
    title(['CP map',datname{d}])



    %% plot histogram of quadrant data

    cp_dist = reshape(dist_gauge,[],1);
    g_dist = reshape(dist_gauge,[],1);
    cp_quad = reshape(quadrant_gauge,[],1);
    cp_code = reshape(code_gauge,[],1);
    g_quad = reshape(quadrant_gauge,[],1);
    g_code = reshape(code_gauge,[],1);
    cp_flag = reshape(CP_gauge,[],1);
    d2land = reshape(repmat(d2land_center,1,size(lon_gauge,2)),[],1);

    filt0 = loncp >= LONR1 & loncp <= LONR2 & latcp >= LATR1 & latcp <= LATR2 & d2land > 0;

    cp_dist = cp_dist(filt0 & cp_flag == 1);
    g_dist = g_dist(filt0);
    cp_quad = cp_quad(filt0 & cp_flag == 1);
    cp_code = cp_code(filt0 & cp_flag == 1);
    g_quad = g_quad(filt0);
    g_code = g_code(filt0);



    bins = 0.5:1:4.5; % for quadrant count
    Msize = 1000;

    cl = {my_color('dark steel blue'),my_color('dull light blue'),my_color('deep sky blue')};
    tab = [];


    cfilt1 = isfinite(cp_code);
    cfilt2 = isfinite(g_code);  



    sub = sub+1; subplot(nrow,ncol,sub_order(sub)); hold on;    

    dist_filt1 = abs(cp_dist) <= 500;
    dist_filt2 = abs(g_dist) <= 500;


    filt1 = dist_filt1 & cfilt1;
    filt2 = dist_filt2 & cfilt2;
    dat1 = cp_quad(filt1);
    dat2 = g_quad(filt2);

    [freq,CI] = coldpool2_function_histogram_bootstrap_v2(dat1,dat2,bins,Msize);


    err = CI - freq;                
    bar(1:4,freq,'FaceColor',cl{3})
    errorbar(1:4,freq,err(:,1),err(:,2),'k.');  

    
    ylabel('min.day^-^1');
    ylim([50,140])

    tab = cat(2,tab,freq);   
    
    title(['CP/sampling',datname{d}])

    xlabel('TC quadrant')
    
    set(gca,'Xtick',1:4,'XtickLabel',{'FR','RR','RL','FL'})
    h= quiver(0,500,0,80,'LineWidth',2,'Color','k','MaxHeadSize',1)
    xlim([-600, 600]); ylim([-600, 600])    
    set(gca,'FontSize',14,'Ygrid','on');
    text(0.5,nanmax(freq),fignum{sub_order(sub)},'FontSize',16,'FontWeight','bold')


end


set(gcf,'Color','w'); pause
export_fig([figout,'Land_effects_summary.png'],'-dpng','-r200')    
pause(0.1);
clf;



%% plot sensitivity
figure('Position',[0 0 1200 800]);
nrow = 3; ncol = 2;
sub = 0;
order = [1,3,5,2,4,6];
for d = 1:2
    % load data              
    load([matdir,'Active_gauges_and_CPs_in_TC_coordinates.mat'],'lon_center','lat_center',...
        'time_center','Vmax_center','Vtrans_center','TCI_center','lon_gauge','lat_gauge','dist_gauge',...
        'code_gauge','angle_gauge','quadrant_gauge','sst_gauge','TC_code','phase_gauge','CP_gauge','CP_presst',...
        'd2land_center')

    load([matdir,'NDBC_offshore_buoys_raw_data.mat'],'land_IDs','ocean_IDs');
    land_IDs = land_IDs+13; % code of NDBC starting from 14
    ocean_IDs = ocean_IDs + 13;
    if d == 1
        ind = cat(2,coastal_IDs, offshore_IDs);
    else
        ind = offshore_IDs;
    end


    lon_gauge = lon_gauge(:,ind);
    lat_gauge = lat_gauge(:,ind);
    dist_gauge = dist_gauge(:,ind);
    code_gauge = code_gauge(:,ind);
    angle_gauge = angle_gauge(:,ind);
    quadrant_gauge = quadrant_gauge(:,ind);
    phase_gauge = phase_gauge(:,ind);
    CP_gauge = CP_gauge(:,ind);


 
    loncp = reshape(lon_gauge,[],1);
    latcp = reshape(lat_gauge,[],1);
    cp_dist = reshape(dist_gauge,[],1);
    g_dist = reshape(dist_gauge,[],1);
    cp_quad = reshape(quadrant_gauge,[],1);
    cp_code = reshape(code_gauge,[],1);
    g_quad = reshape(quadrant_gauge,[],1);
    g_code = reshape(code_gauge,[],1);
    cp_flag = reshape(CP_gauge,[],1);
    d2land = reshape(repmat(d2land_center,1,size(lon_gauge,2)),[],1);
    
    Msize = 1000;
    filt0 = loncp >= LONR1 & loncp <= LONR2 & latcp >= LATR1 & latcp <= LATR2 & g_dist <= 500 & d2land > 0;
 
    
    [nt,nb] = size(lon_gauge);
    Vmax = reshape(repmat(Vmax_center,[1, nb]),[],1);
    Vtrans = reshape(repmat(Vtrans_center,[1, nb]),[],1);
    TCI = reshape(repmat(TCI_center,[1, nb]),[],1);

    cl1 = {my_color('dull light blue'),my_color('deep sky blue'),my_color('dark steel blue')};
    cl2 = {my_color('hot pink'),my_color('magenta'),my_color('violet')};
    bins = 0.5:1:4.5; % for quadrant count


    % plot CP vs distance

    for opt = 1:3
        if opt == 1
            filt1 = Vmax < 64;
            filt2 = Vmax >= 64;
            cats = {'TD-TS','Cat1-5'};
            tit = 'CPs vs. TC intensity';
        elseif opt == 2
            filt1 = Vtrans >= 6.2;
            filt2 = Vtrans < 6.2;
            cats = {'fast6.2','slow6.2'};
            tit = 'CPs vs. TC translation speed';
        elseif opt == 3
            filt1 = TCI <=0 ;
            filt2 = TCI > 0;   
            cats = {'weaken','strengthen'};
            tit = 'CPs vs. TC intensity change';
        end

        cp_quad1 = reshape(quadrant_gauge(filt0 & filt1 & cp_flag == 1),[],1);
        cp_quad2 = reshape(quadrant_gauge(filt0 & filt2 & cp_flag == 1),[],1);

        g_quad1 = reshape(quadrant_gauge(filt0 & filt1),[],1);
        g_quad2 = reshape(quadrant_gauge(filt0 & filt2),[],1);


        [freq1,CI1] = coldpool2_function_histogram_bootstrap_v2(cp_quad1,g_quad1,bins,Msize);
        err1 = CI1 - freq1; 

        [freq2,CI2] = coldpool2_function_histogram_bootstrap_v2(cp_quad2,g_quad2,bins,Msize);
        err2 = CI2 - freq2;    
       
        offset = 0.125; 
        width = 2*offset;
        X1 = [1:4] - offset;
        X2 = [1:4] + offset;
        h = [];
        sub=sub+1; subplot(nrow,ncol,order(sub)); hold on;
        h(1) = bar(X1,freq1,0.25,'FaceColor',cl1{opt})
        errorbar(X1,freq1,err1(:,1),err1(:,2),'k.','LineWidth',1.5);  

        h(2) = bar(X2,freq2,0.25,'FaceColor',cl2{opt})
        errorbar(X2,freq2,err2(:,1),err2(:,2),'k.','LineWidth',1.5);      
        legend(h(:),cats)

        ylabel('min.day^-^1'); xlim([0.5 4.5])  
        set(gca,'FontSize',14,'Xtick',1:4,'XtickLabel',{'FR','RR','RL','FL'});
        text(0.5,60,fignum{sub},'FontSize',16,'FontWeight','bold')
        title(tit)
        ylim([40,160]); set(gca,'FontSize',14,'Ygrid','on');
        
    end
    xlabel('TC quadrant');
    

end

set(gcf,'Color','w'); pause
export_fig([figout,'TC_quadrant_coldpool_distribution_vs_Vmax_Vtrans_TCI_coastal_vs_offshore.png'],'-dpng','-r200')



















































