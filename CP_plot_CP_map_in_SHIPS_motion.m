
clearvars; close all;
CP_general_parameters




figout = [figdir,'Cold_pool_R1/Maps_in_TCs/SHIPS_motion/']; 
if ~exist(figout,'dir'); mkdir(figout); end

% load data              
load([matdir,'Active_gauges_and_CPs_in_SHIPS_motion.mat'],'lon_center','lat_center',...
    'time_center','Vmax_center','SHRD_center','SHDC_center','lon_gauge','lat_gauge','dist_gauge',...
    'code_gauge','angle_gauge','quadrant_gauge','TC_code','phase_gauge','CP_gauge','CP_presst',...
    'mdist_500','mdist_200','Dist_min','tDist_min','sst_gauge','shflx_gauge','lhflx_gauge')

load([matdir,'NDBC_offshore_buoys_raw_data.mat'],'land_IDs','ocean_IDs');
land_IDs = land_IDs+13; % code of NDBC starting from 14
ocean_IDs = ocean_IDs + 13;

opt = 2;
if opt == 2
    figout = [figdir,'Cold_pool_R1/Maps_in_TCs/SHIPS_offshore/']; 
    if ~exist(figout,'dir'); mkdir(figout); end    
    
    ind = setdiff(1:size(lon_gauge,2),land_IDs)
    lon_gauge = lon_gauge(:,ind);
    lat_gauge = lat_gauge(:,ind);
    dist_gauge = dist_gauge(:,ind);
    code_gauge = code_gauge(:,ind);
    angle_gauge = angle_gauge(:,ind);
    quadrant_gauge = quadrant_gauge(:,ind);
    phase_gauge = phase_gauge(:,ind);
    CP_gauge = CP_gauge(:,ind);
    CP_presst = CP_presst(:,ind);
    mdist_500 = mdist_500(:,ind);
    mdist_200 = mdist_200(:,ind);
    Dist_min = Dist_min(:,ind);
    tDist_min = tDist_min(:,ind);
    sst_gauge = sst_gauge(:,ind);
    shflx_gauge = shflx_gauge(:,ind);
    lhflx_gauge = lhflx_gauge(:,ind);
    tag = '_motion_offshore';
else
    tag = '_motion';
end
    

%% plot against data type & period
loncp = reshape(lon_gauge,[],1);
latcp = reshape(lat_gauge,[],1);
R = abs(reshape(dist_gauge,[],1));
phi = reshape(angle_gauge,[],1);
type = reshape(code_gauge,[],1);
mind500 = reshape(mdist_500,[],1);
cp_flag = reshape(CP_gauge,[],1);
quad = reshape(quadrant_gauge,[],1);


filt0 = loncp >= LONR1 & loncp <= LONR2 & latcp >= LATR1 & latcp <= LATR2 & mind500 == 1 & cp_flag == 1;

R = R(filt0);
phi = phi(filt0);
type = type(filt0);
quad = quad(filt0);

phic = linspace(0,2*pi,721);
ind1 = find(phic == 0);
ind2 = find(phic == pi);
ind3 = find(phic == pi/2);
ind4 = find(phic == 3*pi/2);

cl = {'b','r','m'};
codename = {'PIRATA','NDBC','saildrone','All'};
figure('Position',[0 0 1100 1000]);
nrow = 2, ncol = 2;
for d = 1:length(codename)

    subplot(nrow,ncol,d); hold on;            
    filt1 = R <= 500;
    Rc = [100:200:500];
    
    if d < 4
        filt2 = type == d;
    else
        filt2 = isfinite(type);
    end
    filt = filt1 & filt2;


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
    text(-400,400,fignum{d},'FontSize',16,'FontWeight','bold')
    set(gca,'FontSize',14);
           
    title(codename{d})

end
set(gcf,'Color','w'); pause
export_fig([figout,'CP_locations_in_VWS',tag,'.png'],'-dpng','-r200')

%% plot histogram of quadrant data
datname = {'PIRATA','NDBC','Saildrone','All'};
datname2 = {'CP time','Sampling time', 'CP/Sampling time'};
cp_dist = reshape(dist_gauge,[],1);
g_dist = reshape(dist_gauge,[],1);
cp_quad = reshape(quadrant_gauge,[],1);
cp_code = reshape(code_gauge,[],1);
g_quad = reshape(quadrant_gauge,[],1);
g_code = reshape(code_gauge,[],1);
cp_flag = reshape(CP_gauge,[],1);
mind500 = reshape(mdist_500,[],1);

filt0 = loncp >= LONR1 & loncp <= LONR2 & latcp >= LATR1 & latcp <= LATR2 & mind500 == 1;

cp_dist = cp_dist(filt0 & cp_flag == 1);
g_dist = g_dist(filt0);
cp_quad = cp_quad(filt0 & cp_flag == 1);
cp_code = cp_code(filt0 & cp_flag == 1);
g_quad = g_quad(filt0);
g_code = g_code(filt0);



% TC_days = sum(filt0)*600/86400;

Rplot = [500, 200];

figure('Position',[0 0 1200 900]);
nrow = 4; ncol = 3;
subnum = [1 4 7 10 2 5 8 11 3 6 9 12];

bins = 0.5:1:4.5; % for quadrant count
Msize = 1000;
sub = 0;
cl = {my_color('dark steel blue'),my_color('dull light blue'),my_color('deep sky blue')};
tab = [];
for opt = 1:3
    
    for d = 1:length(datname)

        if d < 4
            cfilt1 = cp_code == d;
            cfilt2 = g_code == d;        
        else
            cfilt1 = isfinite(cp_code);
            cfilt2 = isfinite(g_code);  
        end    

        for i = 1%:length(Rplot)
            sub = sub+1; subplot(nrow,ncol,subnum(sub)); hold on;    


            dist_filt1 = abs(cp_dist) <= 500;
            dist_filt2 = abs(g_dist) <= 500;


            if opt == 1
                filt = dist_filt1 & cfilt1;
                dat = cp_quad(filt);
                [freq,CI] = coldpool2_function_histogram_bootstrap(dat,bins,Msize);
                freq = freq*10;
                CI = CI*10;
                err = CI - freq;                
                bar(1:4,freq,'FaceColor',cl{opt})
                errorbar(1:4,freq,err(:,1),err(:,2),'k.');  
                if i == 1
                    ylabel({['\fontsize{16}',datname{d}],'\fontsize{12}min'});
                end
                tab = cat(2,tab,freq);

            elseif opt == 2
                filt = dist_filt2 & cfilt2;
                dat = g_quad(filt);
                
                [freq,CI] = coldpool2_function_histogram_bootstrap(dat,bins,Msize);
                freq = freq*600/86400;
                CI = CI*600/86400;
                err = CI - freq;                
                bar(1:4,freq,'FaceColor',cl{opt})
                errorbar(1:4,freq,err(:,1),err(:,2),'k.');  
                if i == 1
                    ylabel('day');
                end
                tab = cat(2,tab,freq);
            else
                filt1 = dist_filt1 & cfilt1;
                filt2 = dist_filt2 & cfilt2;
                dat1 = cp_quad(filt1);
                dat2 = g_quad(filt2);
                
                [freq,CI] = coldpool2_function_histogram_bootstrap_v2(dat1,dat2,bins,Msize);
                
                
                err = CI - freq;                
                bar(1:4,freq,'FaceColor',cl{opt})
                errorbar(1:4,freq,err(:,1),err(:,2),'k.');  

                if i == 1
                    ylabel('min.day^-^1');
                end
                
                tab = cat(2,tab,freq);   
            end
            
            if d == 4
                xlabel('TC quadrant'); 
            end
            
            if d == 1
                title(datname2{opt},'FontWeight','normal','FontSize',16)
            end
            set(gca,'Xtick',1:4,'XtickLabel',{'FR','RR','RL','FL'},'FontSize',12)
            text(0.5,nanmax(freq),fignum{sub},'FontSize',16,'FontWeight','bold')
        end            
    end

    
end
 

set(gcf,'Color','w'); pause
export_fig([figout,'CP_in_VWS_quadrants_barplot',tag,'.png'],'-dpng','-r200')    
pause(0.1);
clf;


%% plot against Vmax, Vtrans, TCI
loncp = reshape(lon_gauge,[],1);
latcp = reshape(lat_gauge,[],1);
cp_dist = reshape(dist_gauge,[],1);
g_dist = reshape(dist_gauge,[],1);
cp_quad = reshape(quadrant_gauge,[],1);
cp_code = reshape(code_gauge,[],1);
g_quad = reshape(quadrant_gauge,[],1);
g_code = reshape(code_gauge,[],1);
cp_flag = reshape(CP_gauge,[],1);

Msize = 1000;
filt0 = loncp >= LONR1 & loncp <= LONR2 & latcp >= LATR1 & latcp <= LATR2 & g_dist <= 500;
[nt,nb] = size(lon_gauge);
Vmax = reshape(repmat(Vmax_center,[1, nb]),[],1);
SHRD = reshape(repmat(SHRD_center,[1, nb]),[],1)/10;% convert to knots
SHDC = reshape(repmat(SHDC_center,[1, nb]),[],1)/10;



cl1 = {my_color('dull light blue'),my_color('deep sky blue'),my_color('dark steel blue')};
cl2 = {my_color('hot pink'),my_color('magenta'),my_color('violet')};
bins = 0.5:1:4.5; % for quadrant count

figure('Position',[0 0 600 800]);
nrow = 3; ncol = 1;
% plot CP vs distance
sub = 0;
for opt = 1:3
    if opt == 1
        filt1 = Vmax < 64;
        filt2 = Vmax >= 64;
        leg = {'TD-TS','Cat1-5'};
        tit = 'CPs vs. TC intensity';
    elseif opt == 2
        filt1 = SHRD < prctile(SHRD,50);
        filt2 = SHRD >= prctile(SHRD,50);
        leg = {'<15.9kt','>15.9kt'};
        tit = 'CPs vs. SHRD';
    elseif opt == 3
        filt1 = SHDC < prctile(SHDC,50);
        filt2 = SHDC >= prctile(SHDC,50); 
        leg = {'<13.9kt','>13.9kt'};
        tit = 'CPs vs. SHDC';
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
    subplot(nrow,ncol,opt); hold on;
    h(1) = bar(X1,freq1,0.25,'FaceColor',cl1{opt})
    errorbar(X1,freq1,err1(:,1),err1(:,2),'k.','LineWidth',1.5);  
    
    h(2) = bar(X2,freq2,0.25,'FaceColor',cl2{opt})
    errorbar(X2,freq2,err2(:,1),err2(:,2),'k.','LineWidth',1.5);      
    legend(h(:),leg)
    
    ylabel('min.day^-^1'); xlim([0.5 4.5]) ; ylim([0, 150])
    set(gca,'FontSize',14,'Xtick',1:4,'XtickLabel',{'FR','RR','RL','FL'});
    text(0.5,60,fignum{opt},'FontSize',16,'FontWeight','bold')
    title(tit)
         
end
xlabel('TC quadrant');
set(gcf,'Color','w'); pause
export_fig([figout,'Coldpool_distribution_vs_Vmax_and_shear',tag,'.png'],'-dpng','-r200')



















































