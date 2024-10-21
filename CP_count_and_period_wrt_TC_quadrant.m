clearvars; close all;
CP_general_parameters

% tag = '';
tag = '_2prctile';
% tag = '_0.5prctile';

figout = [figdir,'Cold_pool4/Maps_in_TCs',tag,'/']; 
if ~exist(figout,'dir'); mkdir(figout); end
                
load([matdir,'Active_gauges_and_CPs_in_TC_coordinates',tag,'.mat'],'lon_center','lat_center',...
    'time_center','Vmax_center','Vtrans_center','TCI_center','lon_gauge','lat_gauge','dist_gauge',...
    'code_gauge','angle_gauge','quadrant_gauge','TC_code','phase_gauge','CP_gauge','CP_presst',...
    'mdist_500','mdist_200','Dist_min','tDist_min','sst_gauge','shflx_gauge','lhflx_gauge',...
    'R34_center','R50_center','R64_center','Rmax_center','USA_R_meaning')   

%% identify CP period

dist = reshape(dist_gauge,[],1);
quad = reshape(quadrant_gauge,[],1);
time = reshape(repmat(time_center,1,size(lon_gauge,2)),[],1);
CP_flag = reshape(CP_gauge,[],1);
code = reshape(code_gauge,[],1);

filt = CP_flag == 1 & dist <= 500;
dist = dist(filt);
quad = quad(filt);
time = time(filt);
code = code(filt);
cp_period = [];
cp_quad = [];
cp_code = [];
ncp = 0;
Dmax = nanmax(Vtrans_center)*600/1000;
k1 = 1;
for i = 1:length(dist)-1 
    k2 = i + 1;
    if (time(k2) - time(i))*24*60 > 10.1 || dist(k2)-dist(i) > Dmax || code(k2) ~= code(i)
        ncp = ncp + 1;
        cp_period(ncp) = (k2-k1)*10;       
        cp_quad(ncp) = mode(quad(k1:k2-1)); 
        cp_code(ncp) = unique(code(k1:k2-1));
        k1 = k2;
        if cp_code(ncp) == 2
            cp_period(ncp);
        end
        
    end
end

datname = {'PIRATA','NDBC','Saildrone','All'};
bins = 0.5:1:4.5; % for quadrant count
Msize = 100;

cl1 = {my_color('dull light blue'),my_color('deep sky blue'),my_color('orange')};
cl2 = {my_color('chartreuse'),my_color('dark green'),my_color('violet red')};

nbin = length(bins);

figure('Position',[0 0 800 600]); nrow = 2; ncol = 2;
for d = 1:length(datname)
    
    if d < 4
        codefilt = cp_code == d;        
    else
        codefilt = isfinite(cp_code);  
    end 
   
            
    [period, period_ci] = coldpool3_function_histogram_bootstrap(cp_period(codefilt),cp_quad(codefilt),bins,Msize);

    period_err = period_ci - period;
    subplot(nrow,ncol,d); hold on; 
    bar(1:4,period(1:4),'FaceColor',cl1{1})

    errorbar(1:4,period(1:4),period_err(1:4,1),period_err(1:4,2),'k.'); 

    ylabel('CP period (min)'); xlabel('TC quadrant');
    set(gca,'Xtick',1:4,'XtickLabel',{'FR','RR','RL','FL'},'FontSize',12)
    text(0.8,10,fignum{d},'FontSize',16,'FontWeight','bold')

    title(datname{d})
    
end

axes( 'Position', [0, 0.96, 1, 0.04] ) ;
set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
text( 0.5, 0, 'Mean CP period in TC quadrants', 'FontSize', 20, 'FontWeight', 'Bold', ...
  'HorizontalAlignment', 'Center', 'VerticalAlignment', 'bottom','Color','r' ) ;
set(gcf,'Color','w');   pause
export_fig([figout,'CP_periods_in_TC_quadrants_barplot.png'],'-dpng','-r200')
clf;    
 

%% plot cp counts
figure('Position',[0 0 800 600]); nrow = 2; ncol = 2;
for d = 1:length(datname)
    
    if d < 4
        codefilt = cp_code == d;        
    else
        codefilt = isfinite(cp_code);  
    end 
    cp_freq = [];
    for q = 1:4
        cp_freq(q) = sum(cp_quad(codefilt) == q);
    end

    subplot(nrow,ncol,d); hold on; 
    bar(1:4,cp_freq(1:4),'FaceColor',cl1{1})

 

    ylabel('number of CPs'); xlabel('TC quadrant');
    set(gca,'Xtick',1:4,'XtickLabel',{'FR','RR','RL','FL'},'FontSize',12)
    text(0.8,10,fignum{d},'FontSize',16,'FontWeight','bold')

    title(datname{d})
    
end

axes( 'Position', [0, 0.96, 1, 0.04] ) ;
set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
text( 0.5, 0, 'Mean CP period in TC quadrants', 'FontSize', 20, 'FontWeight', 'Bold', ...
  'HorizontalAlignment', 'Center', 'VerticalAlignment', 'bottom','Color','r' ) ;
set(gcf,'Color','w');   pause
export_fig([figout,'CP_frequency_in_TC_quadrants_barplot.png'],'-dpng','-r200')
clf;  