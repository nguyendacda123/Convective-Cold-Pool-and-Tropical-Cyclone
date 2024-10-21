
clearvars; close all;
CP_general_parameters

figout = [figdir,'Cold_pool4/definition/']; 
if ~exist(figout,'dir'); mkdir(figout); end

load([matdir,'Cold_pools_grad_vs_diff_2018_2023.mat'],'TW_slope','time_slope','method','datname')
nd = length(datname);

%% plot histogram


figure('Position',[0, 0, 1000,900]);

nrow = 3; ncol = 2;
Tcrit = -0.7;
for k = 1:nd
    sub = 0;
    for i = 1:nrow
        dat = TW_slope(:,k,i);  
        sub = sub+1; subplot(nrow,ncol,sub); hold on;
        hist(dat,[-3:0.2:3]);     
        xlabel('Ta change');ylabel('count')
        grid on; xlim([-6,6]); % ylim([0,7e5]);
        text(-4,1e5,['Total: ',num2str(sum(isfinite(dat)))],'Color','r','FontSize',12)
        set(gca,'FontSize',12);
        title(method{i});

        sub = sub+1; subplot(nrow,ncol,sub); hold on;
        if k == 1
            filt = dat < Tcrit;
            hist(dat(filt),[-7:0.7:-0.7]); 
            xlim([-6,0]); %ylim([0 1800])
            text(-4,40,['Total: ',num2str(sum(filt))],'Color','r','FontSize',12)
            title([method{i},' < -0.7 ^oC']); 
        else
            filt = dat > -Tcrit;
            hist(dat(filt),[0.7:0.7:7]); 
            xlim([0,6]);
            text(4,40,['Total: ',num2str(sum(filt))],'Color','r','FontSize',12)
            title([method{i},' > 0.7 m.s^-^1']); 
        end
           
        xlabel([method{k},' change']);ylabel('count'); grid on;        
        set(gca,'FontSize',12);

    end
    set(gcf,'Color','w'); pause
    export_fig([figout,datname{k},'_distribution_grad_vs_diff_2021_2023.png'],'-dpng','-r200') 
    pause(0.1); clf
end

%% plot percentile pdf
figure('Position',[0 0 700 1000]);
nrow = 3; ncol = 2;

unit = {' (^oC.10min^-^1)',' (m.s^-^1.10min^-^1)'};

sub = 0;
for m = 1:length(method)
    dTa = squeeze(TW_slope(:,1,m));
    dWnd = squeeze(TW_slope(:,2,m));
    dQ = squeeze(TW_slope(:,3,m));
    
    filt = isfinite(dTa);
    dTa = dTa(filt);
    N = length(dTa);


    p = [0:0.1:1,1.5:0.5:5, 6:95, 95.5:0.5:99, 99.2:0.2:100];
    p_mid = 0.5*(p(2:end) + p(1:end-1));
    

    Tcrit = nan(length(p),1);
    for k = 1:length(p)
        Tcrit(k) = prctile(dTa,p(k));
    end
    
    Wcrit = nan(length(p_mid),1);
    Qcrit = nan(length(p_mid),1);
    for k = 1:length(p_mid)
        filt = dTa >= Tcrit(k) & dTa <= Tcrit(k+1);
        Wcrit(k) = nanmedian(dWnd(filt));
        Qcrit(k) = nanmedian(dQ(filt));

    end
        

    x1 = [0, 0, 90];
    x2 = [100, 10, 100];
    for k = 1:ncol
        sub = sub+1; subplot(nrow,ncol,sub); hold on;

        plot(p,Tcrit,'b-','LineWidth',2)
        xlabel('Percentile'); 
        ylabel(['dTa',unit{1}])  
        xlim([x1(k),x2(k)])
        
        if k == 1
            set(gca,'Xtick',0:20:100); xlim([-1 101]); ylim([-7, 5]);
            title(method{m});
            text(5,Tcrit(p == 5),fignum(sub),'FontWeight','bold','FontSize',14,'Color','k')
        elseif k == 2
            pc1 = Tcrit(p==0.5);
            plot(p,ones(length(p),1)*pc1,'--k');
            text(2,pc1-0.5,['lowest 1%: ',num2str(pc1,'%.2f')],'FontSize',12,'Color','k','FontWeight','bold');
            set(gca,'Xtick',0:10); ylim([-7, 0]);
            title([method{m},' (0-10%)'])
            text(2,Tcrit(p == 2),fignum(sub),'FontWeight','bold','FontSize',14,'Color','k')
        else
            yyaxis left;
            plot(p_mid,Wcrit,'k-');
            yyaxis right;
            plot(p_mid,Qcrit,'k-');    
            set(gca,'Xtick',0:10)
            title([method{m},' Wind & humidity'])
        end
        grid on;
        set(gca,'FontSize',12);

    end
end
%     axes( 'Position', [0, 0.96, 1, 0.04] ) ;
%     set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
%     text( 0.5, 0, ['Distribution of rate of change of ',datname{d}], 'FontSize', 14', 'FontWeight', 'Bold', ...
%       'HorizontalAlignment', 'Center', 'VerticalAlignment', 'bottom','Color','k' ) ;      

set(gcf,'Color','w'); pause
export_fig([figout,datname{d},'_slope_distribution_v3.png'],'-dpng','-r200');
pause(0.1); clf;



%% plot scatter
load([matdir,'Cold_pools_grad_vs_diff_2018_2023.mat'],'TW_slope','time_slope','method','datname')
nd = length(datname);

Tcrits = [-0.70, -0.57, -0.65];
% Tcrit = [-0.68, -0.68, -0.68];
figure('Position',[0, 0, 1000,400]);
nrow = 1; ncol = 2;
sub = 0;
green = my_color('dark green');
violet = my_color('violet');
cl = {green,violet};
for i = 1:ncol 
    T_slope = squeeze(TW_slope(:,1,:));
       
    h = [];
    k1 = 1; 
    k2 = i+1;
    filt0 = isfinite(T_slope(:,k1) + T_slope(:,k2));
    dat1 = T_slope(filt0,k1);
    dat2 = T_slope(filt0,k2);


    % plot diff 10
    sub = sub+1; subplot(nrow,ncol,sub); hold on;
    h(1) = plot(dat1,dat2,'.','Color',cl{i}); lsline 
    xlabel([method{k1},' (^oC.10min^-^1)']); xlim([-7,3])
    ylabel([method{k2},' (^oC.10min^-^1)']); ylim([-7,3])

    [R2, RSE, R, pval, bias] = function_R2_RSE_R_bias(dat1,dat2);
    
    
    text(-6,-6,{['R = ',num2str(R,'%.2f'), ', pval = ',num2str(pval,'%.2f')];
             ['bias = ',num2str(bias,'%.2f')]},'FontSize',12); 
    set(gca,'FontSize',12); grid on;
    
    
end
set(gcf,'Color','w'); pause
export_fig([figout,'Scatter_grad_vs_diff.png'],'-dpng','-r300') 
    


    







































