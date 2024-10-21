

clearvars; close all;
CP_general_parameters

figout = [figdir,'Cold_pool4/Maps_in_TCs/']; 
if ~exist(figout,'dir'); mkdir(figout); end

%% test TC coordinate computation

Xc1 = [-1, -1, 1, 1];
Yc1 = [-1, 1, -1, 1];

Xc2 = [-2, -2, 2, 2];
Yc2 = [-4, 4, -4, 4];
cname = {'SW','NW','SE','NE'};
coord_name = {'Distance to TC center','Angle relative to TC direction','TC quadrant'};
Xp = randsample(-500:1:500,1000);
Yp = randsample(-500:1:500,1000);

for opt = 3
    sub = 0;
    figure('Position',[0 0 800 400]); nrow = 1; ncol = 2;
    for k = [2,4]%1:length(Xc1)
        center1 = [Xc1(k), Yc1(k)];
        center2 = [Xc2(k),Yc2(k)];
        quad = nan(length(Xp),1);
        R =  nan(length(Xp),1);
        phi =  nan(length(Xp),1);
        Np = [center1(1), center1(2) + 2]; % North point
        for i = 1:length(Xp)
    %         quad(i) = HR_function_TC_quadrant(center1,center2,[Xp(i),Yp(i)]);
            [R(i), phi(i), quad(i)] = HR_function_TC_coordinates(center1,center2,[Xp(i),Yp(i)]); 
        end

%         subplot(nrow,ncol,k); hold on;
        a = (center2(2)-center1(2))/(center2(1)-center1(1));
        b = -1;

        XM = [-5,5];
        YM = b/a*(XM-center1(1)) + center1(2);  

        sub = sub+1; subplot(nrow,ncol,sub); hold on;
        if opt == 1
            scatter(Xp,Yp,6,R,'filled');
            cb = colorbar; xlabel(cb,'Distance (km)');  
            colormap(jet)
        elseif opt == 2
            scatter(Xp,Yp,6,phi,'filled');
            cb = colorbar; xlabel(cb,'Angle (deg)');
            colormap(my_cmap('summer_spring'))
        else
            scatter(Xp,Yp,10,quad,'o','filled');
%             cb = colorbar('yLim',[1, 4],'Ytick',1:4); xlabel(cb,'Quadrant'); 
            colormap(my_cmap('parula_autumn'))
        end
        
    %     scatter(Xp,Yp,3,quad,'filled');
    
%         plot([center1(1),center2(1)],[center1(2),center2(2)],'k-','LineWidth',1.5);
        U = center2(1)-center1(1);
        V = center2(2)-center1(2);
        h = quiver(center1(1),center1(2),U,V,'m-','LineWidth',1.5);
%         h.AutoScale = 0;
%         h.MaxHeadSize = 3;
        
        plot(XM,YM,'k--','LineWidth',1.5);
        U = Np(1)-center1(1);
        V = Np(2)-center1(2);
%         plot([center1(1),Np(1)],[center1(2),Np(2)],'m-','LineWidth',1.5);
        h = quiver(center1(1),center1(2),U,V,'m-','LineWidth',1.5);
%         h.AutoScale = 0;
%         h.MaxHeadSize = 3;
        text(Np(1),Np(2),'True North','Color','m','FontSize',12,'BackgroundColor','w');
            

        text(center1(1),center1(2),' center(t)','Color','m','FontSize',12,'BackgroundColor','w');
        text(center2(1),center2(2),' center(t+dt)','Color','m','FontSize',12,'BackgroundColor','w');
        text([4,4,-4,-4],[4,-4,-4,4],{'FR','RR','RL','FL'},'FontSize',12,'FontWeight','bold')
        xlim([-5, 5])
        ylim([-5, 5])
        title([cname{k}, ' storm'])
        set(gca,'FontSize',12);


    end
%     axes( 'Position', [0, 0.96, 1, 0.04] ) ;
%     set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
%     text( 0.5, 0, [coord_name{opt}], 'FontSize', 16', 'FontWeight', 'Bold', ...
%       'HorizontalAlignment', 'Center', 'VerticalAlignment', 'bottom','Color','k' ) ;     

    set(gcf,'Color','w'); pause
    if opt == 1
        fname = 'TC_distance_example.png';
    elseif opt == 2
        fname = 'TC_angle_example.png';
    else
        fname = 'TC_quadrant_example.png';
    end
    
    export_fig([figout,fname],'-dpng','-r200')
    
    
end

       
                















































