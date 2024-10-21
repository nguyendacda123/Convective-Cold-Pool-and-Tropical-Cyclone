
clearvars; close all;
CP_general_parameters

datadir1 = [orig_dir,'Saildrone/sd_202310/']; % data updated on Oct 2023
datadir2 = [orig_dir,'Saildrone/sd_all_sfc/'];
figout = [figdir,'Quality_control/'];
if ~exist(figout,'dir'); mkdir(figout); end

%% General info
% SD_id = [31,32,36,40,41,45,48,57,59,60,64,65,68,69,78,83,84]; % number 31 doesn't have wind speed
SD_id = 1:84;

SD_Y = 2017:2023;
reftime = datenum(1970,1,1,0,0,0);


lon = [];
lat = [];
time = [];
sID = [];
Y_label = [];

% air params
Ta = [];
Rh = [];
Pa = [];
PAR = []; 
Wspd = [];
Wdir = [];
Uwind = [];
Vwind = [];
Gust = [];

% sea surface params
SST = []; % sea surface temp
SST_skin = []; % sea surface skin temperature
SSS = []; % sea surface salt
Wav_T = []; % Wave dominant period
Hs = [];
C_O2 = []; % Concentration of O2
S_O2 = []; % Oxygen saturation

% upper-layer params
Cspd = []; % current speed
Cdir = []; % current direction
CDom = [];
Chla = [];



figure('Position',[0 0 1200 1200]);
nrow = 4; ncol = 1;
verif = 0;

for id = SD_id
    for Y = 2019:2023
        
        flist1 = dir([datadir1,'sd10',num2str(id,'%.2d'),'*',num2str(Y),'*.nc']);
        flist2 = dir([datadir2,'sd10',num2str(id,'%.2d'),'*',num2str(Y),'*.nc']);
        if ~isempty(flist1) || ~isempty(flist2)
            if ~isempty(flist1)
                fname = [datadir1,flist1.name];
            else
                fname = [datadir2,flist2.name];
            end
            info = ncinfo(fname);
            varlist = {};
            for k = 1:length(info.Variables)
                varlist= [varlist, info.Variables(k).Name];
            end
            
            disp(['Processing: ',fname])
            
 
            

            lon = cat(1,lon,double(ncread(fname,'longitude')));
            lat = cat(1,lat,double(ncread(fname,'latitude')));
            t = double(ncread(fname,'time'))/86400 + reftime;
            time = cat(1,time,t);
            N = length(t);
            
            % air params
            
            Pa = cat(1,Pa,double(ncread(fname,'BARO_PRES_MEAN')));
            
%             Uwind = cat(1,Uwind,double(ncread(fname,'UWND_MEAN')));
%             Vwind = cat(1,Vwind,double(ncread(fname,'VWND_MEAN')));
%             Gust = cat(1,Gust,double(ncread(fname,'GUST_WND_MEAN')));
            Ta_p = double(ncread(fname,'TEMP_AIR_MEAN'));
            Rh_p = double(ncread(fname,'RH_MEAN'));
            
            
            % wind speed -----------------------
            if any(strcmp(varlist,'WIND_SPEED_MEAN'))           
                Wspd_p = double(ncread(fname,'WIND_SPEED_MEAN'));
            elseif any(strcmp(varlist,'wind_speed'))           
                Wspd_p = double(ncread(fname,'wind_speed'));
            elseif any(strcmp(varlist,'UWND_MEAN'))
                U = double(ncread(fname,'UWND_MEAN'));
                V = double(ncread(fname,'VWND_MEAN'));
                Wspd_p = sqrt(U.*U + V.*V);
            else              
                Wspd_p = nan(N,1);
            end
            
            % wind direction ------------------------
            if any(strcmp(varlist,'WIND_FROM_MEAN'))           
                Wdir = cat(1,Wdir,double(ncread(fname,'WIND_FROM_MEAN')));
            elseif any(strcmp(varlist,'wind_dir'))           
                Wdir = cat(1,Wdir,double(ncread(fname,'wind_dir')));                
            else
                Wdir = cat(1,Wdir,nan(N,1));    
            end

            
            % sea surface temperature --------------------
            if any(strcmp(varlist,'TEMP_CTD_MEAN'))               
                SST_p = double(ncread(fname,'TEMP_CTD_MEAN'));
            elseif any(strcmp(varlist,'TEMP_SBE37_MEAN')) 
                SST_p = double(ncread(fname,'TEMP_SBE37_MEAN'));
            else
                SST_p = nan(N,1);
            end
            
%             SST_skin = cat(1,SST_skin,double(ncread(fname,'TEMP_IR_SEA_WING_UNCOMP_MEAN')));
            % sea surface salinity --------------------
            if any(strcmp(varlist,'SAL_CTD_MEAN'))               
                SSS = cat(1,SSS,double(ncread(fname,'SAL_CTD_MEAN')));
            elseif any(strcmp(varlist,'SAL_SBE37_MEAN')) 
                SSS = cat(1,SSS,double(ncread(fname,'SAL_SBE37_MEAN')));
            else
                SSS = cat(1,SSS,nan(N,1));
            end
            
            % Wave period ------------------------
            if any(strcmp(varlist,'WAVE_DOMINANT_PERIOD'))           
                Wav_T = cat(1,Wav_T,double(ncread(fname,'WAVE_DOMINANT_PERIOD')));          
            else
                Wav_T = cat(1,Wav_T,nan(N,1));    
            end            
            
            % Wave height ------------------------
            if any(strcmp(varlist,'WAVE_SIGNIFICANT_HEIGHT'))           
                Hs = cat(1,Hs,double(ncread(fname,'WAVE_SIGNIFICANT_HEIGHT')));          
            else
                Hs = cat(1,Hs,nan(N,1));    
            end              

                       
%             C_O2 = cat(1,C_O2,double(ncread(fname,'O2_CONC_SBE37_MEAN'))); % unit = micromol/L
%             S_O2 = cat(1,S_O2,double(ncread(fname,'O2_SAT_SBE37_MEAN')));  % unit = % 
            
            
            % upper layer
            try
                PAR = cat(1,PAR,double(ncread(fname,'PAR_AIR_MEAN'))); % unit = micromol/s/m2
            catch
                PAR = cat(1,PAR,nan(N,1)); 
                disp(['PAR data not exist in :',fname])               
            end
            
            try
                Cspd = cat(1,Cspd,double(ncread(fname,'WATER_CURRENT_SPEED_MEAN')));
                Cdir = cat(1,Cdir,double(ncread(fname,'WATER_CURRENT_DIRECTION_MEAN')));
            catch
                Cspd = cat(1,Cspd,nan(N,1));
                Cdir = cat(1,Cdir,nan(N,1));  
                disp(['Current data not exist in :',fname])
                
            end   
            
            try 
                Chla = cat(1,Chla,double(ncread(fname,'CHLOR_WETLABS_MEAN'))); % unit = microgram/L
            catch
                Chla = cat(1,Chla,nan(N,1));
                disp(['Chla not exist in :',fname])
            end
            
            try 
                CDom = cat(1,CDom,double(ncread(fname,'CDOM_MEAN')));  % unit = ppb   

            catch
                CDom = cat(1,CDom,nan(N,1));
                disp(['Cdom not exist in :',fname])
            end
                                   
            label = str2num(['10',num2str(id,'%.2d'),num2str(Y)]);
            sID = cat(1,sID,ones(N,1)*label);
            datp = cat(2,Ta_p,Wspd_p,SST_p,Rh_p);
            datname = {'Ta','Wspd','SST','Rh'};
            unit = {'^oC','m.s^-^1','^oC','%'};             
            if verif == 1 % visual inspection for Ta, Wspd, SST, Rh
                               
                for k = 1:length(datname)
                    
                    subplot(nrow,ncol,k); hold on;
                    scatter(t,datp(:,k),3,'.');
                    xtick = linspace(min(t),max(t),10);
                    xtickl = datestr(xtick);
                    set(gca,'Xtick',xtick,'XtickLabel',xtickl(:,1:6),'FontSize',12);
                    ylabel(unit{k});
                    title(['SD10',num2str(id,'%.2d'),'-',num2str(Y),' - ',datname{k}])
                    text(xtick(1),nanmax(datp(:,k)),fignum{k},'FontWeight','bold','FontSize',14);
                    edit_flag = input('edit ? 1=yes, 0=no');
                    if edit_flag == 1
                        bound = input('left, right, upper, below');
                        if bound{1} ~= 'n'
                            t_cutoff = datenum([bound{1},num2str(Y)],'ddmmyyyy');
                            datp(t<t_cutoff,k) = nan;
                        end
                        if bound{2} ~= 'n'
                            t_cutoff = datenum([bound{2},num2str(Y)],'ddmmyyyy');
                            datp(t > t_cutoff,k) = nan;
                        end 
                        if bound{3} ~= 'n'
                            datp(datp(:,k) > bound{3},k) = nan;
                        end  
                        if bound{4} ~= 'n'
                            datp(datp(:,k) < bound{3},k) = nan;
                        end 
                    end
                        
                end

                set(gcf,'Color','w'); 
                export_fig([figout,'sd10',num2str(id,'%.2d'),'_',num2str(Y),'.png'],'-dpng','-r200');
                pause(0.1); clf

            end 
            Ta = cat(1,Ta,datp(:,1));
            Wspd = cat(1,Wspd,datp(:,2));
            SST = cat(1,SST,datp(:,3)); 
            Rh = cat(1,Rh,datp(:,4));
            
        end
    end
end
SWR =2*PAR/4.57; % short wave radiation from PAR: micromol/s/m2 -> J/s/m2 ~ w/m2

save([matdir,'Saildrone_rawdata.mat'],'lon','lat','time','date0','sID','Y_label','SD_id','SD_Y',...
                                       'Ta','Rh','Pa','PAR','SWR','Wspd','Wdir','Uwind','Vwind','Gust',...
                                       'SST','SST_skin','SSS','Wav_T','Hs','C_O2','S_O2', ...
                                       'Cspd','Cdir','Chla','CDom')

%% plot locations
load([matdir,'Saildrone_rawdata.mat'],'lon','lat','time','date0','Ta','Rh','Pa','SST',...
    'SSS','Cspd','Cdir','Wav_T','Hs','Wspd','Wdir','sID','Y_label','SD_id','SD_Y')

figure('Position',[0 0 1000 800]); hold on;
plot_projection_ROI(nanmin(lon),nanmax(lon),nanmin(lat),nanmax(lat))


idlist = unique(sID);
cl = jet(length(idlist));
for i = 1:length(idlist)
    id = idlist(i);


    filt = sID == id;
    lon1 = lon(filt); lat1 = lat(filt);
    m_plot(lon1,lat1,'.','Color',cl(i,:));
    k = find(isfinite(lon1),1,'first');
    k2 = find(isfinite(lon1),1,'last');
    m_plot(lon1(k),lat1(k),'ko','MarkerFaceColor','k');
    m_plot(lon1(k2),lat1(k2),'ks','MarkerFaceColor','k');
    m_text(nanmean(lon1),nanmean(lat1),['SD', num2str(id)],'FontSize',12,'Color',cl(i,:));

end

export_fig([figout,'SD_tracks.png'],'-dpng','-r300')


%% plot different variables
varname = {'Pa','Rh','Ta','SST','SSS','Hs','Wspd','Cspd','Wdir','Cdir'};
unit = {'hPa','%','^oC','^oC','psu','m','m.s^-^1','m.s^-^1','deg','deg'};
data = cat(2,Pa,Rh,Ta,SST,SSS,Hs,Wspd,Cspd,Wdir,Cdir);
sub = [1 3 5 7 9 2 4 6 8 10]; 

time_tick = [];
time_label = {};
for Y = SD_Y
    time_label = [time_label,'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'];
    time_tick = cat(1,time_tick,[datenum(Y,1:12,1)-date0]');
end

nvar = length(varname);
cl = jet(nvar);

t1 = datenum(2022,9,21,11,0,0) - date0;
t2 = datenum(2022,9,23,5,0,0) - date0;
nrow = 5; ncol = 2;
IDs = unique(sID)
for i = 1:length(IDs)

    for Y = SD_Y
        filt = sID == IDs(i);
        if sum(filt) > 0
            figure('Position',[0 0 1500 1600]);
            dat = data(filt,:);
            t = time(filt);
            for v = 1:nvar
                subplot(nrow,ncol,sub(v)); hold on;
                filt2 = isfinite(dat(:,v) + t);
                plot(t(filt2),dat(filt2,v),'-', 'Color',cl(v,:));
                legend(varname{v},'Location','best');
                ylabel(unit{v})
                set(gca,'FontSize',12,'Color',[1 1 1]*0.8);
                set(gca,'Xtick',time_tick,'XtickLabel',time_label); 
                xlim([t1, t2]);%xlim([min(t),max(t)])
                
            end
            
            idstr = num2str(id);
            axes( 'Position', [0, 0.95, 1, 0.05] ) ;
            set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
            text( 0.5, 0, ['SD',idstr(1:4),' - ',idstr(5:end)], 'FontSize', 14', 'FontWeight', 'Bold', ...
              'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ; 
            set(gcf,'Color','w')
            export_fig([figout,'SD',idstr,'_data_skimming.png'],'-dpng','-r300');
         
            pause
        end

    end

end
            
    




















































