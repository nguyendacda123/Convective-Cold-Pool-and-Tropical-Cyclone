
clearvars; close all;
CP_general_parameters

datadir = [orig_dir, 'NDBC_buoys/'];
figout = [figdir,'cold_pool_R1/land_effect/'];
if ~exist(figout,'dir'); mkdir(figout); end

%% compute Ta slopes
flist = dir([datadir,'*.nc']);
N0 = length(flist);
lons = nan(N0,1);
lats = nan(N0,1);
res = nan(N0,1);
minres = nan(N0,1);
tlength = nan(N0,1);
Ymin = nan(N0,1);
Mmin = nan(N0,1);
Dmin = nan(N0,1);
Ymax = nan(N0,1);
Mmax = nan(N0,1);
Dmax = nan(N0,1);
sname = cell(N0,1);


date0 = datenum(1970,1,1);
%%
p = gcp('nocreate');
if isempty(p)
    parpool(4)
end

time_ref = nan(N0,1);
parfor i = 1:N0
    i
    fname = [datadir,flist(i).name];
%     info = ncinfo(fname);
%     varlist = {};
%     for k = 1:length(info.Variables)
%         varlist= [varlist, info.Variables(k).Name];
%     end
    lons(i) = double(ncread(fname,'longitude'));
    lats(i) = double(ncread(fname,'latitude'));
    tvec = double(ncread(fname,'time'));
    
    t0 = ncreadatt(fname,'time','units')
    time_ref(i) = datenum(t0(15:end-4));     
    

    if (tvec(end) - tvec(1))/86400 >= 31 % at least 1 month of data
        [Ymin(i),Mmin(i),Dmin(i)] = datevec(date0+nanmin(tvec)/86400);
        [Ymax(i),Mmax(i),Dmax(i)] = datevec(date0+nanmax(tvec)/86400); 
        dt = tvec(2:end) - tvec(1:end-1);
        res(i) = mode(dt)/60;
        tlength(i) = datenum(Ymax(i),Mmax(i),Dmax(i)) - datenum(Ymin(i),Mmin(i),Dmin(i));
        sname{i} = fname(end-12:end-8);
    end
    
end
save([matdir,'NDBC_buoys_metadata.mat'],'Ymin','Mmin','Dmin','Ymax','Mmax','Dmax','sname','res','tlength','lons','lats','time_ref')




%% plot map of locations
load([matdir,'NDBC_buoys_metadata.mat'],'Ymin','Mmin','Dmin','Ymax','Mmax','Dmax','sname','res','tlength','lons','lats')



figure('Position',[0 0 1200 800]); hold on;
plot_projection_ROI(-100,-48,12,35);
mk = {'rs','bo','m^'};
leg = {'10min','30min','60min'};
r = [10,30,60];
h = [];
for k = 1%:length(r)
    
%     filt = res == r(k);
    filt = res <= 10 & lons >= -100 & lons <= 0 & lats >= 0 & lats <= 35;
    h(k) = m_scatter(lons(filt),lats(filt),60,mk{k},'filled');
    
    

    
end
legend(h(:),leg,'Location','Best','FontSize',14)
title('Location & time resolution of NDBC buoys')
set(gcf,'Color','w'); pause
export_fig([figout,'NDBC_buoys_locations_v2.png'],'-dpng','-r300')
close all;


%% making sort list (resolution of 10 min or less)
load([matdir,'NDBC_buoys_metadata.mat'],'Ymin','Mmin','Dmin','Ymax','Mmax','Dmax','sname','res','tlength','lons','lats')
ind = find(res <= 10 & lons >= -100 & lons <= 0 & lats >= 0 & lats <= 35);
sname_10 = sname(ind);
lons_10 = lons(ind);
lats_10 = lats(ind);

figure('Position',[0 0 1200 800]); hold on;
plot_projection_ROI(-100,-48,12,35);
filt = nan(length(lons_10),1);

for i = 1:length(lons_10)
    
    if i > 1
        m_scatter(lons_10(1:i-1),lats_10(1:i-1),60,'bo','filled');
    end
    m_scatter(lons_10(i),lats_10(i),60,'ro','filled');
    filt(i) = input('press 1 = yes; 0 = no');
end

clf
plot_projection_ROI(-100,-48,12,35); hold on;
m_scatter(lons_10(filt == 1),lats_10(filt ==1),60,'ro','filled');
title('NDBC buoys filtered','FontSize',16); set(gcf,'Color','w')
export_fig([figout,'NDBC_buoys_locations_v3_filtered.png'],'-dpng','-r300')

offshore_flag = filt;
save([matdir,'NDBC_buoys_metadata.mat'],'offshore_flag','-append')

save([matdir,'NDBC_buoys_cold_pools_offshore.mat'],'lons_10','lats_10','sname_10')



%% process NDBC raw data

datadir = [orig_dir, 'NDBC_buoys/'];

load([matdir,'NDBC_buoys_cold_pools_offshore.mat'],'lons_10','lats_10','sname_10')
ref_date = datenum(1970,1,1);
N_gtime = [];
N_gTa = [];
N_gID = [];
N_glon = [];
N_glat = [];
N_SST = [];
N_RH = [];
N_WSPD = [];

figure('Position',[0 0 1200 1200]);
nrow = 4; ncol = 1;
verif = 0;
for s = 1:length(sname_10)
    s   
    fname = [datadir, sname_10{s},'h9999.nc']
    t = squeeze(double(ncread(fname,'time')))/86400 + ref_date;
    N_gtime = cat(1,N_gtime,t);
    
    SST = squeeze(double(ncread(fname,'sea_surface_temperature')));
    Wspd = squeeze(double(ncread(fname,'wind_spd')));
    Ta = squeeze(double(ncread(fname,'air_temperature')));
    Td = squeeze(double(ncread(fname,'dewpt_temperature'))); 
    
    RH = 100-5*(Ta-Td);
    

    N_gID = cat(1,N_gID,ones(length(t),1)*s);
    N_glon = cat(1,N_glon,ones(length(t),1)*lons_10(s));
    N_glat = cat(1,N_glat,ones(length(t),1)*lats_10(s));

    
    datp = cat(2,Ta,Wspd,SST,RH);
    datname = {'Ta','Wspd','SST','Rh'};
    unit = {'^oC','m.s^-^1','^oC','%'};             
    if verif == 1 % visual inspection for Ta, Wspd, SST, Rh

        for k = 1:length(datname)

            subplot(nrow,ncol,k); hold on;
%             scatter(t,datp(:,k),2,'.');
            plot(t,datp(:,k),'-');
            xtick = linspace(min(t),max(t),10);
            xtickl = datestr(xtick);
            set(gca,'Xtick',xtick,'XtickLabel',xtickl(:,1:11),'FontSize',12);
            ylabel(unit{k});
            title(['NDBC - ',sname_10{s},' - ',datname{k}])
            text(xtick(1),nanmax(datp(:,k)),fignum{k},'FontWeight','bold','FontSize',14);
            edit_flag = input('edit ? 1=yes, 0=no');
            if edit_flag == 1
                bound = input('left, right, upper, below');
                if bound{1} ~= 'n'
                    t_cutoff = datenum(bound{1},'ddmmyyyy');
                    datp(t<t_cutoff,k) = nan;
                end
                if bound{2} ~= 'n'
                    t_cutoff = datenum(bound{2},'ddmmyyyy');
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
        export_fig([figout,'NDBC_',sname_10{s},'.png'],'-dpng','-r200');
        pause(0.1); clf

    end 
    N_gTa = cat(1,N_gTa,datp(:,1));
    N_WSPD = cat(1,N_WSPD,datp(:,2));
    N_SST = cat(1,N_SST,datp(:,3)); 
    N_RH = cat(1,N_RH,datp(:,4));    
    
    
end

save([matdir,'NDBC_offshore_buoys_raw_data.mat'],'N_gtime','N_gTa','N_RH',...
    'N_gID','N_glon','N_glat','N_SST','N_WSPD','lons_10','lats_10','sname_10','-append')



%% label buoys that are too close to land

load([matdir,'NDBC_offshore_buoys_raw_data.mat'],'N_gtime','N_gTa','N_RH',...
    'N_gID','N_glon','N_glat','N_SST','N_WSPD','lons_10','lats_10','sname_10')
lons_10 = lons_10 + 360;
figure; hold on;
plot_projection_NWAtlantic
land_IDs = [];
ocean_IDs = [];
for i = 1:length(sname_10)
    m_scatter(lons_10(i),lats_10(i),'ro')
    m_text(lons_10(i),lats_10(i),sname_10{i});
    answer = input('Remove this one: 1=yes, 0=no?');
    if answer == 1
        land_IDs = [land_IDs, i];
    else
        ocean_IDs = [ocean_IDs,i];
    end
end

save([matdir,'NDBC_offshore_buoys_raw_data.mat'],'land_IDs','ocean_IDs','-append')


%% plot

figure('Position',[0 0 1200 800]); hold on;
plot_projection_ROI(260,315,8,40)
h(1) = m_scatter(lons_10(land_IDs),lats_10(land_IDs),60,'ro','filled');
h(2) = m_scatter(lons_10(ocean_IDs),lats_10(ocean_IDs),60,'bo','filled');

legend(h(:),{'coastal','offshore'},'Location','Best','FontSize',14)
title('NDBC buoys')
set(gcf,'Color','w'); pause
export_fig([figout,'NDBC_buoys_locations.png'],'-dpng','-r200')





