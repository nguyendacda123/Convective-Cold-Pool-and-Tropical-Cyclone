
clearvars; close all;
CP_general_parameters

datadir = [orig_dir, 'PIRATA_buoys/'];
figout = [figdir,'Cold_pool_R1/Quality_control/']; 
if ~exist(figout,'dir'); mkdir(figout); end

%% compute Ta slopes
prefix = {'airt','w','sst','rain','rh'}; % 4 filename
varlist = {'AT_21','WS_401','WD_410','T_25','RN_485','RH_910'}; % 7 variables
varname = {'Ta','wspd','wdir','sst','rain','Rh'};
unit = {'^oC','m.s^-^1','deg','^oC','mm.h^-^1','%'};
Qlist = {'QAT_5021','QWS_5401','QWD_5410','QT_5025','QRN_5485','QRH_5910'};
vcode = [1, 2, 2, 3, 4, 5];

flist = dir([datadir,'sst*_10m.cdf']);
station_id = {};
lon = [];
lat = [];
for k = 1:length(flist)
    fname = flist(k).name;
    station_id{k} = fname(4:end-8);
    lon(k) = double(ncread([datadir,fname],'lon'));
    lat(k) = double(ncread([datadir,fname],'lat'));
    tref = ncreadatt([datadir,fname],'time','units');
end
ref_time = datenum(1998,11,8,22,50,00);

filt = lon >= LONR1 & lon <= LONR2 & lat >= LATR1 & lat <= LATR2;
station_id = station_id(filt);
lon = lon(filt);
lat = lat(filt);

ns = length(station_id);


%%
data = [];
datcount = nan(ns,length(varlist),6);
datmean = nan(ns,length(varlist),6);
verif = 0;
vplot = [1,2,4,6];
figure('Position',[0 0 1200 1200]);
nrow = 4; ncol = 1;
for s = 1:ns
    dat = [];
    sub = 0; 
    flag = [];
    for v = 1:length(varlist)
        k = vcode(v);
        fname = [datadir, prefix{k},station_id{s},'_10m.cdf'];
        
        t = double(squeeze(ncread(fname,'time')))/(24*60);
        tref = ncreadatt(fname,'time','units');
        t_v = t + datenum(tref(15:end)); 
        var = double(squeeze(ncread(fname,varlist{v})));
        flag = double(squeeze(ncread(fname,Qlist{v})));

        
        for k = 1:6
            datcount(s,v,k) = sum(flag == k-1);
            datmean(s,v,k) = nanmean(var(flag==k-1));
        end
        
        
        if v == 1
            tvec = t_v;
            dat_v = var;
            dat = cat(2,dat,tvec,dat_v);
        else
            if length(t_v) ~= length(tvec) || sum(t_v - tvec) ~= 0 % variables do not have the same time vector
                dat_v = interp1(t_v,var,tvec); % interpolate other vars to T-air
%                 disp('interp')
            else
                dat_v = var;
            end            
            dat = cat(2,dat,dat_v);            
        end
        % verification
        if verif == 1 && sum(v == vplot) == 1
            sub = sub+1; subplot(nrow,ncol,sub); hold on;
            plot(tvec,dat_v);
            xtick = linspace(min(tvec),max(tvec),10);
            xtickl = datestr(xtick);
            set(gca,'Xtick',xtick,'XtickLabel',xtickl(:,4:11),'FontSize',12);
            ylabel(unit{v});
            title(['PIRATA - ',station_id{s},' - ',varname{v}])
            text(xtick(1),nanmax(dat_v),fignum{sub},'FontWeight','bold','FontSize',14);
            if sub == 4
                set(gcf,'Color','w'); pause
                export_fig([figout,'PIRATA_',station_id{s},'.png'],'-dpng','-r200');
                pause(0.1); clf
            end
        end
        
         
    end
    dat = cat(2,dat,ones(size(dat,1),1)*s);
    data = cat(1,data,dat);
end

time = data(:,1);
Ta = data(:,2);
wspd = data(:,3);
wdir = data(:,4);
sst = data(:,5);
rain = data(:,6);
rh = data(:,7);
sID = data(:,8);
sname = station_id;



% data filtering based on graph
Ta(sID == 1 & Ta < 20) = nan;
Ta(sID == 3 & Ta > 30) = nan;

    

save([matdir,'PIRATA_data.mat'],'time','Ta','wspd','wdir','sst','rain','rh','sID','sname','lon','lat')








%% data quality control verification

prefix = {'airt','w','sst','rain','rh'}; % 4 filename
varlist = {'AT_21','WS_401','WD_410','T_25','RN_485','RH_910'}; % 7 variables
qlist = {'QAT_5021','QWS_5401','QWD_5410','QT_5025','QRN_5485','QRH_5910'};
vcode = [1, 2, 2, 3, 4, 5];

flist = dir([datadir,'sst*_10m.cdf']);
station_id = {};
lon = [];
lat = [];
for k = 1:length(flist)
    fname = flist(k).name;
    station_id{k} = fname(4:end-8);
    lon(k) = double(ncread([datadir,fname],'lon'));
    lat(k) = double(ncread([datadir,fname],'lat'));
    tref = ncreadatt([datadir,fname],'time','units');
end
ref_time = datenum(1998,11,8,22,50,00);


ns = length(station_id);

data = [];

for s = 1:ns
    dat = [];
    for v = 1:length(varlist)
        k = vcode(v);
        fname = [datadir, prefix{k},station_id{s},'_10m.cdf'];
        
        t = double(squeeze(ncread(fname,'time')))/(24*60);
        tref = ncreadatt(fname,'time','units')
        t_v = t + datenum(tref(15:end)); 
        var = double(squeeze(ncread(fname,varlist{v})));
        Qual = double(squeeze(ncread(fname,qlist{v})));
        pause
    end
end
        








