%load in data (example given)
%---------------------------

clearvars; close all;
CP_general_parameters

figout = [figdir,'Cold_pool2/definition/']; 
if ~exist(figout,'dir'); mkdir(figout); end


for opt = 1:3
    if opt == 1
        load([matdir,'Saildrone_rawdata.mat'],'time','date0','Rh','Pa','Ta','PAR','SWR','Wspd','SST','SST_skin')
        ws = Wspd;
        rh = Rh;
        at = Ta;
        sst = SST;
        bp = Pa;
        swr = SWR;
        time = time + date0;
        prefix = 'saildrone_';
    elseif opt == 2
        load([matdir,'PIRATA_data.mat'],'time','Ta','wspd','wdir','sst','rain','rh','sID','sname','lon','lat')
        ws = wspd;
        rh = rh;
        at = Ta;
        sst = sst;
        bp = nan(size(ws));
        swr = nan(size(ws));
        prefix = 'pirata_';
    elseif opt == 3      
        load([matdir,'NDBC_offshore_buoys_raw_data.mat'],'N_gtime','N_gTa','N_RH',...
        'N_gID','N_glon','N_glat','N_SST','N_WSPD') 
    
        ws = N_WSPD;
        rh = N_RH;
        at = N_gTa;
        sst = N_SST;
        bp = nan(size(ws));
        swr = nan(size(ws));    
        time = N_gtime;
        prefix = 'ndbc_';
    end

    %%

    % tday=ncread([fdir sdname],'time')/60/60/24 + datenum(1970,1,1); %day number
    % ws=ncread([fdir sdname],'WIND_SPEED_MEAN');
    % rh=ncread([fdir sdname],'RH_MEAN');
    % at=ncread([fdir sdname],'TEMP_AIR_MEAN');
    % sst=ncread([fdir sdname],'TEMP_DEPTH_HALFMETER_MEAN');
    % bp=ncread([fdir sdname],'BARO_PRES_MEAN');
    %---------------------------

    lwr=370*ones(1,length(sst));    %downward longwave radiation (does not affect results)
    pr=0*ones(1,length(sst));       %precipitation (does not affect results)
    bp(isnan(bp))=1015;             %barometric pressure (does not affect results)

    %time variable (does not affect results)
    % for i=1:length(sst)
    % 	t=datevec(tday(i));
    % 	time(i)=t(1)*10^8 + t(2)*10^6 + t(3)*10^4 + t(4)*10^2 + t(5)*10^0;
    % end

    %Input to COARE algorithm:

    %     u = relative wind speed (m/s) at height zu(m)
    %     zu 
    %     t = bulk air temperature (degC) at height zt(m)
    %     zt 
    %    rh = relative humidity (%) at height zq(m)
    %     zq 
    %     P = surface air pressure (mb) (default = 1015)
    %    ts = water temperature (degC) see jcool below
    %    Rs = downward shortwave radiation (W/m^2) (default = 150) 
    %    Rl = downward longwave radiation (W/m^2) (default = 370)
    %   lat = latitude (default = +45 N)
    %    zi = PBL height (m) (default = 600m)
    %  rain = rain rate (mm/hr)
    %    cp = phase speed of dominant waves (m/s)  
    %  sigH =  significant wave height (m)

    %put input data into array
    data=-9999*ones(length(sst),15);
    data(:,1)=ws';
    data(:,2)=3.4;
    data(:,3)=at';
    data(:,4)=2.3;
    data(:,5)=rh';
    data(:,6)=2.3;
    data(:,7)=bp';
    data(:,8)=sst';
    data(:,9)=swr';
    data(:,10)=lwr';
    data(:,11)=20;
    data(:,12)=600;
    data(:,13)=pr';
    data(:,14)=nan;
    data(:,15)=nan;

    % %write COARE input data to file (example given)
    % %----------------------------------------------
    % pth=[matdir,'HF/'];
    % if ~exist(pth,'dir'); mkdir(pth); end
    % fid=fopen([pth 'lhf_input_saildrone.txt'],'w');
    % fprintf(fid,'%17.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %6.2f %6.2f %6.2f %6.2f %5.2f \r\n',data');
    % fclose(fid);
    % %----------------------------------------------
    % 
    % %%
    % fname = 'lhf_input_saildrone.txt';
    dat = coldpool3_HF_coare35vn(data);

    lhflx = dat(:,4);
    shflx = dat(:,3);

    save([matdir,'HF/',prefix,'HF.mat'],'lhflx','shflx','time')
end
