
% general parameters

homedir = '/media/da/DA5T/Cold_Pool/'; % working directory
orig_dir = [homedir,'Datasets/']; % original data directory
matdir = [homedir,'matdir/']; % processed data directory
figdir = [homedir,'figures/']; % output figure

fignum = {'a)', 'b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m','n'};


reftime = datenum(1970,1,1,0,0,0);
date0 = datenum(2020,1,1,0,0,0);

LONMIN = 250;
LONMAX = 380;
LATMIN = -10;
LATMAX = 70;

LONR1 = 360-98;
LONR2 = 360-34;
LATR1 = 11;
LATR2 = 37;

Tcrits = [-0.70, -0.57, -0.65]; % 1 percentile


% Tcrits = [-0.41, -0.33, -0.38]; % 2 percentile
% tag = '_2prctile';

% Tcrits = [-1.09, -0.85, -0.98]; % 0.5 percentile

