function data = get_sctaDay2(baseDir,dateN)
% Gets a days worth of SCTA data with options to only get specific channels
%
% Usage
%   data = get_sctaDay(date,baseDir)
%
% Inputs
%   BASEDIR  - EITHER the base directory for where to fine the file so that
%              '/Users/wilcock/test/2017-08-02/SCTA-Accel_2017-08-24.nc'
%              would have a base directory of 
%              '/Users/wilcock/test'
%              OR the actual file name in which case DATEN is ignored
%   DATEN    - date number or string for day to get ifs BASEDIR is the base
%              directory and not the file to load
%
% Outputs
%   DATA    - Stucture with data in column vectors/matrices of the following fields
%               t - Time for acceleration 1
%               a  - Acceleration in 3 columns of x, y, z
%               as - Scalar acceleration (rss or root of sum of squares);
%               T - Temperture of the accelerometer 1

%% Preliminaries
if nargin~=2 && nargin~=1
  error('get_sctaDay requires 1-2 input arguements')
end

if isdir(baseDir)
  if ~strcmp(baseDir(end),'/')
    baseDir = [baseDir '/'];
  end
  
  if ischar(dateN) 
    dateN = datenum(dateN);
  end
  dateString = [datestr(dateN,'yyyy-mm-dd')];
  fprintf(['Loading SCTA data for ' dateString '\n']);
  fileA = [baseDir dateString '/SCTA-Accel2_' dateString '.nc'];
  fileAbackup = [baseDir dateString '/SCTA-Accel2_' dateString ' (1).nc'];
else
  fileA = baseDir;
end
  

data = [];

%% Load Acceleration Data
fprintf('  Loading Acceleration Data\n') 
if ~exist(fileA,'file')
  temp = fileA;
  fileA = fileAbackup;
  if ~exist(fileA,'file')
    disp(['get_sctaDay - No acceleration file found! - ' fileA]);
    data.t = [];
    data.a = [];
    data.as = [];
    data.T = [];
    fprintf('\n')
    return
  end
end
t = ncread(fileA,'time');
t = t/1000;
x = [];
ax = ncread(fileA,'accelX');
ay = ncread(fileA,'accelY');
az = ncread(fileA,'accelZ');
rss = ncread(fileA,'rss');
x = [ax ay az rss];
Ta = ncread(fileA,'temp');
x = [x Ta];
fprintf(['    ' int2str(length(t)) ' samples\n']);
data.t = datenum([repmat([1970 1 1 0 0],length(t),1) t]);
data.a = x(:,1:3);
data.as = x(:,4);
data.T = x(:,5);

%% Sort if times not in order
if any(diff(data.t)<0)
  disp('Accelerometer data not time sorted - Corrected')
  [data.t,i] = sort(data.t);
  data.a = data.a(i,:);
  data.as = data.as(i);
  data.T = data.T(i);
end

fprintf('\n')


      
      
      
      








