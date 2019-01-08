function data = get_tiltDay(baseDir,baseName,dateN)
% Gets a days worth of SCTA data with options to only get specific channels
%
% Usage
%   data = get_tiltDay(date,baseName,dateN)
%
% Inputs
%   BASEDIR  - EITHER the base directory for where to fine the file so that
%              '/Users/wilcock/test/2017-08-02/SCTA-Accel_2017-08-24.nc'
%              would have a base directory of 
%              '/Users/wilcock/test'
%              OR the actual file name in which case BASENAME and DATEN is ignored
%   BASENAME - Base name for the tilt data file of the form 'SCTA-Tilt_20_19_' where the
%              full file name for August 9 2018 is SCTA-Tilt_20_19_2018-08-09.nc
%   DATEN    - date number or string for day to get ifs BASEDIR is the base
%              directory and not the file to load
%
% Outputs
%   DATA    - Stucture with data in column vectors/matrices of the following fields
%               t - Time for acceleration 1
%               a  - Acceleration in 3 columns of x, y, z
%               T - Temperture of the accelerometer 1

%% Preliminaries
if nargin~=3 && nargin~=1
  error('get_sctaDay requires 11-2 input arguements')
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
  fileA = [baseDir dateString '/' baseName dateString '.nc'];
else
  fileA = baseDir;
end
  

data = [];

%% Load Acceleration Data
fprintf('  Loading Acceleration Data\n') 
if ~exist(fileA,'file')
  disp(['get_sctaDay - No acceleration file found! - ' fileA]);
  data.t = [];
  data.a = [];
  data.as = [];
  data.T = [];
else
  t = ncread(fileA,'time');
  t = t/1000;
  x = [];
  ax = ncread(fileA,'accelX');
  ay = ncread(fileA,'accelY');
  x = [ax ay];
  Ta = ncread(fileA,'temp');
  x = [x Ta];
  fprintf(['    ' int2str(length(t)) ' samples\n']);
  data.t = datenum([repmat([1970 1 1 0 0],length(t),1) t]);
  data.a = x(:,1:2);
  data.T = x(:,3);
end

fprintf('\n')


      
      
      
      








