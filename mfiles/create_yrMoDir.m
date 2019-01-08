function [plotDir,plotSubDir] = create_yrMoDir(parentDir,d,varargin)
% Creates a Year/Month directory structure with a path to current day to store plots
%
% Usage 
%   plotSubDir = create_yrMoDir(parentDir,d,varargin)
%
% Input
%   parentDir - Parent directory with or without trailing '/'
%   d         - Date string for current date
%   varargin  - Names of one or more sub-subdirectories to create without '/'
%
% Output
%    plotDir    - Parent directory with trailing '/'
%    plotSubDir - Subdirectory for plot files ending in 'yyyy/mm/' format

if ~isempty(parentDir) && parentDir(end)=='/'
  plotDir = parentDir(1:end-1);
else
  plotDir = parentDir;
end
mmString = datestr(d,'mm');
if ~exist(plotDir,'dir')
  error('Plot directory does not exist')
else  
  yrString = datestr(d,'yyyy');
  plotSubDir = [plotDir '/' yrString];
  if ~exist(plotSubDir,'dir')
    mkdir(plotSubDir)
  end
  mnString = datestr(d,'mm');
  plotSubDir = [plotDir '/' yrString '/' mnString];
  if ~exist(plotSubDir,'dir')
    mkdir(plotSubDir)
  end
  plotSubDir = [plotSubDir '/'];
  for i = 1:length(varargin)
    plotSubSubDir = [plotSubDir varargin{i}];
    if ~exist(plotSubSubDir,'dir')
      mkdir(plotSubSubDir)
    end
  end
end
plotDir = [plotDir '/'];