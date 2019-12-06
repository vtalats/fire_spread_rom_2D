function [ cmap ] = mycolormap( color_idx )
%MYCOLORMAP Output a custom colormap
%   Define custom color maps
%
% Input
%   color_idx :
%     0 - gray scale
%     1 - reverse gray scale
%     2 - burning field
%     default = 0
% Ouput
%   cmap : custom colormap
%------------------------------------------------------------------
% Author: Alan Lattimer, JENSEN HUGHES 
% Date:   July 1, 2016
%------------------------------------------------------------------

  if nargin < 1
    color_idx = 0;
  end
  
  switch color_idx
    case 0
      g = linspace(0,1,1000)';
      cmap = [g,g,g];
    case 1
      g = linspace(1,0,1000)';
      cmap = [g,g,g];
    case 2
      r = [0.2*ones(1,5), linspace(1.0,1.0,1000)]';
      g = [0.4*ones(1,5), linspace(1.0,0.0,1000)]';
      b = [0.0*ones(1,5), linspace(0.6,0.0,1000)]';
      cmap = [r,g,b];
    otherwise
      g = linspace(0,1,1000)';
      cmap = [g,g,g];
  end
  
      
end

