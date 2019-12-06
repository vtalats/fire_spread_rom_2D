function [ ] = gen_movie_2D( U, xy, N, frate, fname)
%GEN_MOVIE Create a movie of the 2D solution
%   This creates a movie fname.avi of the solution over time given a frame
%   rate of frate


% % Clear up space 
% close all;

% Get the total number of frames
num_frames = size(U,1);

% Set the axis size the same on every loop to prevent auto re-sizing
z_buff = 0.1*max(max(U));
z_axis_lim = [0.9.*min(min(U)) 1.1*max(max(U))];

xplt = reshape(xy(:,1),N(1),N(2));
yplt = reshape(xy(:,2),N(1),N(2));

wildfire = mycolormap(2);

% Plot the initial conditions
fid = figure;
surf(xplt,yplt,reshape(U(1,:)',N(1),N(2)),'EdgeColor','none')
colormap(wildfire)
zlim(z_axis_lim);
title('Advection-Diffusion-Reaction - 2D');

 
% Set up the movie.
writerObj = VideoWriter([fname '.avi']); % Name it.
writerObj.FrameRate = frate; % How many frames per second.
open(writerObj); 
 
for k=1:num_frames 
  % Makes sure you use your desired frame.
  figure(fid); 

  % Plot each frame
  surf(xplt,yplt,reshape(U(k,:)',N(1),N(2)),'EdgeColor','none')
  colormap(wildfire)
  zlim(z_axis_lim);
  title('Heat Eqn - 2D');

  % 'gcf' can handle if you zoom in to take a movie.
  frame = getframe(gcf); 
  writeVideo(writerObj, frame);
end

close(writerObj); % Saves the movie.


end

