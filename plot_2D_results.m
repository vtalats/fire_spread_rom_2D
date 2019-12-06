function [ ] = plot_2D_results( T, Tred, xy, N, ts, t )
%PLOT_2D_RESULTS Plot the 2 D results
%   Detailed explanation goes here


% % Clear up space 
% close all;

% Get the total number of frames
snaptime = round(t(ts));

% Set the axis size the same on every loop to prevent auto re-sizing
z_axis_lim = [0.9.*min(min(T)) 1.1*max(max(T))];

xplt = reshape(xy(:,1),N(1),N(2));
yplt = reshape(xy(:,2),N(1),N(2));

wildfire = mycolormap(2);

% Plot the FOM
figure('Name', 'FOM Fire-Spread 2D');
surf(xplt,yplt,reshape(T(ts,:)',N(1),N(2)),'EdgeColor','none')
colormap(wildfire)
zlim(z_axis_lim);
xlabel('x (m)');
ylabel('y (m)');
zlabel('Temp (K)');
title(sprintf('Fire-Spread 2D - FOM - %ds',snaptime));
saveas(gcf,sprintf('fs_2d_fom_%d.png',snaptime))

% Plot the POD/DEIM
figure('Name', 'POD/DEIM Fire-Spread 2D');
surf(xplt,yplt,reshape(Tred(ts,:)',N(1),N(2)),'EdgeColor','none')
colormap(wildfire)
zlim(z_axis_lim);
xlabel('x (m)');
ylabel('y (m)');
zlabel('Temp (K)');
title(sprintf('Fire-Spread 2D - POD/DEIM - %ds',snaptime));
saveas(gcf,sprintf('fs_2d_poddeim_%d.png',snaptime))

 
% % Set up the movie.
% writerObj = VideoWriter([fname '.avi']); % Name it.
% writerObj.FrameRate = frate; % How many frames per second.
% open(writerObj); 
%  
% for k=1:num_frames 
%   % Makes sure you use your desired frame.
%   figure(fid); 
% 
%   % Plot each frame
%   surf(xplt,yplt,reshape(T(k,:)',N(1),N(2)),'EdgeColor','none')
%   colormap(wildfire)
%   zlim(z_axis_lim);
%   title('Heat Eqn - 2D');
% 
%   % 'gcf' can handle if you zoom in to take a movie.
%   frame = getframe(gcf); 
%   writeVideo(writerObj, frame);
% end
% 
% close(writerObj); % Saves the movie.
% 
% 
% end


end

