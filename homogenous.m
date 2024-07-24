clc; clear; close all;

% This script file is solving wave equation on square plate using FE
% Check https://www.mathworks.com/help/pde/ug/wave-equation.html for
% details

tic; 

%% define the parameters for the wave equation

% c1 = c2 = (14.3e+3)/3 

each = 1; % can add loop if running for a range of values 

data = load("velocity3D.mat");
tlist = load("t_array.mat");

% make sure these parameters match the actual variable name 
uIntrp = data.velocity3D;
yq = data.uniqueY;
xq = data.uniqueX;
tlist=tlist.time_steps;

yq = yq * 10^-3; % changing units to m from mm
xq = xq * 10^-3;

C1_x = 15e-3; % where the center of inclusion could be

%% video 
figure;
umax = max(uIntrp(:));
umin = min(uIntrp(:));

videoFileNameBase = "DispAnimation3D_200_homogenous_";
videoFileName = strcat(videoFileNameBase, num2str(each), ".avi");

if exist(videoFileName, 'file')
    delete(videoFileName);
end

v = VideoWriter(videoFileName);
open(v);
disp(videoFileName);

radius = 5e-3;
width = 60.e-3;
totalLength = 60.e-3;

for i = 1:size(uIntrp, 3) % iterate over time
    surf(xq, yq, uIntrp(:, :, i));
    shading interp;
    axis([-totalLength totalLength -width width umin umax]);

    clim([umin umax]*0.75); 
    colormap('jet');
    colorbar;
    view(2); 

    % displays inclusion circle on video:

    % hold on;
    % x_center = 0.65;
    % y_center = 0.5;
    % radius_norm = radius / (2 * totalLength);
    % annotation('ellipse', [x_center - radius_norm, ...
    %     y_center - radius_norm, 2*radius_norm-0.012, 2*radius_norm], ...
    %     'Color', 'r', 'LineWidth', 0.8);
    % hold off;

    xlabel('x');
    ylabel('y');
    zlabel('u');

    frame = getframe(gcf); 
    writeVideo(v, frame); 
end

close(v);

%% calculate shear wave speed
getGroupVel(uIntrp, yq', xq', tlist, each, C1_x);
