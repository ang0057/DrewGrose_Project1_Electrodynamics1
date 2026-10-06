clear;
close all;
clc;

%% ========================================================================
% General space setup and variable declarations

% --- make the simulation space ---
% The goal: make a 0.5m x 0.02m 2D cross section of the full space and 
% place the field mill in the center. In this cross-sectional view, the 
% mill will have a 0.2m x 2.4mm footprint. 
dx = 0.001; % dx = 1mm
x = 0:dx:0.5; % 0.5m
dz = 0.0001; % dz = 0.1mm (finer since E varies in z)
z = 0:dz:0.02; % 0.02m or 20mm
%0.3mm = 0.0003m

% meshgrid
[Z, X] = meshgrid(z, x);
% ndgrid
% [X, Z] = ndgrid(x, z);
% ---------------------------------


% constitutive parameters
eps_0 = 8.854e-12; % same for copper since its eps_r = 1
sigma_0 = 0;
sigma_cu = 5.8e7;

% initialize V
V = zeros(size(X));

% implement boundary conditions to create the 100 V/m electric field
V(:, end) = 100;    % upper boundary in z
V(:, 1) = 0;        % lower boundary in z


% --- build field mill plates ---
% all plates will be 0.5mm thick (5 grid cells in z)
% d_shutter will be twice d shown on figure 1 of the instructions pdf
% d is assumed to be 0.3mm (3 grid cells in z)
% thus, d_shutter = 0.6mm (6 grid cells in z), and the total plate assembly
% is 24 grid cells tall in z

% grounded shutter (initial position)
V(150:250,120:124) = 10;

% sense plates 
% (these are 100V right now to visualize where they are in the space)
V(150:245,108:112) = 100;
V(255:350,108:112) = 100;

% grounded ground plate
V(150:350,100:104) = 10;

% -------------------------------


% --- simulate shutter motion ---
% TODO
% -------------------------------

%% ========================================================================
% MAIN SIMULATION LOOP

%scale variable to make it more simple
scale = ((dx*dz)^2)/(2*((dx^2)+(dz^2)));

tic
%iterate (solve) phi equation
for ii = 1:1000
    V(2:end-1, 2:end-1) = scale*(((V(3:end, 2:end-1) + ...
        V(1:end-2, 2:end-1))/(dx^2)) + ...
        ((V(2:end-1, 3:end) + V(2:end-1, 1:end-2))/(dz^2)));

    % reinforce grounded plate conditions
    V(150:250,120:124) = 10; % shutter
    V(150:350,100:104) = 10; % big ground plate
end
toc

%% ========================================================================
% Plotting

% --- auburn colormap generation ---
% function from TRACE repo:

% RGB values from [0:255]
orange = [232,97,0];
blue   = [11,35,65];
white  = [255,255,255];

%keep as an even number
length_colormap = 100;

% make blue to white
auburn_color_map = zeros(length_colormap,3);
auburn_color_map(1:length_colormap/2,1) = (interp1([1,length_colormap/2], [blue(1) white(1)], 1:length_colormap/2))';
auburn_color_map(1:length_colormap/2,2) = (interp1([1,length_colormap/2], [blue(2) white(2)], 1:length_colormap/2))';
auburn_color_map(1:length_colormap/2,3) = (interp1([1,length_colormap/2], [blue(3) white(3)], 1:length_colormap/2))';

% make white to orange
auburn_color_map(length_colormap/2+1:length_colormap,1) = (interp1([length_colormap/2+1,length_colormap], [white(1) orange(1)], length_colormap/2+1:length_colormap))';
auburn_color_map(length_colormap/2+1:length_colormap,2) = (interp1([length_colormap/2+1,length_colormap], [white(2) orange(2)], length_colormap/2+1:length_colormap))';
auburn_color_map(length_colormap/2+1:length_colormap,3) = (interp1([length_colormap/2+1,length_colormap], [white(3) orange(3)], length_colormap/2+1:length_colormap))';

% RGB values normalized from [0:1]
auburn_color_map = auburn_color_map/255;

% ----------------------------------

% Figure 1
figure;
imagesc(x, z, V'); 
colormap(auburn_color_map); colorbar;
axis xy;
title('Voltage as a Function of the 2-D Space');
xlabel('Length (m)'); ylabel('Height (m)');