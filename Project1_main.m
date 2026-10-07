clear;
close all;
clc;

%% ========================================================================
% General space setup and variable declarations
auburn_color_map = Generate_Auburn_Colormap_v0();

% --- make the simulation space ---
% The goal: make a 0.5m x 0.02m 2D cross section of the full space and 
% place the field mill in the center. In this cross-sectional view, the 
% mill will have a 0.2m x 2.4mm footprint. 
dx = 0.001; % dx = 1mm
x = 0:dx:0.5; % 0.5m
dz = 0.0001; % dz = 0.1mm (finer since E varies in z)
z = 0:dz:0.02; % 0.02m or 20mm

% z << x allows us to ignore fringing fields

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
% V = zeros(size(X));

% implement boundary conditions to create the 100 V/m electric field
% convention: V(x, z); will be transposed when plotted

% since the space has a height of 0.02m, 2V at the top of it and 0V at the
% bottom results in the 100 V/m electric field
V = repmat(linspace(0, 2, size(Z,2)), size(Z,1), 1);

% old:
% V(:, end) = 100;    % upper boundary in z
% V(:, 1) = 0;        % lower boundary in z


% --- build field mill plates ---
% all plates will be 0.5mm thick (5 grid cells in z)
% d_shutter will be twice d shown on figure 1 of the instructions pdf
% d is assumed to be 0.3mm (3 grid cells in z)
% thus, d_shutter = 0.6mm (6 grid cells in z), and the total plate assembly
% is 24 grid cells tall in z

% grounded shutter (initial position)
V(150:250,120:124) = 0;
% V(150:250,114:119) = 0; % shadow (shielding the sense plate)

% sense plates 
% (these still need to be held at zero)
V(150:245,109:113) = 0;
V(255:350,109:113) = 0;

% d
% V(150:350,106:108) = 0;

% grounded ground plate
V(150:350,101:105) = 0;

% -------------------------------

%% ========================================================================
% MAIN SIMULATION LOOP

%scale variable term (handles denominator of full update equation)
scale_nonuniform = ((dx*dz)^2)/(2*((dx^2)+(dz^2)));
scale = 1/(2*(dz/dx) + 2*(dx/dz));

% --- simulate shutter motion ---

% Figure 1
figure;
imagesc(x, z, V'); 
colormap(auburn_color_map); 
colorbar;
axis xy;
title('Voltage as a Function of the 2-D Space');
xlabel('Length (m)'); ylabel('Height (m)');

tic
for ii = 1:101

    % scoot shutter to the right in x by 1 each iteration
    % (equates to adding ii-1 to the indices)
    % the -1 allows for the initial position to be used first
    V(150+(ii-1):250+(ii-1),120:124) = 0;

    % covering old indices with a different voltage for testing motion
    V(150+(ii-2),120:124) = 1.2; 

    % tic
    % solve update equation (non-uniform scale)
    for jj = 1:1000
        V(2:end-1, 2:end-1) = scale_nonuniform*(((V(3:end, 2:end-1) + ...
            V(1:end-2, 2:end-1))/(dx^2)) + ...
            ((V(2:end-1, 3:end) + V(2:end-1, 1:end-2))/(dz^2)));

        % reinforce grounded plate conditions
        % V(150:250,120:124) = 0; % shutter
        V(150+(ii-1):250+(ii-1),120:124) = 0;
        V(150:350,100:104) = 0; % big ground plate
        % sense plates
        V(150:245,109:113) = 0;
        V(255:350,109:113) = 0;
    end
    % toc

    % tic
    % solve update equation (uniform scale)
    % for ii = 1:1000
    %     % V(2:end-1, 2:end-1) = scale*(((V(3:end, 2:end-1) + ...
    %     %     V(1:end-2, 2:end-1))/(dx^2)) + ...
    %     %     ((V(2:end-1, 3:end) + V(2:end-1, 1:end-2))/(dz^2)));
    %
    %     V(2:end-1, 2:end-1) = scale*( ...
    %         (dz/dx)*(V(1:end-2, 2:end-1) - V(3:end, 2:end-1)) + ...
    %         (dx/dz)*(V(2:end-1, 1:end-2) - V(2:end-1, 3:end)) );
    %
    %     % reinforce grounded plate conditions
    %     V(150:250,120:124) = 0; % shutter
    %     V(150:350,100:104) = 0; % big ground plate
    %     % sense plates
    %     V(150:245,109:113) = 0;
    %     V(255:350,109:113) = 0;
    % end
    % toc

    imagesc(x, z, V');
    drawnow;
end
toc
% -------------------------------

%% ========================================================================
% Plotting
% auburn_color_map = Generate_Auburn_Colormap_v0();

% Figure 2
figure;
imagesc(x, z, V'); 
colormap(auburn_color_map); 
colorbar;
axis xy;
title('Voltage as a Function of the 2-D Space');
xlabel('Length (m)'); ylabel('Height (m)');

%% ========================================================================
% Functions

% --- auburn colormap generation ---
% function from TRACE repo:

function auburn_color_map = Generate_Auburn_Colormap_v0()
% This function generates the Auburn colors colormap where
% blue is low value and orange is high value. RGB values
% according to https://ocm.auburn.edu/brand-center/_assets/pdf/au-colorpalette-primary-supporting.pdf

% INPUTS: None.

% OUTPUTS: 100 by 3 array of Auburn blue to orange

% DEPENDENCIES: none

% Author: Clint Snider
% last edited: 02/05/2025

% DESIRED UPDATES: Maybe make length 256.

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
end

% ----------------------------------