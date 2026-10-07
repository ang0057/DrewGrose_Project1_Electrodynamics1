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

% Transimpedance amplifier feedback resistor (assuming 100 kOhm)
R_f = 1e5;

% shutter plate constant velocity: v = dx / dt
% assuming a physical shutter speed of 1 m/s (dt = dx / v = 1ms)
v_shutter = 1.0; 
dt = dx / v_shutter;

% implement boundary conditions to create the 100 V/m electric field
% convention: V(x, z); will be transposed when plotted

% since the space has a height of 0.02m, 2V at the top of it and 0V at the
% bottom results in the 100 V/m electric field
V = repmat(linspace(0, 2, size(Z,2)), size(Z,1), 1);


% --- build field mill plates ---
% -all plates will be 0.5mm thick (5 grid cells in z)
% -d_shutter, labeled d on figure 1 of the instructions pdf, will be twice
%  the distance between the sense plates and the big ground plate (d_lower)
% -d_lower is assumed to be 0.3mm (3 grid cells in z)
%   -thus, d_shutter = 0.6mm (6 grid cells in z), and the total plate 
%    assembly is 24 grid cells tall in z

% define plate locations/dimensions as index ranges
z_shutter = 120:124;
z_sense   = 109:113;
z_ground  = 101:105;

x_sense1  = 150:245; % sense plate 1 (left)
x_sense2  = 255:350; % sense plate 2 (right)
x_ground  = 150:350; % bottom ground plate

% now, enforce the grounded potential of each plate

% grounded shutter (initial position)
V(150:250, z_shutter) = 0;

% for reference: d_shutter is at z indices 114:119

% sense plates (these still need to be held at zero)
V(x_sense1, z_sense)  = 0;
V(x_sense2, z_sense)  = 0;

% for reference: d_lower is at z indices 106:108

% grounded ground plate
V(x_ground, z_ground) = 0;
% -------------------------------


% plot the space with initial conditions (example purposes)
fig1 = figure;
imagesc(x, z, V');
colormap(auburn_color_map);
c1 = colorbar; c1.Label.String = '(V/m)'; c1.Label.FontSize = 14;
axis xy;
title(sprintf('Voltage Distribution of the Space at t=0s'));
xlabel('Length (m)'); ylabel('Height (m)');

%% ========================================================================
% MAIN SIMULATION LOOP

% scale factor term of full update equation
scale_nonuniform = ((dx*dz)^2)/(2*((dx^2)+(dz^2)));
% scale = 1/(2*(dz/dx) + 2*(dx/dz));


% --- pre-relax at the initial position to fix artifacts ---
% (this worked somewhat. the initial spike for current and voltage changed
% from -6 to ~-3, so that will be good enough)
for kk = 1:1000
    V(2:end-1, 2:end-1) = scale_nonuniform * ( ...
        (V(3:end, 2:end-1) + V(1:end-2, 2:end-1)) / (dx^2) + ...
        (V(2:end-1, 3:end) + V(2:end-1, 1:end-2)) / (dz^2) );

    % reinforce grounded plate conditions and space boundary conditions
    V(150:250, z_shutter) = 0;
    V(x_sense1, z_sense)  = 0;
    V(x_sense2, z_sense)  = 0;
    V(x_ground, z_ground) = 0;
    V(:, 1)               = 0;   % bottom space boundary
    V(:, end)             = 2; % top space boundary
end
% ---------------------------------------------------------


% --- simulate shutter motion ---
% motion setup
num_timesteps = 101;
Q_sense1 = zeros(1, num_timesteps);
Q_sense2 = zeros(1, num_timesteps);

% set up Figure 2 to animate the motion during runtime
fig2 = figure;
colormap(auburn_color_map); 

% tic % the whole loop takes ~30s to run
for ii = 1:num_timesteps
    % clear previous shutter position
    if ii > 1
        prev_x_shutter = (150 + ii - 2) : (250 + ii - 2);
        % set V back to initial value at those grid cells
        V(prev_x_shutter, z_shutter) = repmat(z(z_shutter) * 100, ...
            length(prev_x_shutter), 1);
    end

    % scoot shutter to the right in x by 1 each iteration
    % (equates to adding ii-1 to the indices)
    % the -1 allows for the initial position to be used first
    curr_x_shutter = (150 + ii - 1) : (250 + ii - 1);
    V(curr_x_shutter, z_shutter) = 0;

    % covering old indices with a different voltage for testing motion
    V(150+(ii-2),120:124) = 1.2; 

    % solve Laplace update equation (non-uniform scale)
    for jj = 1:1000
        V(2:end-1, 2:end-1) = scale_nonuniform * ( ...
            (V(3:end, 2:end-1) + V(1:end-2, 2:end-1)) / (dx^2) + ...
            (V(2:end-1, 3:end) + V(2:end-1, 1:end-2)) / (dz^2) );

        % reinforce grounded plate conditions and space boundary conditions
        V(curr_x_shutter, z_shutter) = 0;
        V(x_sense1, z_sense)         = 0;
        V(x_sense2, z_sense)         = 0;
        V(x_ground, z_ground)        = 0;
        V(:, 1)                      = 0;   % bottom space boundary
        V(:, end)                    = 2; % top space boundary
    end

    % --- calculate charge (Q) on sense plates ---
    % E_z at top surface of sense plate
    z_above = 114; % and z = 113 is the top of the sense plates
    E_z_sense1 = (V(x_sense1, z_above) - V(x_sense1, 113)) / dz;
    E_z_sense2 = (V(x_sense2, z_above) - V(x_sense2, 113)) / dz;

    % surface charge density rho_s = eps_0 * E_z
    rho_s1 = eps_0 * E_z_sense1;
    rho_s2 = eps_0 * E_z_sense2;

    % total induced charge Q (assuming square plates)
    L_y = length(x_sense1) * dx; % ?
    Q_sense1(ii) = sum(rho_s1) * dx * L_y;
    Q_sense2(ii) = sum(rho_s2) * dx * L_y;
    % --------------------------------------------
    

    % update animation plot
    imagesc(x, z, V');
    colormap(auburn_color_map);
    c2 = colorbar; c2.Label.String = 'Electric Field Intensity (V/m)';
    c2.Label.FontSize = 14;
    axis xy;
    title(sprintf('Voltage Distribution (Step %d / %d)', ...
                                             ii, num_timesteps));
    xlabel('Length (m)'); ylabel('Height (m)');
    drawnow;
end
% toc
% -------------------------------

%% ========================================================================
% Transimpedance Amplifier Output

% calculate current: i(t) = dQ / dt
i_sense1 = [0, diff(Q_sense1) / dt];
i_sense2 = [0, diff(Q_sense2) / dt];

% output for differential transimpedance configuration
V_out = -R_f * (i_sense1 - i_sense2);

% build time vector to use in plots
time = (0:num_timesteps-1) * dt;

%% ========================================================================
% Comparison of Analytical and Numerical Solutions

% physical constants and parameters
E_0 = 100;
W = length(x_sense1) * dx;
v = v_shutter;

% analytical V_out calculation: V_out = -2 * R_f * eps_0 * E_0 * W * v
V_out_analytical = -2 * R_f * eps_0 * E_0 * W * v;

% extract numerical output result at steady-state (timestep index 50)
V_out_numerical = V_out(50);

% calculate percent difference relative to analytical baseline
percent_diff = abs(V_out_numerical - ...
                        V_out_analytical) / abs(V_out_analytical) * 100;

% display comparison results in command window
fprintf('\n=======================================================\n');
fprintf('   ELECTRIC FIELD MILL SOLUTION COMPARISON RESULTS\n');
fprintf('=======================================================\n');
fprintf('Analytical Output Voltage  : %10.4f uV\n', V_out_analytical * 1e6);
fprintf('Numerical Output Voltage   : %10.4f uV (at steady-state)\n', V_out_numerical * 1e6);
fprintf('Absolute Difference        : %10.4f uV\n', abs(V_out_numerical - V_out_analytical) * 1e6);
fprintf('Percent Difference         : %10.2f %%\n', percent_diff);
fprintf('=======================================================\n\n');

%% ========================================================================
% Plotting

fig3 = figure; % Figure 3
subplot(3,1,1);
plot(time, Q_sense1 * 1e9, 'b', time, Q_sense2 * 1e9, 'r', 'LineWidth', 1.5);
title('Induced Charge on Sense Plates');
xlabel('Time (s)'); ylabel('Charge (nC)');
legend('Sense Plate 1', 'Sense Plate 2'); grid on;

subplot(3,1,2);
plot(time, i_sense1 * 1e6, 'b', time, i_sense2 * 1e6, 'r', 'LineWidth', 1.5);
title('Current Flowing to Amplifier Inputs');
xlabel('Time (s)'); ylabel('Current (\mu A)');
legend('i_1(t)', 'i_2(t)'); grid on;

subplot(3,1,3);
plot(time, V_out, 'k', 'LineWidth', 1.5);
title('Transimpedance Amplifier Output Voltage V_{out}(t)');
xlabel('Time (s)'); ylabel('Voltage (V)');
grid on;

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