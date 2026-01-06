clc; clear; close all;

% Define grid for 3D space
[x, y, z] = meshgrid(linspace(-2, 2, 20), linspace(-2, 2, 20), linspace(-2, 2, 20));

% Create Sphere (r = 1)
r = 1;
r2=0.9;
sphere_vol = (x.^2 + y.^2 + z.^2) <= r^2; 
cube_vol = (abs(x) <= r) & (abs(y) <= r) & (abs(z) <= r);
cube_vol2 = (abs(x) <= r2) & (abs(y) <= r2) & (abs(z) <= r2);
% Create Ellipsoid (Different radii)
a = 1.5;  % Radius along x-axis
b = 1;    % Radius along y-axis
c = 0.7;  % Radius along z-axis
ellipsoid_vol = (x.^2 / a^2 + y.^2 / b^2 + z.^2 / c^2) <= 1; 

% Plot Sphere
figure;
isosurface(x, y, z, sphere_vol, 0.5);
xlabel('X'); ylabel('Y'); zlabel('Z');
title('3D Sphere');
axis equal; grid on;
camlight; lighting phong;

% Plot Ellipsoid
figure;
isosurface(x, y, z, ellipsoid_vol, 0.5);
xlabel('X'); ylabel('Y'); zlabel('Z');
title('3D Ellipsoid');
axis equal; grid on;
camlight; lighting phong;

% Plot cube1
figure;
isosurface(x, y, z, cube_vol, 0.5);
xlabel('X'); ylabel('Y'); zlabel('Z');
title('3D Cube');
axis equal; grid on;
camlight; lighting phong;

% Plot cube2
figure;
isosurface(x, y, z, cube_vol2, 0.5);
xlabel('X'); ylabel('Y'); zlabel('Z');
title('3D Cube');
axis equal; grid on;
camlight; lighting phong;

% Save volumetric data
save('sphere_data.mat', 'sphere_vol');
save('ellipsoid_data.mat', 'ellipsoid_vol');
save('cube_data.mat', 'cube_vol');
save('cube_data2.mat', 'cube_vol2');