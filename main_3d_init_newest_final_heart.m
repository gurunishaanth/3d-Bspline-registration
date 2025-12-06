clc; 
clear;
close all;

% Load the required functions
funHandles = main_functions_3D_b_4();  % Get the function handles structure
ImageRegistrationH = funHandles.ImageRegistration3D;

% Define B-spline grid resolution
k = [5, 5, 5]; % Control point spacing for 3D

% Load 3D volumes
data = load('usSlabs_Diastole_Systole_reducedSize.mat'); % Load 3D volume data

% Check what variables are inside the .mat file
vars = fieldnames(data);  
disp('Loaded Variables:');
disp(vars);

% Assuming the image is stored as a variable (e.g., 'elep')
original = data.(vars{1}); % Load the first variable dynamically
original = im2double(original);  % Convert to double if needed
R = (original - min(original(:))) / (max(original(:)) - min(original(:)));



% Assuming the image is stored as a variable (e.g., 'elep')
rotated = data.(vars{2}); % Load the first variable dynamically

rotated = im2double(rotated);  % Convert to double if needed
T = (rotated - min(rotated(:))) / (max(rotated(:)) - min(rotated(:)));


[x_dim, y_dim, z_dim] = size(R);


% --- Register ---
[beta, Registered, Dssd, U, V, W]  = ImageRegistrationH(R, T, k);

% --- Show Middle Slice & Difference ---
midSlice = round(size(R, 3) / 2);
figure;
subplot(1,3,1); imshow(R(:,:,midSlice), []); title('Reference');
subplot(1,3,2); imshow(Registered(:,:,midSlice), []); title('Registered');
subplot(1,3,3); imshow(Registered(:,:,midSlice) - R(:,:,midSlice), []); title('Difference');

% --- Deformation Field Visualization ---
[m, n, p] = size(R);

% Subsample grid for quiver
[xg, yg, zg] = ndgrid(1:m, 1:n, 1:p);

% Extract displacement field from new_u


% You actually already get U,V,W from ImageRegistration3D:
% [beta, Registered, Dssd, U, V, W] = ImageRegistrationH(R, T, k);

% For quiver3, you probably want a subsampled grid to avoid too many arrows
step = 1;  % every 5 voxels
figure;
quiver3(xg(1:step:end,1:step:end,1:step:end), ...
        yg(1:step:end,1:step:end,1:step:end), ...
        zg(1:step:end,1:step:end,1:step:end), ...
        U(1:step:end,1:step:end,1:step:end), ...
        V(1:step:end,1:step:end,1:step:end), ...
        W(1:step:end,1:step:end,1:step:end), 0);
hold on;

% Optional: mark active voxels
[rx, ry, rz] = ind2sub(size(R), find(R));
[gx, gy, gz] = ind2sub(size(T), find(T));
plot3(rx, ry, rz, 'rx');
plot3(gx, gy, gz, 'gx');

xlabel('X'); ylabel('Y'); zlabel('Z');
title('Deformation Field (B-spline)');
axis tight; grid on; view(3);


% Optional: mark active voxels
[rx, ry, rz] = ind2sub(size(R), find(R));
[gx, gy, gz] = ind2sub(size(T), find(T));
plot3(rx, ry, rz, 'rx');
plot3(gx, gy, gz, 'gx');

xlabel('X'); ylabel('Y'); zlabel('Z');
title('Deformation Field (B-spline)');
axis tight; grid on; view(3);

%%

% Pick a middle slice along z
midSlice = round(p/2);

% Extract slice data
Ux = U(:,:,midSlice);
Vy = V(:,:,midSlice);
Rslice = R(:,:,midSlice);
Tslice = T(:,:,midSlice);

[m, n] = size(Rslice);
[xg, yg] = meshgrid(1:n, 1:m);

% Subsample for clarity
step = 1;

figure; hold on;

% --- Plot all deformation vectors in light gray ---
quiver(xg(1:step:end,1:step:end), ...
       yg(1:step:end,1:step:end), ...
       Ux(1:step:end,1:step:end), ...
       Vy(1:step:end,1:step:end), 0, 'Color', [0.7 0.7 0.7]);

% --- Highlight reference voxel arrows in red ---
[rx, ry] = find(Rslice);
for i = 1:length(rx)
    quiver(ry(i), rx(i), Ux(rx(i), ry(i)), Vy(rx(i), ry(i)), 0, 'r', 'LineWidth', 1.5);
end

% --- Highlight target voxel arrows in green ---
[gx, gy] = find(Tslice);
for i = 1:length(gx)
    quiver(gy(i), gx(i), Ux(gx(i), gy(i)), Vy(gx(i), gy(i)), 0, 'g', 'LineWidth', 1.5);
end

axis equal tight;
set(gca,'YDir','reverse');
title(['2D Deformation Field + Active Voxels, z = ', num2str(midSlice)]);
xlabel('X'); ylabel('Y');
legend({'All Deformation','Reference Voxels','Target Voxels'});


[x, y, z] = meshgrid(linspace(-2, 2, y_dim), linspace(-2, 2, x_dim), linspace(-2, 2, z_dim));
% Visualization
figure('Name', '3D Image Registration', 'NumberTitle', 'Off');
subplot(1, 3, 1);
isosurface(x, y, z, R, 0.5);
title('sphere');
axis equal; grid on;
camlight; lighting phong;


subplot(1, 3, 2);
isosurface(x, y, z, T, 0.5);
title('cube');
axis equal; grid on;
camlight; lighting phong;


subplot(1, 3, 3);
isosurface(x, y, z, full(Registered), 0.5);
title('Registered Volume');
axis equal; grid on;
camlight; lighting phong;