clear;
close all;
clc;

funHandles = main_functions_2d();  % Get the function handles structure
BSplineTransformationH= funHandles.BSplineTransformation;
ImageRegistrationH= funHandles.ImageRegistration;
k = [100,100];

% Load 3D volumes
data = load('twoConsecutiveVolumes.mat'); % Load 3D volume data

% Check what variables are inside the .mat file
vars = fieldnames(data);  
disp('Loaded Variables:');
disp(vars);

% Assuming the image is stored as a variable (e.g., 'elep')
original = data.(vars{1}); % Load the first variable dynamically
original = original( :,:,100);  % Convert to double if needed




% Assuming the image is stored as a variable (e.g., 'elep')
rotated = data.(vars{2}); % Load the first variable dynamically
rotated = rotated( :,:,100);  % Convert to double if needed

%rotated = [rotated,zeros(x,100)];
%shift
%original = imtranslate(original,[100, 0]);
%rotated = imtranslate(rotated,[100, 0]);
[x,y] = size(original);

[beta,regi,Dssd] = ImageRegistrationH(original,rotated,k);

z = [0;0];
z(1) = ceil(x/k(1))+3;
z(2) = ceil(y/k(2))+3;
U = zeros(x,y);
V = zeros(x,y);
for i = 1:x
    for j = 1:y
        new_u = BSplineTransformationH([i,j],beta,k,z);
        U(i,j) = new_u(1);
        V(i,j) = new_u(2);
    end
end
figure
streamslice(1:y,1:x,V,U)
%Plotting
hold all
figure('name', 'Images', 'NumberTitle', 'Off');
%circle
subplot(1,2,1);
imshow(original);
title('frame 1');

%elep Image
subplot(1,2,2);
imshow(rotated);
title('frame 2');
figure
%After Image registration
imshow(full(regi));
title('registered image');
%Differenz vorher
figure
subplot(1,2,1)
imshow(imabsdiff(original,rotated));
title('frame 1 - frame 2');
%Differenz nachher
subplot(1,2,2)
imshow(imabsdiff(original,full(regi)));
title('frame 1 - registered image');