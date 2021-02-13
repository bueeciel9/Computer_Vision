% Solution - problem 1, Homework 2
%
% Chroma key compositing algorithm.
% Replaces a green screen background in an image with a new background.
%
function image_out = chroma(fname_main, fname_bkg, thresh1, thresh2)
%
% Inputs:
%fname_main = character string containing filename for main image
%fname_bkg = character string containing filename for background image
% Output:
%image_out = output array containing composite output image
% Initialize variables
nsamples = 3; %number of points for sampling background color
rgb_samp = zeros(nsamples,3);
%sampled RGB colors
A = [[.299,.587,.114];
%define color transformation matrix...
[-.169,-.331,.5];
[.5,-.419,-.081]];
b = [0,128,128]';
%...and offset
% Read main image and background image
image_main = imread(fname_main);
image_bkg = imread(fname_bkg);
%Find dimensions of main image
n1 = size(image_main,1);
n2 = size(image_main,2);
N = n1 * n2;


%Resize background image to same as main
image_bkg = imresize(image_bkg,[n1,n2]);
figure;
imagesc(image_bkg);
axis image;
title('Background image');
% Display the main image
figure;
imagesc(image_main);
%display image
axis image;
%make image shape correct
title('Main image');
hold on;


%allows future plotting steps to overlay on image
% Sample color of green screen background at 3 points, displaying dot at
% each point, as you click
disp('Use mouse to select background points');
for i = 1:nsamples
[x,y] = ginput(1); %get x,y coordinates at a point
x = round(x);
%round the sampled coordinates
y = round(y);
%round the sampled coordinates
rgb_samp(i,:) = image_main(y,x,:); %Note that 'ginput' uses (x,y) notation
%whereas matlab arrays are indexed
%like matrices
plot(x,y,'.');
%draw dot at (x,y) to show your selection
end


%Find average RGB color of sampled background points
bckval_RGB = mean(rgb_samp,1)';
%Convert RGB color to YCbCr color
bckval_YCbCr = A * bckval_RGB + b;

%Convert main image to YCbCr
%NOTE: The following steps can be done one pixel at a time using loops,
%but it's generally better to use matlab's matrix calculation capabilities
lex_RGB = double(reshape(image_main,N,3))'; %rearrange pixels as list
bb = repmat(b,1,N);
%replicate b into a matrix
lex_YCbCr = A * lex_RGB + bb;
%convert whole image in one command
image_main_YCbCr = reshape(lex_YCbCr',n1,n2,3); %put it back to original array
format
%Calculate every pixel's distance from the average background (in (Cb,Cr)
%space)
dist = sqrt((image_main_YCbCr(:,:,2) - bckval_YCbCr(2)).^2. +(image_main_YCbCr(:,:,3) - bckval_YCbCr(3)).^2.);
%Display distance image
figure;
%open a new figure
colormap(gray);
%set color map to grayscale
axis image;
%make axes proportionate to image
imagesc(dist);
%display distance map
title('Distance image');
%Calculate alpha mask
%(I got lazy here and used loops, which you should try to avoid)
alpha = zeros(n1,n2);
for i=1:n1
for j=1:n2
    if dist(i,j) > thresh2 
        alpha(i,j) = 1.;
elseif dist(i,j) > thresh1
    alpha(i,j) = (dist(i,j) - thresh1) / (thresh2-thresh1);
    end
end
end
%Display the alpha mask
alpha = double(alpha);
figure;
colormap(gray);
axis image;
imagesc(alpha);
title('alpha mask');
%Make composite image
image_out = zeros(n1,n2,3);
for i=1:3
image_out(:,:,i) = alpha .* double(image_main(:,:,i)) + (1 - alpha).* double(image_bkg(:,:,i));
end
image_out=uint8(round(image_out));




