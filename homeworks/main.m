%Main program to solve Problem 1, Homework 2
fname_main='lady.jpg';
fname_bkg='amalfi.jpg';
thresh1 = 10;
thresh2 = 50; % thresh 차이에 따라 어떤 것이 달라지는지 확인해보기. 

%name of main image
%name of background image
%lower threshold
%upper threshold

%Call the compositing function (see chroma.m)
image_out = chroma(fname_main, fname_bkg, thresh1, thresh2);
%Display the result
figure;
imagesc(image_out);
axis image;
title('Composite image');
