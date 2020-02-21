close all; clear; clc;
mex mainSS.c

imageName   = 'test.png';
salientName = 'test_s.png';
resultName  = 'output.png';
salNumb     = 3;

%parametr Criminisi
p_r  = 5;
alfa = 0.2;

s_r = 10000; %% dla s_r > 9000 w algorytmie przyjmuj? ca?y obraz

SI  = double(imread(salientName));
I   = im2double(imread(imageName));

[nx,ny,nz] = size(I);

mask = double(1-((I(:,:,1) == 0 ) & ...
                (  I(:,:,2) == 1) & ...
                (  I(:,:,3) == 0)));
C = mask;


tic
Ir = mainSS(nx,ny,nz,I(:),SI(:),mask(:),C(:),p_r,s_r,alfa,6);
t = toc;

I = reshape(Ir,[nx,ny,nz]);

imwrite(I, [resultName 'pr_' num2str(p_r) 'sr_' num2str(s_r) 'alfa_' num2str(alfa) 't_' num2str(t) '.png']);

figure
imshow(I)
