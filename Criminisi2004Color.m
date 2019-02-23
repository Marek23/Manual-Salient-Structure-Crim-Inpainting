close all; clear; clc;

mex main.c

images = dir('C:\MAREK\MAGISTERKA\Obrazy\imgmask\*.png');

for image=1:length(images)
    
images(image).name

p_r=3;
while p_r < 24

% bo potrzebuje wymiaru
I = im2double(imread(['C:\MAREK\MAGISTERKA\Obrazy\imgmask\' images(image).name])); %images(image).name

s_r = ceil(sqrt(size(I,1)*0.02*size(I,2)*0.02))
while s_r < 10000 %% dla s_r > 9000 w algorytmie przyjmuj? ca?y obraz

clearvars -except image images alfa p_r s_r
%%parametry
alfa = 0.2;

BrokenAreaColor = 0.95;

I   = im2double(imread(['C:\MAREK\MAGISTERKA\Obrazy\imgmask\' images(image).name])); %images(image).name
%I    = im2double(imread('triangle.png')); %images(image).name
C   = im2double(imread(['C:\MAREK\MAGISTERKA\Obrazy\conf\'    images(image).name])); %images(image).name
%C    = ones([size(I,1) size(I,2)]);

mask = double(1-((I(:,:,1) < 0.03) & ...
                   (    I(:,:,2) > BrokenAreaColor) & ...
                   (    I(:,:,3) < 0.03)));

[nx,ny,nz] = size(I);

tic
Ir = main(nx,ny,nz,I(:),mask(:),C(:),p_r,s_r,alfa);
t = toc;

I = reshape(Ir,[nx,ny,nz]);

imwrite(I, ['C:\MAREK\MAGISTERKA\Obrazy\crimtest\' images(image).name 'pr_' num2str(p_r) 'sr_' num2str(s_r) 'alfa_' num2str(alfa) 't_' num2str(t) '.png']);

%figure
%imshow(I)
s_r = s_r + 8000;
end

p_r = p_r+4;
end

end
