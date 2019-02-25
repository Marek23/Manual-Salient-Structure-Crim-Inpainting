close all; clear; clc;

mex mainSS.c

images = dir('C:\MAREK\MAGISTERKA\Obrazy\imgmaskSS\*.png');

for image=1:length(images)
    
images(image).name

test = strsplit(images(image).name,'_');

salNumb = str2double(test(1));
p_r     = str2double(test(2));

% bo potrzebuje wymiaru
I = im2double(imread(['C:\MAREK\MAGISTERKA\Obrazy\imgmask\' test{3}])); %images(image).name

s_r = ceil(sqrt(size(I,1)*0.02*size(I,2)*0.02));
while s_r < 10000 %% dla s_r > 9000 w algorytmie przyjmuj? ca?y obraz

clearvars -except image images alfa p_r s_r test salNumb
%%parametry
alfa = 0.2;

BrokenAreaColor = 0.95;

I   = im2double(imread(['C:\MAREK\MAGISTERKA\Obrazy\imgmask\'  test{3}])); %images(image).name
%I    = im2double(imread('triangle.png')); %images(image).name
C   = im2double(imread(['C:\MAREK\MAGISTERKA\Obrazy\conf\'     test{3}])); %images(image).name
SI  = double(imread(['C:\MAREK\MAGISTERKA\Obrazy\imgmaskSS\'   images(image).name])); %images(image).name
%C    = ones([size(I,1) size(I,2)]);

mask = double(1-((I(:,:,1) < 0.03) & ...
                   (    I(:,:,2) > BrokenAreaColor) & ...
                   (    I(:,:,3) < 0.03)));

[nx,ny,nz] = size(I);

tic
Ir = mainSS(nx,ny,nz,I(:),SI(:),mask(:),C(:),p_r,s_r,alfa,6);
t = toc;

I = reshape(Ir,[nx,ny,nz]);

imwrite(I, ['C:\MAREK\MAGISTERKA\Obrazy\salstruct\' images(image).name 'pr_' num2str(p_r) 'sr_' num2str(s_r) 'alfa_' num2str(alfa) 't_' num2str(t) '.png']);

%figure
%imshow(I)
s_r = s_r + 8000;
end

end
