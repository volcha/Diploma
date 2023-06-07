close all; clear all;

orig_image = imread('88.jpg'); % считывание исходного изображения
orig_image = rgb2gray(orig_image); % перевод изображения в полутон. вид
orig_image = im2double(orig_image);
figure(1); imshow(orig_image,[]); % вывод исходного изображения

[n_size,m_size] = size(orig_image); % размеры изображения
number_of_pixels = n_size * m_size; % количество пикселей

% усиление размытости
filtr = fspecial('gaussian',[7 7],4);
blur_image = imfilter(orig_image,filtr,'same');
figure(2); imshow(blur_image);
%blur_image = orig_image;

wnm = 'db10'; % выбор типа вейвлета
N = 5; %уровень разложения
[c,s] = wavedec2(blur_image,N,wnm); % декомпозиция
s
[mc,nc]=size(c);
nzs = 0;
no_blur = blur_image;
for i = 1:1:N
    Cn=c;
    if (i == 1) 
        nzs = nzs + s(i,1)*s(i,2);
    else
        nzs = nzs + s(i,1)*s(i,2)*3;
    end
    Cn(1:nzs) = zeros(1,nzs);
    X = waverec2(Cn,s,wnm); % реконструкция
    X = imadjust(X);
    no_blur = no_blur + blur_image + X;
end
no_blur = no_blur -  min(no_blur(:)); %нормировка
no_blur = no_blur/max(no_blur(:));
no_blur = imadjust(no_blur); 
figure(44); imshow(no_blur,[]);
subtraction = orig_image - no_blur;
pow = subtraction.^2;
mse = sum(pow(:)) / number_of_pixels

W = [ 1 1 1; 1 -8 1; 1 1 1];% фильтр Лапласиана
G = imfilter(blur_image,W,'replicate'); % фильтрация
GL = blur_image-G; % увеличение резкости
GL = imadjust(GL); 
figure(5); imshow(GL,[]);
subtraction = orig_image - GL;
pow = subtraction.^2;
mse_filtr = sum(pow(:)) / number_of_pixels