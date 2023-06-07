close all; clear all;

%% �������� �����������
orig_image = imread('4.png'); % ���������� ��������� �����������
orig_image = rgb2gray(orig_image); % ������� ����������� � �������. ���
orig_image = im2double(orig_image);
figure(1); imshow(orig_image,[]); % ����� ��������� �����������

%% ��� ���� "���� & �����"
salt_pepper_noise = imnoise(orig_image,'salt & pepper',0.4); % ��������� ���� ���� "���� � �����"
figure(2); imshow(salt_pepper_noise,[]); % ����� ����������� �����������

%% ������� ��� ���� ��������� ����������
%Wavefuncs = {'haar', 'sym2', 'sym3', 'sym4', 'sym5', 'sym6', 'sym7', 'sym8', 'bior 1.1', 'bior 1.3', 'bior 1.5', 'bior 2.2', 'bior 2.4', 'bior 2.6', 'bior 2.8', 'bior 3.1', 'bior 3.3', 'bior 3.5', 'bior 3.7', 'bior 3.9', 'bior 4.4', 'bior 5.5', 'bior 6.8', 'db2', 'db3', 'db4', 'db5', 'db6', 'db7', 'db10', 'db20', 'coif2',  'coif3',  'coif4',  'coif5', 'rbio1.1', 'rbio 1.3', 'rbio 1.5', 'rbio 2.2', 'rbio 2.4', 'rbio 2.6', 'rbio 2.8', 'rbio 3.1', 'rbio 3.3', 'rbio 3.5', 'rbio 3.7', 'rbio 3.9', 'rbio 4.4', 'rbio 5.5', 'rbio 6.8', 'dmey'}; % ��� ��������
Wavefuncs = {'rbio2.2'}; % ��� ��������
[si1, si2] = size(Wavefuncs);
N = 1; % ������� ����������
[n_size,m_size] = size(orig_image); % ������� �����������
number_of_pixels = n_size * m_size; % ���������� ��������

min_err = 1;
best_wave = 'a';
for wave = 1:1.0:si2 % ������������� ������ ��������
    Wavefun = Wavefuncs{wave};
    [C,S] = wavedec2(salt_pepper_noise,N,Wavefun); % ������������
    [mc,nc]=size(C);
    Cn = C;
    S
    nk = 0;
    for j = 1:N
        nk = nk + S(N+2-j,1) * S(N+2-j,2) * 3;
        nf = nc - nk;
        Cn(nf:end) = zeros(1,nc-nf+1);
        no_salt_pepper = waverec2(Cn,S,Wavefun); % �������������
        no_salt_pepper = no_salt_pepper -  min(no_salt_pepper (:));
        no_salt_pepper  = no_salt_pepper /max(no_salt_pepper (:));
        no_salt_pepper = imadjust(no_salt_pepper); 
        figure(4); imshow(no_salt_pepper,[]);
        subtraction = orig_image - no_salt_pepper;
        pow = subtraction.^2;
        err = sum(pow(:)) / number_of_pixels
        %pause;
        if (min_err > err) % ���� ����������� ������
            min_err = err;
            best_wave = Wavefun;
        end
    end
    Wavefun
    %pause;
end
min_err
best_wave

%% ��������� � ������� ��������������� ������� �����������
FT_img = fft2(double(salt_pepper_noise));
n = 2; % ����� �� ���� �����
D0 = 200; % ��� �������� ����� ��������
u = 0:(n_size-1);
v = 0:(m_size-1);
idx = find(u > n_size/2);
u(idx) = u(idx) - n_size;
idy = find(v > m_size/2);
v(idy) = v(idy) - m_size;
[V, U] = meshgrid(v, u);
D = sqrt(U.^2 + V.^2);
H = 1./(1 + (D./D0).^(2*n));
G = H.*FT_img;
output_image = real(ifft2(double(G))); 
output_image = output_image -  min(output_image(:)); %����������
output_image = output_image/max(output_image(:)); %����������
figure(6); output_image = imadjust(output_image);
imshow(output_image, [ ]);

subtraction = orig_image - output_image;
pow = subtraction.^2; 
bfnch_err = sum(pow(:)) / number_of_pixels

%% ��������� � ������� ���������� �������
med = medfilt2(salt_pepper_noise);
med = medfilt2(med);
med = med -  min(med(:)); %����������
med = med/max(med(:)); %����������
figure(10); med = imadjust(med);
imshow(med, [ ]);

subtraction = orig_image - med;
pow = subtraction.^2; 
med_err = sum(pow(:)) / number_of_pixels