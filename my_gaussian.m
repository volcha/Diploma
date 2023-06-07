close all; clear all;

%% �������� �����������
orig_image = imread('8.png'); % ���������� ��������� �����������
orig_image = rgb2gray(orig_image); % ������� ����������� � �������. ���
orig_image = im2double(orig_image);
figure(1); imshow(orig_image,[]); % ����� ��������� �����������

%% ������� ���
gaus_noise = imnoise(orig_image,'gaussian',0,0.1); % ������� ���
figure(2); imshow(gaus_noise,[]);

%% ������� ��� ���� ��������� ����������
%Wavefuncs = {'haar', 'sym2', 'sym3', 'sym4', 'sym5', 'sym6', 'sym7', 'sym8', 'bior 1.1', 'bior 1.3', 'bior 1.5', 'bior 2.2', 'bior 2.4', 'bior 2.6', 'bior 2.8', 'bior 3.1', 'bior 3.3', 'bior 3.5', 'bior 3.7', 'bior 3.9', 'bior 4.4', 'bior 5.5', 'bior 6.8', 'db2', 'db3', 'db4', 'db5', 'db6', 'db7', 'db10', 'db20', 'coif2',  'coif3',  'coif4',  'coif5', 'rbio1.1', 'rbio 1.3', 'rbio 1.5', 'rbio 2.2', 'rbio 2.4', 'rbio 2.6', 'rbio 2.8', 'rbio 3.1', 'rbio 3.3', 'rbio 3.5', 'rbio 3.7', 'rbio 3.9', 'rbio 4.4', 'rbio 5.5', 'rbio 6.8', 'dmey'}; % ��� ��������
Wavefuncs = {'bior2.6'}; % ��� ��������
[si1, si2] = size(Wavefuncs);
N = 1; % ������� ����������
[n_size,m_size] = size(orig_image); % ������� �����������
number_of_pixels = n_size * m_size; % ���������� ��������

min_err = 1;
best_wave = 'a';
for wave = 1:1.0:si2 % ������������� ������ ��������
    Wavefun = Wavefuncs{wave};
    [C,S] = wavedec2(gaus_noise,N,Wavefun); % ������������
    [mc,nc]=size(C);
    Cn = C;
    S
    nk = 0;
    for j = 1:N
        nk = nk + S(N+2-j,1) * S(N+2-j,2) * 3;
        nf = nc - nk;
        Cn(nf:end) = zeros(1,nc-nf+1);
        no_gaus = waverec2(Cn,S,Wavefun); % �������������
        no_gaus = no_gaus -  min(no_gaus(:)); %����������
        no_gaus = no_gaus/max(no_gaus(:)); %����������
        no_gaus = imadjust(no_gaus); 
        figure(4); imshow(no_gaus,[]);
        subtraction = orig_image - no_gaus;
        pow = subtraction.^2;
        mse = sum(pow(:)) / number_of_pixels
        %pause;
        if (min_err > mse) % ����� ����������� ������
            min_err = mse;
            best_wave = Wavefun;
        end
    end
    Wavefun
    %pause;
end
min_err
best_wave

%% ��������� � ������� ��������������� ������� �����������
% ��������� �������������� ����� ����������� �����������
% ������������� ������������ ������� MATLAB fft2 (2D ������� �������������� �����)
FT_img = fft2(double(gaus_noise));
% Assign the order value
n = 2; % ����� �� ���� �����
 
% ����� ������� �����
D0 = 200; % ��� �������� ����� ��������
 
% �������� �������
u = 0:(n_size-1);
v = 0:(m_size-1);
idx = find(u > n_size/2);
u(idx) = u(idx) - n_size;
idy = find(v > m_size/2);
v(idy) = v(idy) - m_size;

% ������� ���������� MATLAB meshgrid(v, u) ����������
% 2D-�����, ���������� ���������� ��������
% v � u. ������� V, ������ ������ ������� �������� ������ v
% � ������� U, ��� ������ ������� �������� ������ u
[V, U] = meshgrid(v, u);

% ���������� ��������� ����������
D = sqrt(U.^2 + V.^2);

% ����������� ����� ����������
H = 1./(1 + (D./D0).^(2*n));
  
% ������� ����� ���������������� �����
% ����������� � �����
G = H.*FT_img;

% ��������� ��������������� ����������� �������� ��������������� �����
% �� ���������� ����������� � �������������� ������������ ������� MATLAB
% ifft2 (2D �������� ������� �������������� �����)
output_image = real(ifft2(double(G))); 

% ����������� ����������� ����������� � ���������������� �����������
output_image = output_image -  min(output_image(:)); %����������
output_image = output_image/max(output_image(:)); %����������
figure(6); output_image = imadjust(output_image);
imshow(output_image, [ ]);

subtraction = orig_image - output_image;
pow = subtraction.^2; 
filtr_mse = sum(pow(:)) / number_of_pixels