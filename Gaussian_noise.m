close all; clear all;

%% �������� �����������
orig_image = imread('8.png'); % ���������� ��������� �����������
orig_image = rgb2gray(orig_image); % ������� ����������� � �������. ���
orig_image = im2double(orig_image);
figure(1); imshow(orig_image,[]); % ����� ��������� �����������

%% ������� ��� ���� ��������� ����������
%Wavefuncs = {'sym2', 'sym3', 'sym4', 'sym5', 'sym6', 'sym7', 'sym8', 'bior 1.1', 'bior 1.3', 'bior 1.5', 'bior 2.2', 'bior 2.4', 'bior 2.6', 'bior 2.8', 'bior 3.1', 'bior 3.3', 'bior 3.5', 'bior 3.7', 'bior 3.9', 'bior 4.4', 'bior 5.5', 'bior 6.8', 'db2', 'db3', 'db4', 'db5', 'db6', 'db7', 'db10', 'db20', 'coif2',  'coif3',  'coif4',  'coif5', 'rbio1.1', 'rbio 1.3', 'rbio 1.5', 'rbio 2.2', 'rbio 2.4', 'rbio 2.6', 'rbio 2.8', 'rbio 3.1', 'rbio 3.3', 'rbio 3.5', 'rbio 3.7', 'rbio 3.9', 'rbio 4.4', 'rbio 5.5', 'rbio 6.8', 'dmey'}; % ��� ��������
Wavefuncs = {'rbio2.6'}; % ��� ��������
[si1, si2] = size(Wavefuncs);
N1 = 3; N2 = 1; % ������ ����������
Thr_min = 0.0; Thr_max = 0.7; % ������ ��� �������� ����
[n_size,m_size] = size(orig_image); % ������� �����������
number_of_pixels = n_size * m_size; % ���������� ��������

%% ����������� ���
gaus_noise = imnoise(orig_image,'gaussian',0,0.1); % ����������� ���
figure(2); imshow(gaus_noise,[]);

min_err = 1;
best_Thr = 0.0; % ������ �����
best_SorH = 's'; % ������ ��� ������
best_N=N1; % ������ ������� ����������
best_Wavefun = 'a';
for wave = 1:1.0:si2 % ������������� ������ ��������
    Wavefun = Wavefuncs{wave};
    for N = N1:-1.0:N2
        [Csp,Ssp]=wavedec2(gaus_noise,N,Wavefun); % �������-����������
        Wavefun
        Ssp % ����� ������� �������� �����, ������� ������ ���������� � �������� ������   
                % ���������������� � �������������� ������������� ���� �������
        for j = 0.0:1.0:1.0 % ��� ��������� ���������
            if (j ==0)
                SorH = 's';
            else
                SorH = 'h';
            end
            Thr = 0.0;
            for i = Thr_min:0.1:Thr_max % ����� ������� ������
                [no_gaus,Csp1,Ssp1,perf0,perfl2] = wdencmp('gbl',Csp,Ssp,Wavefun,N,i,SorH,1); % ��������
                    % ���� � ������� ���������� �������
                no_gaus = no_gaus -  min(no_gaus(:)); %����������
                no_gaus = no_gaus/max(no_gaus(:)); %����������
                no_gaus = imadjust(no_gaus); % �����-���������
                % ������� �������� �������� ������� ��������, �� ���������, ����� ��
                % ����� ��������, ������� ������
                subtraction = orig_image - no_gaus;
                pow = subtraction.^2;
                err = sum(pow(:)) / number_of_pixels;
                if (min_err > err) % ���� ����������� ������
                    min_err = err;
                    Thr = i;
                    best_Thr = i;
                    best_SorH = SorH;
                    best_N = N;
                    best_Wavefun = Wavefun;
                    figure(3); no_gaus = imadjust(no_gaus); imshow(no_gaus,[]);
                    %pause;
                    SorH
                    Wavefun
                    N
                    Thr
                    min_err
                end
            end
            if (Thr > 0.0)
                for i = Thr-0.1:0.01:Thr+0.1 % ����� ������� ������
                    [no_gaus,Csp1,Ssp1,perf0,perfl2] = wdencmp('gbl',Csp,Ssp,Wavefun,N,i,SorH,1); % ��������
                        % ���� � ������� ���������� �������
                    no_gaus = no_gaus -  min(no_gaus(:)); %����������
                    no_gaus = no_gaus/max(no_gaus(:)); %����������
                    no_gaus = imadjust(no_gaus); % �����-���������
                    % ������� �������� �������� ������� ��������, �� ���������, ����� ��
                    % ����� ��������, ������� ������
                    subtraction = orig_image - no_gaus; 
                    pow = subtraction.^2;
                    err = sum(pow(:)) / number_of_pixels;
                    if (min_err > err) % ���� ����������� ������
                        min_err = err;
                        best_Thr = i;
                        best_SorH = SorH;
                        best_N = N;
                        best_Wavefun = Wavefun;
                        figure(4); no_gaus = imadjust(no_gaus); imshow(no_gaus,[]);
                        %pause;
                        SorH
                        Wavefun
                        N
                        best_Thr
                        min_err
                    end
                end
            end
        end
    end
end
best_SorH
best_Thr
best_N
best_Wavefun
min_err
[Csp,Ssp]=wavedec2(gaus_noise,best_N,best_Wavefun); % �������-����������
[no_gaus,Csp1,Ssp1,perf0,perfl2] = wdencmp('gbl',Csp,Ssp,best_Wavefun,best_N,best_Thr,best_SorH,1);
figure(5); no_gaus = imadjust(no_gaus); imshow(no_gaus,[]);

%% ��������� � ������� ��������������� ������� �����������
FT_img = fft2(double(gaus_noise));
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
filtr_err = sum(pow(:)) / number_of_pixels