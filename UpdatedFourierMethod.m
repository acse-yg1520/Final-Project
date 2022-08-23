clearvars;
close all;
clc;

img = imread('/Users/ashleyguan/Desktop/bubble-sizer-master/water/IMG_9380.JPG');  % Read the image
img = imresize(img, 0.5);      % Resize the image for faster processing
img = rgb2gray(img);           % Convert RGB image to grayscale
%imshow(img)

bknd_img = imread('/Users/ashleyguan/Desktop/bubble-sizer-master/water/IMG_9399.JPG');  % Read the background image
bknd_img = rgb2gray(bknd_img);      % Convert RGB image to grayscale
T = adaptthresh(bknd_img, 0.4 , 'ForegroundPolarity', 'dark'); %use background correction to help the gray threshold
BW = imbinarize(img,imresize(T, 0.5)); % Image binarization

%Use morphological operations to help make the bubbles solid
se = strel('disk', 5); %strel object to perform binary operations
B = imclose(~BW,se);
B = imfill(B,'holes');

% Remove small objects < 10 pixels
im = bwareaopen(B, 40);
% Operate along x-axis
im( ~any(im,2), : ) = [];  % Remove empty lines
[mx,nx] =size(im);
fft_x = fft(hann(mx).*(im-mean(im,2)),[],2);% FFT along the rows
FourierMean_x = 20*log10(mean(abs(fft_x).^2, 1)); % Taking average along the y axis
FourierMean_x = FourierMean_x - max(FourierMean_x);
FourierMean_x = FourierMean_x(:, 1:nx/2); % Half the dimension due to symmetry

% Operate along y-axis
im = bwareaopen(B, 40);
im( :, ~any(im,1) ) = [];  % Remove empty lines
[my,ny] = size(im);
fft_y = fft(hann(my).*(im-mean(im,1)),[],1);% FFT along the columns
FourierMean_y = 20*log10(mean(abs(fft_y).^2, 2)); % Taking average along the x axis
FourierMean_y = FourierMean_y - max(FourierMean_y);
FourierMean_y = FourierMean_y(1:my/2, :); % Half the dimension due to symmetry

Fs = 183;
freq_x = linspace(0, Fs/(2*0.5), nx/2);
freq_y = linspace(0, Fs/(2*0.5), my/2);

figure();
plot(freq_x, FourierMean_x);
axis([0 20 -20 0]);
title('FFT along x-axis');
xlabel('Frequency(pxl/mm)');
figure();
plot(freq_y,FourierMean_y);
axis([0 20 -40 0]);
title('FFT along y-axis');
xlabel('Frequency(pxl/mm)');

D32_x = 3.69./freq_x.^1.1;
figure;
plot(freq_x, D32_x);
axis([0 10 0 10])


D32_y = 3.69./freq_y.^1.1;
figure;
plot(freq_y, D32_x);
axis([0 10 0 10])
