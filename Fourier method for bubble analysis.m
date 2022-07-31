img = imread('IMG_9380.JPG');  % Read the image
img = imresize(img, 0.5);      % Resize the image for faster processing
img = rgb2gray(img);           % Convert RGB image to grayscale
%imshow(img)

bknd_img = imread('IMG_5784.JPG');  % Read the background image
bknd_img = rgb2gray(bknd_img);      % Convert RGB image to grayscale
T = adaptthresh(bknd_img, 0.4 , 'ForegroundPolarity', 'dark'); %use background correction to help the gray threshold
BW = imbinarize(img,imresize(T, 0.5)); % Image binarization

%Use morphological operations to help make the bubbles solid
se = strel('disk', 5); %strel object to perform binary operations
B = imclose(~BW,se);
B = imfill(B,'holes');

% Process image size information
[M N] = size(B); 

% Plot the pulse train obtained from the image
plot(B(101,:))

%% Compute normalized power spectrum
img_power = fftshift(fft2(B));
img_power_nor = (abs(img_power)/(M*N)).^2; % compute normalised psd                                              % Normalize

%% Apply zero padding to adjust non-square dimensions
diff = abs(M-N);  % difference of rows and columns numbers
maxD = max(M,N);  % maximum of dimensions

if M > N                           % More rows than columns
    if (mod(diff,2) == 0)          % Even difference
        imgB =  [zeros(M, diff/2) img_power_nor zeros(M, diff/2)]        % Add columns to match dimensions
    else                           % Odd difference
        imgB = [zeros(M, floor(diff/2)) img_power_nor zeros(M, floor(diff/2) + 1)];
    end
elseif M < N                       % More columns than rows
    if (mod(diff,2) == 0)          % Even difference
        imgB = [zeros(diff/2, N); img_power_nor; zeros(diff/2, N)];         % Add rows to match dimensions
    else
        imgB = [zeros(floor(diff/2), N); img_power_nor; zeros(floor(diff/2) + 1, N)];
    end
end
                                     
%% Compute radially average power spectrum
halfDim = floor(maxD/2) + 1;  % Consider half the dimensions due to symmetry
pf = zeros(halfDim, 1);
count = zeros(halfDim, 1);

for i = 1 : maxD
  for j = 1 : maxD
    radius = sqrt((i - halfDim) ^ 2 + (j - halfDim) ^ 2);
    Index = ceil(radius) + 1;
    if Index <= halfDim
       pf(Index) = pf(Index) + imgB(i, j);
       count(Index) = count(Index) + 1;
    end
  end
end
pf = pf./ count; % Get average

%% Plot the normalised psd
figure; semilogy(pf)
title('Average Radial Profile', 'FontSize', 10);
xlabel('Frequency', 'FontSize', 10);
ylabel('Normalised PSD', 'FontSize', 10);

% Model D32 vs BW
% D32 = 3.7./(sf.^1.1);
