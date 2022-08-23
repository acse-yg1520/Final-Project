%% Read the image

clear;
img = imread('/Users/ashleyguan/Desktop/bubble-sizer-master/water/IMG_9397.JPG');  % Read the image
img = imresize(img, 0.5);      % Resize the image for faster processing
img = rgb2gray(img);           % Convert RGB image to grayscale
%imshow(img)

%% We need the dark side as background and bright side as our target objects

opposite=ones(size(img,1),size(img,2))*255;
opposite=uint8(opposite);
opposite=opposite-img;
I=opposite;
I = imfill(I,'holes');

%% Reconstruct the image by opening-closing
se = strel('disk',20);
Ie = imerode(I, se);
Iobr = imreconstruct(Ie, I);

Iobrd = imdilate(Iobr, se);
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
figure;
%imshow(Iobrcbr);

%% Mark the foreground objects

fgm = imregionalmax(Iobrcbr);% Find regional maxima of Opening-Closing by Reconstruction

% Modify the regional maxima
se2 = strel(ones(5,5));
fgm2 = imclose(fgm, se2);
fgm3 = imerode(fgm2, se2);
fgm4 = bwareaopen(fgm3, 20)

%% Mark the background objects
bw=imbinarize(Iobrcbr,graythresh(Iobrcbr)); % Binarize the Opening-Closing by Reconstruction image

%% Compute the watershed ridge line
D = bwdist(bw);
DL = watershed(D);
bgm = DL == 0;

%% Compute the gradient magnitude as segmentation function

% Use prewitt operator as edge detection
hy = fspecial('prewitt');
hx = hy';
Iy = imfilter(double(I), hy, 'replicate');
Ix = imfilter(double(I), hx, 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);

%% Compute watershed sgementation based on the gradient magnitude

gradmag2 = imimposemin(gradmag, bgm | fgm);
L = watershed(gradmag2);
bw(L==0)=0;
imshow(bw);
