
clear;
img = imread('/Users/ashleyguan/Desktop/bubble-sizer-master/frother/IMG_7587.JPG');  % Read the image
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


[bool,idx,props, boundaries] = ifoverlap(B);
while bool >= 1

  for i = idx
        thisBoundary = boundaries{i};
	    x = thisBoundary(:,2);
	    y = thisBoundary(:,1);
        perimeter = props(i).Perimeter
        splited = separation(x,y,perimeter,B)
        B = splited;
        [bool,idx,props,boundaries] = ifoverlap(B);
  end
end

%% Detect the objects are overlapping or not
function [bool, idx, props, boundaries] = ifoverlap(img)

bool = 0;
idx = [];
n = 1;
boundaries = bwboundaries(img);
numberOfBoundaries = size(boundaries, 1);
props = regionprops(img, 'Centroid', 'Perimeter', 'Solidity','Eccentricity')
for i = 1:numberOfBoundaries
    if props(i).Solidity <= 0.80 && props(i).Eccentricity >= 0.8
        bool = bool + 1;
        idx(n) = i;
        n = n+1;
    end
end

end



%% Separate the overlapping objects by connecting the contour points
function splited = separation(x,y,perimeter,img)

% Find the candidate concave points for split
gap = floor(length(x)/50);
X = x(1:gap:end,:);
Y = y(1:gap:end,:);
cout = [Y,X];
n = 1;
contourp = []; k = 1;
for i = 1:length(X)
    xn = []; yn = [];
    m = 1;
    for theta = 0:5:360
       xn(m) = cout(i,2) + 20*cosd(theta);
       yn(m) = cout(i,1) + 20*sind(theta);
       m = m+1;
    end
    in = inpolygon(xn,yn,X,Y);
    thred = length(find(in~=0))/m;
    if thred >= 0.6
        contourp(k,2) = cout(i,2);
        contourp(k,1) = cout(i,1);
        k = k+1
    end
end

% Determine the exact split points
dic = {};
k = 1;
for i = 1:length(contourp)-1
    for j = i+1: length(contourp)
        %perimeter = props(3).Perimeter;
        p = perimeter/length(X);
        idx = find([cout(:,2)] == contourp(i,2));
        jdx = find([cout(:,2)] == contourp(j,2));

        dist1 = abs((jdx - idx)*p); dist2 = abs(length(X)-max(idx,jdx)+abs(min(idx,jdx)))*p; dist = min(dist1,dist2);
        eu_dist = sqrt((cout(jdx,1)-cout(idx,1))^2 + (cout(jdx,2)-cout(idx,2))^2);
        dic(k,:) = {[contourp(i,2),contourp(j,2)],[contourp(i,1),contourp(j,1)], eu_dist/dist}
        k = k+1
    end
end
index = find([dic{:,3}] == min( cell2mat(dic(1:end,3))));


% Connect the contour points
posx = cell2mat(dic(index,1)); posy = cell2mat(dic(index,2));

if posx(1) < posx(2)
    x1 = posx(1); y1 = posy(1); x2 = posx(2); y2 = posy(2)
elseif posx(1) > posx(2)
    x1 = posx(2); y1 = posy(2); x2 = posx(1); y2 = posy(1)
end

L=max(abs(posx(2)-posx(1)),abs(posy(2)-posy(1)));
k=(y2-y1)/(x2-x1);
X_n=[];Y_n=[];

for i=1:L+1
       X_n(i)=x1;
       Y_n(i)=y1;
       x1=x1+1;
       y1=y1+k;
end
Y_n=fix(Y_n+0.5);


% Get the splited image
splited=B;
for i=1:length(X_n)
    splited(Y_n(i),X_n(i))=0;
    splited(Y_n(i)-1,X_n(i))=0;
    splited(Y_n(i)+1,X_n(i))=0;
    splited(Y_n(i),X_n(i)-1)=0;
    splited(Y_n(i),X_n(i)+1)=0;
    splited(Y_n(i)-1,X_n(i)-1)=0;
    splited(Y_n(i)+1,X_n(i)+1)=0;
    splited(Y_n(i)-1,X_n(i)+1)=0;
    splited(Y_n(i)+1,X_n(i)-1)=0;
end

end
