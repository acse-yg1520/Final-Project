%Template image processing algorithm coming with Bubble Analyser software
%
% Syntax: FindCircles(I, params)
%
% Inputs:
%    I: image, either rgb or grayscale
%    params: parameters needed for certain operations, either set by the
%    user or in the corresponding .config file
%    
% Outputs:
%    D: number array with the equivalent diameter, in mm, of each bubble detected and segmented
%    L: labelled image resulting from the image processing algorithm
%    extra_info: structure containing extra information about bubbles
%    (eccentricity, solidity, etc)
%
%
% Author: Reyes, Francisco; Quintanilla, Paulina; Mesa, Diego
% email: f.reyes@uq.edu.au,  
% Website: https://gitlab.com/frreyes1/bubble-sizer
% Copyright Feb-2021;
%
%This file is part of Bubble Analyser.
%
%    Bubble Analyser is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation version 3 only of the License.
%
%    Bubble Analyser is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with Bubble Analyser. If not, see <https://www.gnu.org/licenses/>.
%
%------------- BEGIN CODE --------------
function [D, L_image, extra_info] = SplitPointsDetection(img, params)

%Collect parameters from the structure
se = strel('disk', params.Morphological_element_size); %strel object to perform binary operations;
nb = params.Connectivity; %neighbourhood used (4 or 8)
marker_size = params.Marker_size; %marker size of the Watershed segmentation, in px
px2mm = params.px2mm; %img resolution
img_resample = params.resample;
bknd_img = params.background_img;
convexity = params.convexity;
eccentricity = params.eccentricity;
solidity = params.solidity;
E_max = params.Max_Eccentricity;
S_min = params.Min_Solidity;
Dmin = params.min_size; %minimum bubble size, in mm!
do_batch = params.do_batch; %Check if we are doing batch processing or justone image

%Resize images for making processing faster
[n, m, k] = size(img);
img = imresize(img, img_resample); %resample img to make process faster
if k>1
    img = rgb2gray(img);
end

%use background correction to help the gray threshold
if ~isempty(bknd_img)
    if size(bknd_img,3)>1
        bknd_img = rgb2gray(bknd_img);
    end
    T = adaptthresh(bknd_img, 0.4 , 'ForegroundPolarity', 'dark');
end
if ~isempty(bknd_img)
    BW = imbinarize(img,imresize(T, img_resample));
else
    BW = imbinarize(img);
end


%Use morphological operations to help make the bubbles solid
B = imclose(~BW,se);
B = imfill(B,'holes');


%Now use watershed to separate bubbles that are overlapping
R = -bwdist(~B);
% mask = imextendedmin(R,nb);
% R2 = imimposemin(R,mask);
R2 = imhmin(R,marker_size,nb); %J = imhmin(I,H,conn) computes the H-minima transform, where conn specifies the connectivity.
Ld2 = watershed(R2);
B(Ld2 == 0) = 0;
CH = bwconvhull(B,'objects');
R3 = -bwdist(~CH);
mask = imextendedmin(R3,nb);
R4 = imimposemin(R3,mask);
Ld3 = watershed(R4,nb);
CH(Ld3==0) = 0;


%Eccentricity: 0 -> circle, 1 -> line; Solidity = Area/ConvexArea
S = regionprops(CH,'EquivDiameter','Perimeter', 'ConvexImage','Solidity','Eccentricity');
boundaries = bwboundaries(CH);

for i = 1:length(S)
    ConvexPerim = regionprops(S(i).ConvexImage, 'Perimeter');
    if (S(i).Perimeter/ConvexPerim.Perimeter) > convexity && (S(i).Solidity < solidity) && (S(i).Eccentricity >= eccentricity)
         Boundary = boundaries{i}; 
	     x = Boundary(:,2); 
	     y = Boundary(:,1); 
         new_img = separation(x,y,S(i).Perimeter,CH);
         CH = new_img;
    end 

end

CH = imclearborder(CH);
%Now list the detected objects and calculate geometric properties
CC = bwconncomp(CH,4);
%Eccentricity: 0 -> circle, 1 -> line; Solidity = Area/ConvexArea
new_S = regionprops(CC,'EquivDiameter','Eccentricity','Solidity');

%Reject abnormal objects, possibly unseparated bubbles
E = [new_S.Eccentricity]'; %column vector with eccentricity
D = [new_S.EquivDiameter]'; %column vector with diameters
S = [new_S.Solidity]';
%!!Remember we scaled down the image by some factor!!
D = D * px2mm * 1/img_resample; %now in mm
idx = E>=E_max | S<=S_min | D<Dmin; %abnormal bubbles: too stretched
%remove abnormal bubbles
D = D(~idx);
%collect extra bubble shape descriptors
extra_info.Eccentricity = E(~idx);
extra_info.Solidity = S(~idx);
        
%Update label image if required. 
if do_batch
    %when doing batch processing we don't need to create a fancy label image 
    L_image = [];
else
    %when processing individual images we can create a nice label image to
    %show the results to the user
    allowableAreaIndexes = ~idx;
    keeperIndexes = find(allowableAreaIndexes);
    keeperBlobsImage = ismember(bwlabel(CH,4), keeperIndexes);
    keeperBlobsImage = imresize(keeperBlobsImage,[n m]); %put it back to original size
    L_image = bwlabel(keeperBlobsImage, nb);

end
%% Separate the overlapping objects by connecting the concave points
function splited = separation(x,y,perimeter,img)

% Find the candidate split points
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

% Determine the exact concave points
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


% Use DDA method to connect two split points
posx = cell2mat(dic(index,1)); posy = cell2mat(dic(index,2));

if posx(1) < posx(2)
    x1 = posx(1); y1 = posy(1); x2 = posx(2); y2 = posy(2);
elseif posx(1) > posx(2)
    x1 = posx(2); y1 = posy(2); x2 = posx(1); y2 = posy(1);
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
splited=img;
for i = 1:length(X_n)
    i = X_n<1; X_n(i) = []; Y_n(i) = [];
    i = X_n>n; X_n(i) = []; Y_n(i) = [];
    i = Y_n<1; X_n(i) = []; Y_n(i) = [];
    i = Y_n>m; X_n(i) = []; Y_n(i) = [];
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

end


