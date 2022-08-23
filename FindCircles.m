
img = imread('/Users/ashleyguan/Desktop/bubble-sizer-master/water/IMG_9386.JPG');  % Read the image
img = imresize(img, 0.5);      % Resize the image for faster processing
img = rgb2gray(img);            % Convert RGB image to grayscale
%imshow(img)

% Piecewise Linear Gray Level Transformation
x1 = img
[L1,L2]=size(x1);
a=100;
b=200;
mf=255;
c=50;
d=230;
mg=255;
f=double(x1);
num1=zeros(1,256);
num2=zeros(1,256);
for i=1:L1
    for j=1:L2
        num1(1,f(i,j))=num1(1,f(i,j))+1;
        if(f(i,j)<a&f(i,j)>=0)
            g(i,j)=round(c*f(i,j)/a);
        elseif(f(i,j)<b&f(i,j)>=a)
            g(i,j)=round(((d-c)/(b-a))*(f(i,j)-a)+c);
        elseif(f(i,j)<=mf&f(i,j)>=b)
            g(i,j)=round(((mg-d)/(mf-b))*(f(i,j)-b)+d);
        end
        num2(1,g(i,j))=num2(1,g(i,j))+1;
    end
end
img = uint8(g);

% Find circles and locate them by centers and radii.
[centers_1, radii_1, metric_1] = imfindcircles(img, ...       
        [5,30], ...                                              % Radius range in pixels.
        'ObjectPolarity', 'dark', ...                            % Find circular objects that are darker than the background, [brighter, darker].
        'method', 'twostage', ...                                % Two methods for finding circles, [phase coding, twostage].
        'sensitivity', 0.8, ...                                  % 'Sensitivity', [0,1], is set to 0.85 by default.
       'edgethreshold', 0.3)                               
                                

[centers_2, radii_2, metric_2] = imfindcircles(img, ...       
        [30,500], ...                                            % Radius range in pixels.
        'ObjectPolarity', 'dark', ...                            % Find circular objects that are darker than the background, [brighter, darker].
        'Method','TwoStage')                                    % A high value (closer to 1) will allow only the strong edges to be included, whereas a low value (closer to 0) includes even the weaker edges.

imshow(img); hold on;
viscircles(centers_1, radii_1, 'edgecolor', 'r'); 
viscircles(centers_2, radii_2, 'edgecolor', 'r');   

