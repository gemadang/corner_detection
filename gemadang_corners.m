%{
Hints:
- For eigenvalue computation: help eig
- For sorting: help sort, help sortrows
%}
% Geri Madanguit
% READ ME: In this code, change the value of A to b or c depending on what
% picture you want to view

clear all;
clc;

b = imread('Building1.jpg');
c = imread('CheckerBoard.jpg');
c = rgb2gray(c);

A = c;

%figure, imshow(A);

% First, apply Gaussian smoothing (with standard deviation ?) 
%   to an input image I, to obtain Is

kerG2r = 3; %gaussian mask row
kerG2c = 3; %gaussian mask column
std = 2.5; %sigma
    
centerx = round((kerG2r+1)/2);
centery = round((kerG2c+1)/2);
    
int_gauss_ker = ones(kerG2r,kerG2c); 
float_ker = zeros(kerG2r,kerG2c);

 for k=1:kerG2c
    for h=1:kerG2r
        float_ker(h,k) = exp((-(((h-centerx)^2+(k-centery)^2)/2*(std)^2)));
    end
 end
    
gmin = float_ker(centerx,centery);
norm_factor = 1/gmin;
int_gauss_ker = round(norm_factor*float_ker);
gauss_filtered_im = conv2(int_gauss_ker,A,'full');

max_gauss_filtered = max(max(gauss_filtered_im));

gauss_filtered_im = uint8((gauss_filtered_im/max_gauss_filtered)*255);

%figure, imshow(gauss_filtered_im);

% Implement the corner detection algorithm (CORNERS), by using Is as input,
%   as described in class

%input: image, two parameters: threshold on lambda2, tau, and linear size
%of a square window (neighborhood), says 2N+1

%1. compute image gradient over entire image I

ex_vector = [-1,0,1];
ey_vector = ex_vector';

[row, col] = size(gauss_filtered_im);

for r=1:row
    ex(r,:) = conv(ex_vector, gauss_filtered_im(r,:));
end

for c=1:col
    ey(:,c) = conv(ey_vector, gauss_filtered_im(:,c));
end

%2. for each image point p:
%(a) form the matrix C of (4.9) over a (2N+1) * (2N+1) neighborhood Q of
%   p;
%(b) compute lambda2, the smallest eigenvalue of C;
%(c) if lambda2>tau, save the coordinates of p into a list, L

th = 70000;
N = 5;  %N=1 --> 3x3
        %N=2 --> 5x5
        %N=3 --> 7x7
        %N=4 --> 9x9

neighbor_size = 2*N+1;

r_eig = 1;
%L = zeros(:,3);
[row_p, col_p] = size(gauss_filtered_im);
for r=N+1:row_p-N
    for c=N+1:col_p-N
        [isEigen, eigen_val] = get_eigen(th,r,c,ex,ey,N);
        if(isEigen)
            %add to matrix_list
            L(r_eig,1) = r;
            L(r_eig,2) = c;
            L(r_eig,3) = eigen_val; 
            r_eig = r_eig+1;
        end            
    end
end

%3. Sort L in decreaseing order of lambda2
sorted_L = sortrows(L,-3); %sortrows in descending order

%4. Scanning the sorted list top to bottom: for each current point, p,
%delete all points appearing further on in the list which belong to the
%neighborhood of p.


flag = 0;
new_r = 1;
flag_mat = zeros(size(sorted_L));
[L_rows, ~] = size(sorted_L);

for r=1:L_rows-1
    
    x_hi = sorted_L(r,1) + N;
    x_lo = sorted_L(r,1) - N;
    y_hi = sorted_L(r,2) + N;
    y_lo = sorted_L(r,2) - N;
    
    flag = flag_mat(r,3);   %flag = 1 when there is redundancy
        
    if(~flag)
        new_L(new_r,1) = sorted_L(r,1);
        new_L(new_r,2) = sorted_L(r,2);
        new_L(new_r,3) = sorted_L(r,3); 
        new_r = new_r+1;
    end
        
    for a=r+1:L_rows
        if(sorted_L(a,1) <= x_hi && sorted_L(a,1) >= x_lo)
            if(sorted_L(a,2)<= y_hi && sorted_L(a,2) >= y_lo)
                flag_mat(a,3) = 1;
            end
        end
    end
    
end

%The output is a list of feature points for which lamba2>tau and whose
%neighborhoods do not overlap

[nR, ~] = size(new_L);
for ro=1:nR
    corners(ro,1) = new_L(ro,1); 
    corners(ro,2) = new_L(ro,2);
end

output_im = gauss_filtered_im;
for r=1:nR
    output_im(corners(r,1),corners(r,2)) = 255;
end

output_im = gray2rgb(output_im);
%imshow(uint8(output_im));
%hold on;

%{
for r=1:nR
    output_im(corners(r,1),corners(r,2),3) = 150;
    %plot(corners(r,1),corners(r,2),'r+', 'MarkerSize', 5, 'LineWidth', 3);
end
%}

imshow(output_im);
% Test your corner detection algorithm on images ?Building1.jpg? and 
%   ?CheckerBoard.jpg?. Try different values of the ?, 
%   the neighborhood size, and the threshold (?) on ?2. 
%   Compare and evaluate your results.
function [Image]=gray2rgb(Image)
%Gives a grayscale image an extra dimension
%in order to use color within it
[m n]=size(Image);
rgb=zeros(m,n,3);
rgb(:,:,1)=Image;
rgb(:,:,2)=rgb(:,:,1);
rgb(:,:,3)=rgb(:,:,1);
Image=rgb/255;
end
function [output, ld2] = get_eigen(thresh,row, col,e_x,e_y,n_val)
    mat_C = zeros(2);
  
    for r=row-n_val:row+n_val
        for c=col-n_val:col+n_val
            mat_C(1,1) = mat_C(1,1) + e_x(r,c)*e_x(r,c);
            mat_C(1,2) = mat_C(1,2) + e_x(r,c)*e_y(r,c);
            mat_C(2,1) = mat_C(2,1) + e_x(r,c)*e_y(r,c);
            mat_C(2,2) = mat_C(2,2) + e_y(r,c)*e_y(r,c);
        end
    end
    e = eig(mat_C);
    
    lambda = e(2);
    %check smaller eig value
    if (e(2)>e(1))
        lambda = e(1);
    end
    
    output = (lambda > thresh);
    ld2 = lambda;
end
