clear all
close all
clc

% number of images on your training set.
num=10;

%read and show images;
IMG=[];   %img matrix
Im = [];

for j=1:10
    for i=1:10
        filename = strcat('s',num2str(j),'\',num2str(i),'.pgm');
        Re = imread(filename); % reading the image into Im
        %Re(i,:,j) = reshape(Im{i,j},[],1);  % reshape to 10x10304 of 10 datasets
        subplot(ceil(sqrt(num)),ceil(sqrt(num)),i)
        imshow(Re)
        Imz(j,i,:,:)=Re(:,:);
        annotation('textbox', [0 0.9 1 0.1], ...
            'String', 'Training', ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center')
    end
    drawnow;
    [irow, icol]=size(Re);    % get the number of rows (N1) and columns (N2)
    temp=reshape(Re',irow*icol,1);
    IMG=[IMG temp];
end

old_m=mean2(IMG);
old_std=std2(IMG);

%mean image;
MeanIm(:,:) = uint8(mean(mean(Imz)));
figure(3);
imshow(MeanIm);
title('Mean Image')

m = mean2(MeanIm);
st = std2(MeanIm);

%Covariance matrices:
A=double(IMG');
L=A*A';
% evr are the eigenvector for L
% evl are the eigenvalue for both L=IMG'*IMG and C=IMG*IMG';
[evr, evl]=eig(L);

d = sort(evl(evl~=0)','descend');
v = fliplr(evr);
 %sort,  will return an ascending sequence
 [B, index]=sort(d);
 ind=zeros(size(index));
 dtemp=zeros(size(index));
 vtemp=zeros(size(v));
 len=length(index);
 for i=1:len
    dtemp(i)=B(len+1-i);
    ind(i)=len+1-index(i);
    vtemp(:,ind(i))=v(:,i);
 end
 d=dtemp;
 v=vtemp;
u=[];
%Normalization of eigenvectors
for i=1:size(v,2)
   v(:,i)=v(:,i)./sqrt(sum(v(:,i).^2));
   temp=sqrt(d(i));
   u=[u (double(IMG)*v(:,i))./temp];
end

%Normalization of eigenvectors
for i=1:size(u,2)
   u(:,i)=u(:,i)./sqrt(sum(u(:,i).^2));
end

% show eigenfaces;
figure;
for i=1:size(u,2)
    Im=reshape(u(:,i),icol,irow);
    Im=Im';
    Im=histeq(Im,255);
    subplot(ceil(sqrt(num)),ceil(sqrt(num)),i)
    imshow(Im)
    drawnow;
    annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Eigenfaces', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
end

%Recognizing
[filename, pathname] = uigetfile({'*.pgm'},'Input Image Selector');
InputImage = imread(strcat(pathname,filename));
figure(4)
subplot(1,2,1)
imshow(InputImage); 
title('Input image')
InImage=reshape(double(InputImage)',irow*icol,1);  
temp=InImage;
me=mean(temp);
st=std(temp);
temp=(temp-me)*old_std/st+old_m;
NormImage = temp;
Difference = temp-m;

p = [];
aa=size(u,2);
for i = 1:aa
    p = [p; dot(NormImage,u(:,i))];
end
ReshapedImage = m + u(:,1:aa)*p;    %m is the mean image, u is the eigenvector
ReshapedImage = reshape(ReshapedImage,icol,irow);
ReshapedImage = ReshapedImage';
%show the reconstructed image.
subplot(1,2,2)
imshow(uint8((ReshapedImage))); 
title('Reconstructed image')