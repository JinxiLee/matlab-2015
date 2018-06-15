%%
clear all

%%
cd('D:\AMO\150915_Bubble identification');
Image_colored = imread('Y1.jpg','jpeg');

Image_Greyed=(sum(Image_colored,3));
imagesc(Image_Greyed);

Image_Green=double(Image_colored(:,:,2));
imagesc(Image_Green);
Image_Red=double(Image_colored(:,:,1));
imagesc(Image_Red);
Image_Blue=double(Image_colored(:,:,3));
imagesc(Image_Blue);

Image_Related=((Image_Red+Image_Green)./(Image_Blue));
Image_Related=Image_Related./max(Image_Related(:))*255;
Image_Related=(Image_Related);
imagesc(Image_Related);
axis equal
%% Special Gaussian for edge detection
Filter_Size=100;  %be Even
Filter_X=fspecial('disk', Filter_Size);    %gaussian mask
Filter_Y=fspecial('disk', Filter_Size);    %gaussian mask
Filter_X(1:(Filter_Size),:)=-Filter_X(1:(Filter_Size),:);
Filter_Y(:,1:Filter_Size)=-Filter_Y(:,1:Filter_Size);

imagesc(Filter_X);

%% Special edge detection
Image_Edge=(filter2(Filter_X,Image_Related).^2+filter2(Filter_Y,Image_Related).^2).^0.5;
Image_Edge(1:Filter_Size,:)=0;
Image_Edge(:,1:Filter_Size)=0;
Image_Edge((size(Image_Edge,1)-Filter_Size+1):end,:)=0;
Image_Edge(:,(size(Image_Edge,2)-Filter_Size+1):end)=0;
Image_Edge=Image_Edge/max(Image_Edge(:));
imagesc(Image_Edge);
%%
Grad_TH=0.01;

Ball_Size=1;
Element=fspecial('disk', Ball_Size);
Element(Element>0)=1;
Image_Eroded=imerode(Image_Related,Element);
Image_Dilated=(imdilate(Image_Related,Element));
Image_Grad=Image_Dilated-Image_Eroded;
Image_Grad=Image_Grad./max(Image_Grad(:));
Image_Grad_BW=Image_Grad;
Image_Grad_BW(Image_Grad>Grad_TH)=1;
Image_Grad_BW(Image_Grad<Grad_TH)=0;

imagesc(Image_Grad_BW);

%% Opening
Element2=fspecial('disk', Ball_Size*3);
Element2(Element2>0)=1;
Image_Opened=imopen(Image_Grad,Element2);

imagesc(Image_Opened);

%% BW 總之不行用BW
Image_Shed=watershed(Image_Edge);
imagesc(Image_Shed);


%% Top&bottom hat
Element_Size=1;

Element=fspecial('disk', Element_Size);    %gaussian mask
Element(Element>0)=1;
Image_Top=imtophat(Image_Related,Element);
imagesc(Image_Top);
Image_Top_BW=Image_Top;
Image_Top_BW(Image_Top_BW>0)=1;
imagesc(Image_Top_BW);
%%
Image_Bot=imbothat(Image_Related,Element);
imagesc(Image_Bot);

%% Distance map
Image_Top_BW=Image_Top;
Image_Top_BW(Image_Top_BW>0)=1;
imagesc(Image_Top_BW);

Image_Top_Dist=bwdist(Image_Top_BW);
imagesc(Image_Top_Dist);
%%
Th_H=0.06;
Th_L=0.03;
%Filter_1=[-1 0 1;-2 0 2; -1 0 1];
%Filter_2=[1 2 1;0 0 0; -1 -2 -1];
%Filter_3=[0 0 0;0 1 0; 0 1 0];

%Image_Filtered=uint8((filter2(Filter_1,Image_Sub,'same').^2+filter2(Filter_2,Image_Sub,'same').^2).^0.5);
Image_Edge=edge(Image_Related,'canny',[Th_L Th_H]);
imagesc(Image_Edge);
colormap(gray);


%% Mask
Level=0.07;
Image_Mask=im2bw(Image_Related,Level);

imagesc(Image_Mask);
colormap(gray);


%% Rolling ball
Ball_Size=200;
Ball=fspecial('gaussian', Ball_Size,Ball_Size/3);    %gaussian mask

Image_Green_Related_Rolled=uint8(filter2(Ball,Image_Green_Related,'same'));
Image_Green_Related_Sub=Image_Green_Related-Image_Green_Related_Rolled;
imagesc(Image_Green_Related_Sub);


%% Filling holes
Ball_Size=1;
Element=fspecial('disk', Ball_Size);
Element(Element>0)=1;
Image_Eroded=imerode(Image_Green_Related,Element);
imagesc(Image_Eroded);

Image_Dilated=(imdilate(Image_Eroded,Element));
imagesc(Image_Dilated);


%%
Image_Masked_Green=Image_Dilated.*Image_Greyed;
imagesc(Image_Masked_Green);


%% BW
Level=0.13;
Image_BW=im2bw(Image_Edge,Level);

imagesc(Image_BW);

%%
N=3;
Element=ones(N,N);
Image_Dilated=imdilate(Image_BW,Element);
imagesc(Image_Dilated);

%% Erode and Dilate

N=2;
Element=ones(N,N);
Image_Eroded=imerode(Image_Dilated,Element);
imagesc(Image_Eroded);

%%



%%
%Level=0.9999;
%Test=im2bw(Image_Filtered_2,Level);
Test_2=-bwdist(Image_Eroded);
Test_3=Test_2;
Test_4=watershed(Test_3);
%imagesc(~Test);
subplot(2,2,1);
imagesc(Image_Eroded);
axis equal
subplot(2,2,2);
imagesc(Test_2);
axis equal
subplot(2,2,3);
imagesc(Test_3);
axis equal
subplot(2,2,4);
imagesc(Test_4);
axis equal

colormap(gray);