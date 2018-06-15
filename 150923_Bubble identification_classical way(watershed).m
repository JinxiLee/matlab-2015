%%
clear all

%%
cd('D:\AMO\150915_Bubble identification');
Image_colored = imread('Y2.jpg','jpeg');

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

%% Grad

%Ball_Size=1;
%Element=fspecial('disk', Ball_Size);
%Element(Element>0)=1;
%Image_Eroded=imerode(Image_Related,Element);
%Image_Dilated=(imdilate(Image_Related,Element));
%Image_Grad=Image_Dilated-Image_Eroded;
%Image_Grad=Image_Grad./max(Image_Grad(:));

Grad_X=diff(Image_Greyed,1,1);
Grad_X(2:size(Grad_X,1)+1,:)=Grad_X(1:size(Grad_X,1),:);
Grad_Y=diff(Image_Greyed,1,2);
Grad_Y(:,2:size(Grad_Y,2)+1)=Grad_Y(:,1:size(Grad_Y,2));

Image_Grad=(Grad_X.^2+Grad_Y.^2).^0.5;

Image_Edge=double(edge(Image_Greyed,'canny'));

imagesc(Image_Edge);
colormap(gray);
%imagesc(Grad);
%% Opening
TH=0.8;
Ball_Size=10;
pixel_TH=5000;
Element2=fspecial('disk', Ball_Size*3);
Element2(Element2>0)=1;
Image_Opened=imopen(Image_Related,Element2);
Image_Opened_inv=imcomplement(Image_Opened);
Image_Opened_inv=(Image_Opened_inv-min(Image_Opened_inv(:)))/(max(Image_Opened_inv(:))-min(Image_Opened_inv(:)));
Image_Opened_inv_BW=im2bw(Image_Opened_inv,TH);
Image_Opened_inv_BW_final=bwareaopen(Image_Opened_inv_BW,pixel_TH);
imagesc(Image_Opened);
%% Erode (讓marker略小一點)
Ball_Size=100;
Element=fspecial('disk', Ball_Size);
Element(Element>0)=1;
Image_Marker=imerode(Image_Opened_inv_BW_final,Element);
%Image_Marker=bwulterode(Image_Opened_inv_BW_final,'chessboard');
imagesc(Image_Marker);

%Image_Marker=zeros(size(Image_Marker,1),size(Image_Marker,2));
%Image_Marker(735:745,885:895)=1;
%imagesc(Image_Marker);

%% to generate BGD marker
Image_Marker_BGD=ones(size(Image_Marker,1),size(Image_Marker,2));
Image_Marker_BGD(2:(size(Image_Marker,1)-1),2:(size(Image_Marker,2)-1))=0;

imagesc(Image_Marker_BGD);
Image_Marker_all=Image_Marker | Image_Marker_BGD;
imagesc(Image_Marker_all);


%% reconstruction & watershed

Image_Grad_2=imimposemin(Image_Grad,Image_Marker_all);
imagesc(Image_Grad_2);

Image_Result=watershed(Image_Grad_2);

imagesc(Image_Result);
axis equal

%%
Image_Show=Image_Greyed;
Image_Show(Image_Result==0)=max(Image_Greyed(:));

imagesc(Image_Show);
axis equal

%% Identify points in the line

Xgrid(1:size(Image_Result,1),1:size(Image_Result,2))=0;
Ygrid(1:size(Image_Result,1),1:size(Image_Result,2))=0;

for p=1:size(Image_Result,2)
    Xgrid(:,p)=p;
end
for q=1:size(Image_Result,1)
    Ygrid(q,:)=q;
end
Xgrid_Sample=Xgrid;
Ygrid_Sample=Ygrid;

Xgrid_Sample(Image_Result>0)=[];
Ygrid_Sample(Image_Result>0)=[];

%% 一個很神的circular fit (Izhak Bucher's) 類似的概念可以用到其他fitting嗎? 

x=Xgrid_Sample(:); y=Ygrid_Sample(:);
a=[x y ones(size(x))]\[-(x.^2+y.^2)];
xc = -.5*a(1);
yc = -.5*a(2);
R  =  sqrt((a(1)^2+a(2)^2)/4-a(3));

%%

imagesc(Image_Show);
axis equal
hold on
viscircles([xc,yc],R,'EdgeColor','r');     %X和Y在顯示上會是反的 所以只好改顯示方向!
hold off


%%
ratio=0.37;
R_hole=ratio*R;

imagesc(Image_Show);
axis equal
hold on
viscircles([xc,yc],R_hole,'EdgeColor','r');     %X和Y在顯示上會是反的 所以只好改顯示方向!
hold off

%%
Mask_Hole=zeros(size(Image_Marker,1),size(Image_Marker,2));

Mask_Hole(((Xgrid-xc).^2+(Ygrid-yc).^2).^0.5<R_hole)=1;

imagesc(Mask_Hole);
axis equal
%% To find the marker inside

Image_Inside_Mask=(max(Image_Greyed(:))-Image_Greyed);
Image_Inside_Mask=abs(Image_Inside_Mask-mean(mean(Image_Inside_Mask(((Xgrid-xc).^2+(Ygrid-yc).^2).^0.5<R_hole))));
Image_Inside_Mask=Image_Inside_Mask.*Mask_Hole;
Image_Inside_Mask(Image_Inside_Mask<mean(mean(Image_Inside_Mask(((Xgrid-xc).^2+(Ygrid-yc).^2).^0.5<R_hole)))*2)=0;
imagesc(Image_Inside_Mask);

Center_X_Bubble=round(sum(sum(Xgrid.*Image_Inside_Mask))/(sum(Image_Inside_Mask(:))));
Center_Y_Bubble=round(sum(sum(Ygrid.*Image_Inside_Mask))/(sum(Image_Inside_Mask(:))));


%%
Marker_Size=10;
Image_Inside_Marker=zeros(size(Image_Marker,1),size(Image_Marker,2));
Image_Inside_Marker(Center_Y_Bubble-Marker_Size:Center_Y_Bubble+Marker_Size,Center_X_Bubble-Marker_Size:Center_X_Bubble+Marker_Size)=1;
imagesc(Image_Inside_Marker)
Image_Inside_Marker_all=Image_Inside_Marker | Image_Marker_BGD;
%% inside watershed
Image_Inside_Mask=(max(Image_Greyed(:))-Image_Greyed);
Image_Inside_Mask=abs(Image_Inside_Mask-mean(mean(Image_Inside_Mask(((Xgrid-xc).^2+(Ygrid-yc).^2).^0.5<R_hole))));
Image_Inside_Mask=Image_Inside_Mask.*Mask_Hole;
Image_Inside_Mask_2=imimposemin(Image_Inside_Mask,Image_Inside_Marker_all);
imagesc(Image_Inside_Mask_2);

Image_Result_Inside=watershed(Image_Inside_Mask_2);

imagesc(Image_Result_Inside);
axis equal


%% Top&bottom hat