clear all

window=100;

A=rand(200,200);

B_old=rand(200+window-1,200+window-1);

Filter=fspecial('gaussian', window,window/3);    %gaussian mask

B=conv2(B_old,Filter,'valid');

B=(B-min(B(:)))/(max(B(:))-min(B(:)));


std_A=zeros(200,200);
std_B=zeros(200,200);

mean_A=zeros(200,200);
mean_B=zeros(200,200);


Sample_size=10;
for p=1:200-Sample_size+1
    for q=1:200-Sample_size+1
        Sample=A(p:p+Sample_size-1,q:q+Sample_size-1);
        std_A(p,q)=std(Sample(:));       
        mean_A(p,q)=mean(Sample(:));      
    end
end
for p=1:200-Sample_size+1
    for q=1:200-Sample_size+1
        Sample=B(p:p+Sample_size-1,q:q+Sample_size-1);
        std_B(p,q)=std(Sample(:));       

        mean_B(p,q)=mean(Sample(:));      
    end
end

C=zeros(200,200);

subplot(5,2,1)
imagesc(A);
title('A');
axis equal

subplot(5,2,2)
imagesc(B);
title('B');
axis equal

subplot(5,2,3)
imagesc(std_A);
title('std_A');
axis equal
caxis([0 0.5])

subplot(5,2,4)
imagesc(std_B);
title('std_B');
axis equal
caxis([0 0.5])

subplot(5,2,5)
hist(std_A,0:0.01:0.5);
title('Histogram of std_A');

subplot(5,2,6)
hist(std_B,0:0.01:0.5);
title('Histogram of std_B');

subplot(5,2,7)
imagesc(mean_A);
title('mean_A');
axis equal
caxis([0 1])

subplot(5,2,8)
imagesc(mean_B);
title('mean_B');
axis equal
caxis([0 1])

subplot(5,2,9)
hist(mean_A,0:0.01:1);
title('Histogram of mean_A');

subplot(5,2,10)
hist(mean_B,0:0.01:1);
title('Histogram of mean_B');


%% C generation

New_A=rand(200,200);


New_B_old=rand(200+window-1,200+window-1);

Filter=fspecial('gaussian', window,window/3);    %gaussian mask

New_B=conv2(New_B_old,Filter,'valid');

New_B=(New_B-min(New_B(:)))/(max(New_B(:))-min(New_B(:)));

mask=zeros(200,200);

center_1=[55 55];

center_2=[145 145];

grid_X=repmat([1:200]',[1 200]);
grid_Y=repmat([1:200],[200 1]); 

mask=double((((grid_X-60).^2+(grid_Y-60).^2).^0.5<50)|((grid_X-140).^2+(grid_Y-140).^2).^0.5<50);
imagesc(mask);

C=mask.*New_A+(1-mask).*New_B;


imagesc(C);
colormap('gray');
axis equal

%% D generation

New_A=rand(200,200);


New_B_old=rand(200+window-1,200+window-1);

Filter=fspecial('gaussian', window,window/3);    %gaussian mask

New_B=conv2(New_B_old,Filter,'valid');

New_B=(New_B-min(New_B(:)))/(max(New_B(:))-min(New_B(:)));

mask=zeros(200,200);

center_1=[60 60];

center_2=[140 140];

grid_X=repmat([1:200]',[1 200]);
grid_Y=repmat([1:200],[200 1]); 

mask=double((((grid_X-100).^2+(grid_Y-100).^2).^0.5<80)&((grid_X-100).^2+(grid_Y-100).^2).^0.5>30);
imagesc(mask);

D=mask.*New_A+(1-mask).*New_B;

imagesc(D);
colormap('gray');
axis equal
%%
cd('D:\');
imwrite(A,'A.png','png');
imwrite(B,'B.png','png');
imwrite(C,'C.png','png');
imwrite(D,'D.png','png');