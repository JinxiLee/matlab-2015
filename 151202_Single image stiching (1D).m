clear all
tic
% Data reading & stacking
folder_path_without_index='D:\AMO\151202_About Pitch and Yaw\';

frame_width=648;
frame_height=488;

X_overlapping=25;
Y_overlapping=21;

cmin=12;%12;%11
cmax_set=16;%16;%60;%120;%20;
start_index_offset_from_max=-30;
G_Ratio=0.95;
Gaussian_X_width=400*G_Ratio;
Gaussian_Y_width=300*G_Ratio;
X_offset=10;
Y_offset=30;
if_notch=0;


frame_width_eff=frame_width-X_overlapping;
frame_height_eff=frame_height-Y_overlapping;


%parts=regexp(folder_path_without_index, '\','split');
%Parent_folder=folder_path_without_index(1:(length(folder_path_without_index)-length(parts{end})));
correction_set=3;
map_method='mean';  %mean std kurtosis skewness diff

Image_Volume=zeros(648,488,13);


for p=1:13
    Image_Volume(:,:,p)=imread(sprintf('%s%d.jpg',folder_path_without_index,p),'jpg')';
end


%%
NN=13;

imagesc(Image_Volume(:,:,NN));


%%

correction_A=ones(frame_width,frame_height);

% left&right bound

for tt=1:X_overlapping
    correction_A(tt,:)=correction_A(tt,:)*((tt-1)/((X_overlapping-1)));
    correction_A(frame_width-tt+1,:)=correction_A(frame_width-tt+1,:)*((tt-1-0.5)/((X_overlapping-1)));
end
for tt=1:Y_overlapping
    correction_A(:,tt)=correction_A(:,tt)*((tt-1)/(Y_overlapping-1));
    correction_A(:,frame_height-tt+1)=correction_A(:,frame_height-tt+1)*((tt-1)/((Y_overlapping-1)));
end


imagesc(correction_A);
colormap('gray');
axis equal
xlim([0 size(correction_A,2)]);
ylim([0 size(correction_A,1)]);


correction_B_X=ones(frame_width,frame_height);
correction_B_Y=ones(frame_width,frame_height);

for tt=1:frame_height
    correction_B_X(:,tt)=gaussmf((1:frame_width),[Gaussian_X_width frame_width/2+X_offset]);
end
for tt=1:frame_width
    correction_B_Y(tt,:)=gaussmf((1:frame_height),[Gaussian_Y_width frame_height/2+Y_offset]);
end
correction_B=1./(correction_B_X.*correction_B_Y);

correction_image=correction_A.*correction_B;


%correction_image(:)=1;
imagesc(correction_image);
colormap('gray');
axis equal
xlim([0 size(correction_B,2)]);
ylim([0 size(correction_B,1)]);
%%

stiched_image=zeros(frame_width_eff*13+X_overlapping,frame_height_eff*1+Y_overlapping);

for p=1:13
    X_FOV_number=13-p;
    Y_FOV_number=0;
    Averaged_frame=Image_Volume(:,:,p).*correction_A;
    stiched_image(((X_FOV_number)*frame_width_eff+1):((X_FOV_number)*frame_width_eff+frame_width),((Y_FOV_number)*frame_height_eff+1):((Y_FOV_number)*frame_height_eff+frame_height))=stiched_image(((X_FOV_number)*frame_width_eff+1):((X_FOV_number)*frame_width_eff+frame_width),((Y_FOV_number)*frame_height_eff+1):((Y_FOV_number)*frame_height_eff+frame_height))+Averaged_frame;
    disp(p);
end

stiched_image=stiched_image/max(stiched_image(:));
%%

imagesc(stiched_image);

axis equal
    
%%
imwrite(stiched_image,[sprintf('%sStiched_image',folder_path_without_index),'.png']);


%% Binary

Xgrid(1:size(stiched_image,1),1:size(stiched_image,2))=0;
Ygrid(1:size(stiched_image,1),1:size(stiched_image,2))=0;

for p=1:size(stiched_image,1)
    Xgrid(p,:)=p;
end
for q=1:size(stiched_image,2)
    Ygrid(:,q)=q;
end
%%

level=0.2;
pixel_TH=5000;
Ball_Size=1;
Element=fspecial('disk', Ball_Size*3);
Element(Element>0)=1;
Image_Closed=imclose(stiched_image,Element);
Image_Edge=double(edge(Image_Closed,'sobel'));

Image_Edge(Ygrid>300)=0;

imagesc(Image_Edge);

%%
Image_Edge_inv=Image_Edge(:,size(Image_Edge,2):-1:1);

%%
Y_Profile=zeros(size(Image_Edge,1),1);
for p=1:size(Image_Edge,1)
    Q=find(Image_Edge_inv(p,:)>0,1,'first');
    if isempty(Q)~=1
        Y_Profile(p)=Q;
    end
end
plot(Y_Profile);

%%

X_Start=651;
X_End=7450;

Y_Profile_Patial=Y_Profile(X_Start:X_End);

FOV_N=(1:length(Y_Profile_Patial))/688;
Micron_N=(1:length(Y_Profile_Patial))/5;
plot(FOV_N,Y_Profile_Patial);
xlabel('FOV (0.84 mm)');
ylabel('Edge Position (micron)');


