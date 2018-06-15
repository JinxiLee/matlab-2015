clear all
tic
% Data reading & stacking
folder_path_without_index='D:\OCT data\151204_Set 3\YZ (txt)\';

frame_width=368;
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


number=0:-1:-12;

Image_Volume=zeros(frame_width,frame_height,length(number));

for p=1:length(number)
    Image_Volume(:,:,p)=dlmread(sprintf('%s%d.txt',folder_path_without_index,number(p)))';
end


%%
NN=4;

imagesc(Image_Volume(:,:,NN));
caxis([0 500]);


%% correction array

correction_A=ones(1,frame_height);

% left&right bound

for tt=1:Y_overlapping
    correction_A(:,tt)=correction_A(:,tt)*((tt-1)/(Y_overlapping-1));
    correction_A(:,frame_height-tt+1)=correction_A(:,frame_height-tt+1)*((tt-1)/((Y_overlapping-1)));
end


imagesc(correction_A);
colormap('gray');
axis equal
xlim([0 size(correction_A,2)]);
ylim([0 size(correction_A,1)]);
%%

stiched_image=zeros(frame_width,frame_height_eff*length(number)+Y_overlapping);

for p=1:length(number)
    X_FOV_number=0;
    Y_FOV_number=length(number)-p;
    Averaged_frame=Image_Volume(:,:,p);
    stiched_image(((X_FOV_number)*frame_width_eff+1):((X_FOV_number)*frame_width_eff+frame_width),((Y_FOV_number)*frame_height_eff+1):((Y_FOV_number)*frame_height_eff+frame_height))=stiched_image(((X_FOV_number)*frame_width_eff+1):((X_FOV_number)*frame_width_eff+frame_width),((Y_FOV_number)*frame_height_eff+1):((Y_FOV_number)*frame_height_eff+frame_height))+Averaged_frame;
    disp(p);
end

stiched_image=stiched_image/max(stiched_image(:));
%%

imagesc(stiched_image);
colormap('gray');

axis equal
    
    xlim([0 size(stiched_image,2)]);
    ylim([0 size(stiched_image,1)]);
    axis off
%%
imwrite(stiched_image,[sprintf('%sStiched_image',folder_path_without_index),'.png']);

%% max searching
Slope_array=zeros(length(number),1);
Max_array=zeros(frame_height_eff*length(number)+Y_overlapping,1);
Mean_Max_array=zeros(length(number),1);

for p=1:length(number)
    Y_FOV_number=length(number)-p;
    Y_array=1.33*(1:size(Image_Volume,2));
    [value max_temp]=max(Image_Volume(:,:,p),[],1);
    Z_array=max_temp*0.2;
    CFIT= fit(Y_array',Z_array','poly1');
    Slope_array(p)=CFIT.p1;
    Max_array(((Y_FOV_number)*frame_height_eff+1):((Y_FOV_number)*frame_height_eff+frame_height))=Max_array(((Y_FOV_number)*frame_height_eff+1):((Y_FOV_number)*frame_height_eff+frame_height))+Z_array'.*correction_A';
    Mean_Max_array(p)=mean(Z_array);
end

%% to calculate the residual of Mean_Max_array
CFIT_mean= fit((1:length(Mean_Max_array))',Mean_Max_array,'poly1');
Mean_Max_array_fit=CFIT_mean.p1*(1:length(Mean_Max_array))'+CFIT_mean.p2;

plot(1:length(Mean_Max_array),Mean_Max_array,1:length(Mean_Max_array),Mean_Max_array_fit);
legend('Interference signal position','fitting');
ylabel('Interference signal position (micron)');
xlabel('FOV number');
plot(1:length(Mean_Max_array),Mean_Max_array-Mean_Max_array_fit);
ylabel('Interference signal position with proper alignment (micron)');
xlabel('FOV number');

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


