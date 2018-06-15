clear all

% Data reading & stacking
frame_width=648;
frame_height=488;

X_overlapping=10;
Y_overlapping=10;


frame_width_eff=frame_width-X_overlapping;
frame_height_eff=frame_height-Y_overlapping;


num_of_frame_per_division=16;

normalization_factor=50;

X_mosaic_number=3;
Y_mosaic_number=3;


division_starting_index=4;
total_division_number=3;

frame_starting_index=21;
total_frame_number=20;




division_number=division_starting_index:(division_starting_index+total_division_number-1);
total_FOV_number=X_mosaic_number*Y_mosaic_number;
k=0;

temp_frame_volume=zeros(frame_width,frame_height,total_division_number*num_of_frame_per_division);

stiched_image=zeros(frame_width_eff*X_mosaic_number+X_overlapping,frame_height_eff*Y_mosaic_number+Y_overlapping);

% correction image generation

% correction 1 (edge summation correction)

correction_A=ones(frame_width,frame_height);

% left&right bound

for tt=1:X_overlapping
    correction_A(tt,:)=correction_A(tt,:)*(tt/(X_overlapping+1));
    correction_A(frame_width-tt+1,:)=correction_A(frame_width-tt+1,:)*(tt/(X_overlapping+1));
end
for tt=1:Y_overlapping
    correction_A(:,tt)=correction_A(:,tt)*(tt/(Y_overlapping+1));
    correction_A(:,frame_height-tt+1)=correction_A(:,frame_height-tt+1)*(tt/(Y_overlapping+1));
end

imagesc(correction_A);
colormap('gray');
axis equal
xlim([0 size(correction_A,2)]);
ylim([0 size(correction_A,1)]);

Gaussian_X_width=400;
Gaussian_Y_width=300;
X_offset=-100;
Y_offset=40;

correction_B_X=ones(frame_width,frame_height);
correction_B_Y=ones(frame_width,frame_height);

for tt=1:frame_height
    correction_B_X(:,tt)=gaussmf((1:frame_width),[Gaussian_X_width frame_width/2+X_offset]);
end
for tt=1:frame_width
    correction_B_Y(tt,:)=gaussmf((1:frame_height),[Gaussian_X_width frame_height/2+Y_offset]);
end
correction_B=1./(correction_B_X.*correction_B_Y);

correction_image=correction_A.*correction_B;
%correction_image(:)=1;
imagesc(correction_B);
colormap('gray');
axis equal
xlim([0 size(correction_B,2)]);
ylim([0 size(correction_B,1)]);

for N=0:(total_FOV_number-1)
    
    X_FOV_number=rem(N,X_mosaic_number); %0~2
    Y_FOV_number=2-floor(N/X_mosaic_number); %0~2

    
    folder_path=sprintf('D:\\OCT data\\150629\\2015_0629_5th_%d\\',N+1);
    cd(folder_path);
    mkdir('divide');


    for p=1:length(division_number)
        file_path=[folder_path sprintf('%08d',division_number(p))];
        fin=fopen(file_path);
        A=fread(fin,[frame_width,frame_height*num_of_frame_per_division],'float32','b');
        if fin ==-1
             k=k+1;
             fclose('all');
        else
            for q=1:num_of_frame_per_division
                temp_frame_volume(:,:,q+((p-1))*num_of_frame_per_division)=A(:,(frame_height*(q-1)+1):frame_height*q);
            end
             fclose('all');
        end
    end
    %temp_frame_volume(:)=1;
    averaged_image=mean(temp_frame_volume(:,:,frame_starting_index:(frame_starting_index+total_frame_number-1)),3).*correction_image;
    
    stiched_image(((X_FOV_number)*frame_width_eff+1):((X_FOV_number)*frame_width_eff+frame_width),((Y_FOV_number)*frame_height_eff+1):((Y_FOV_number)*frame_height_eff+frame_height))=stiched_image(((X_FOV_number)*frame_width_eff+1):((X_FOV_number)*frame_width_eff+frame_width),((Y_FOV_number)*frame_height_eff+1):((Y_FOV_number)*frame_height_eff+frame_height))+averaged_image;

    
end
cmin=8;
cmax=50;


imagesc(stiched_image);
colormap('gray');
caxis([cmin cmax]);
axis equal
xlim([0 size(stiched_image,2)]);
ylim([0 size(stiched_image,1)]);

Normailzed_image=(stiched_image-cmin)/(cmax);
Normailzed_image(Normailzed_image>1)=0;

% Stacking

imwrite(Normailzed_image,[cd,'\divide\','averaged_frame','.png']);
