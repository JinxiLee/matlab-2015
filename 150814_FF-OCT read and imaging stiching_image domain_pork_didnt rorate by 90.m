clear all
tic
% Data reading & stacking
folder_path_without_index='D:\OCT data\150720\2015_0720_150720_pork_5th (hight p 80 micron)';

%parts=regexp(folder_path_without_index, '\','split');
%Parent_folder=folder_path_without_index(1:(length(folder_path_without_index)-length(parts{end})));

frame_width=648;
frame_height=488;

X_overlapping=25;
Y_overlapping=21;

cmin=10;
cmax_set=16;%20;


frame_width_eff=frame_width-X_overlapping;
frame_height_eff=frame_height-Y_overlapping;


num_of_frame_per_division=16;

normalization_factor=50;

X_mosaic_starting_index=2;  %start from 0
Y_mosaic_starting_index=7;


X_mosaic_number=3;%12;%12;
Y_mosaic_number=4;%16;%16;


starting_frame_index=1;   %0.2 micron per frame, start from 1
frame_average_factor=20;    %averaged to stack
total_stack_number=1;

%

total_FOV_number=X_mosaic_number*Y_mosaic_number;
k=0;


%stiched_volume=zeros(frame_width_eff*X_mosaic_number+X_overlapping,frame_height_eff*Y_mosaic_number+Y_overlapping,num_of_frame_per_division);
% correction image generation

% correction 1 (edge summation correction)

correction_A=ones(frame_width,frame_height);

% left&right bound

for tt=1:X_overlapping
    correction_A(tt,:)=correction_A(tt,:)*((tt-1)/((X_overlapping-1)));
    correction_A(frame_width-tt+1,:)=correction_A(frame_width-tt+1,:)*((tt-1)/((X_overlapping-1)));
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

G_Ratio=0.65;
Gaussian_X_width=400*G_Ratio;
Gaussian_Y_width=300*G_Ratio;
X_offset=70;
Y_offset=0;

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

dlmwrite('correction_B.txt',correction_B,'delimiter','\t','newline','pc');

%correction_image(:)=1;
imagesc(correction_B);
colormap('gray');
axis equal
xlim([0 size(correction_B,2)]);
ylim([0 size(correction_B,1)]);

stiched_images=zeros(frame_width_eff*X_mosaic_number+X_overlapping,frame_height_eff*Y_mosaic_number+Y_overlapping,total_stack_number);
for NNN=1:total_stack_number
    
    division_starting_index=floor((starting_frame_index+(NNN-1)*frame_average_factor)/num_of_frame_per_division);   %start from zero
    the_starting_frame_index_in_the_first_division=starting_frame_index+(NNN-1)*frame_average_factor-division_starting_index*num_of_frame_per_division;
    division_end_index=((starting_frame_index+NNN*frame_average_factor)/num_of_frame_per_division);
    division_number=division_starting_index:division_end_index;
    temp_frame_volume=zeros(frame_width,frame_height,num_of_frame_per_division*length(division_number));
    stiched_image=zeros(frame_width_eff*X_mosaic_number+X_overlapping,frame_height_eff*Y_mosaic_number+Y_overlapping);
    
    for N=0:(total_FOV_number-1)
        X_number=rem(N,X_mosaic_number)+X_mosaic_starting_index; %0~2
        Y_number=floor(N/X_mosaic_number)+Y_mosaic_starting_index; %0~2
        X_FOV_number=X_number-X_mosaic_starting_index;
        Y_FOV_number=Y_mosaic_number-1-Y_number+Y_mosaic_starting_index;
        if (X_number<10)&&(Y_number<10)
            folder_path=sprintf('%s_% d_% d\\',folder_path_without_index,X_number,Y_number);
        elseif (X_number>9)&&(Y_number<10)
            folder_path=sprintf('%s_%d_% d\\',folder_path_without_index,X_number,Y_number);
        elseif (Y_number>9)&&(X_number<10)
            folder_path=sprintf('%s_% d_%d\\',folder_path_without_index,X_number,Y_number);
        else
            folder_path=sprintf('%s_%d_%d\\',folder_path_without_index,X_number,Y_number);
        end
        cd(folder_path);
        for NN=1:length(division_number)
            file_path=[folder_path sprintf('%08d',division_number(NN))];
            fin=fopen(file_path);
            A=fread(fin,[frame_width,frame_height*num_of_frame_per_division],'float32','b');
            if fin ==-1
                k=k+1;
                fclose('all');
            else
                for q=1:num_of_frame_per_division
                    %temp_frame_volume(:,:,(NN-1)*num_of_frame_per_division+q)=A(:,(frame_height*(q-1)+1):frame_height*q).*correction_image;
                    temp_frame_volume(:,:,(NN-1)*num_of_frame_per_division+q)=A(:,(frame_height*(q-1)+1):frame_height*q);
                end
                fclose('all');
            end
        end
        Averaged_frame=mean(temp_frame_volume(:,:,the_starting_frame_index_in_the_first_division:(the_starting_frame_index_in_the_first_division+frame_average_factor-1)),3);
        Averaged_frame=(Averaged_frame-cmin)/cmax_set;  %因為無從做max intensity判斷, 只好用固定值
        %Averaged_frame(Averaged_frame>1)=1;
        Averaged_frame(Averaged_frame<0)=0; 
        Averaged_frame=Averaged_frame.*correction_image;        
        stiched_image(((X_FOV_number)*frame_width_eff+1):((X_FOV_number)*frame_width_eff+frame_width),((Y_FOV_number)*frame_height_eff+1):((Y_FOV_number)*frame_height_eff+frame_height))=stiched_image(((X_FOV_number)*frame_width_eff+1):((X_FOV_number)*frame_width_eff+frame_width),((Y_FOV_number)*frame_height_eff+1):((Y_FOV_number)*frame_height_eff+frame_height))+Averaged_frame;
        disp(N);

    end
    %if max(max(stiched_image))<cmax_set
    %    cmax=max(max(stiched_image));
    %else
        %cmax=cmax_set;
    %end
    %Normailzed_image=(stiched_image-cmin)/cmax;
    %Normailzed_image(Normailzed_image>1)=1;
    %Normailzed_image(Normailzed_image<0)=0;
    stiched_image_write=stiched_image;
    stiched_image_write(stiched_image_write>1)=1;
    mkdir(sprintf('%s_stiched_image\\',folder_path_without_index));
    imwrite(stiched_image_write,[sprintf('%s_stiched_image\\',folder_path_without_index),sprintf('stiched_image_%.1f micron_X%d_Y%d',(starting_frame_index+(NNN-1)*frame_average_factor)*0.2,X_mosaic_starting_index,Y_mosaic_starting_index),'.png']);
    fout=fopen([sprintf('%s_stiched_image\\',folder_path_without_index),sprintf('stiched_image_%.1f micron_X%d_Y%d',(starting_frame_index+(NNN-1)*frame_average_factor)*0.2,X_mosaic_starting_index,Y_mosaic_starting_index)],'w+');
    fwrite(fout,stiched_image,'float32','b');
end
imagesc(stiched_image_write);
    colormap('gray');
    caxis([0 1]);
    axis equal
    xlim([0 size(stiched_image,2)]);
    ylim([0 size(stiched_image,1)]);
    
    %imagesc(histeq(stiched_image,[0 1]));
    %colormap('gray');
    %caxis([0 1]);
    %axis equal
    %xlim([0 size(stiched_image,2)]);
    %ylim([0 size(stiched_image,1)]);
    
    
    %fin2=fopen([cd,'\divide\',sprintf('stiched_image_%.1f micron',(starting_frame_index+(NNN-1)*frame_average_factor)*0.2)]);
    %QQQ=fread(fin2,[size(stiched_image,1),size(stiched_image,2)],'float32','b');

    %imagesc(histeq(stiched_image,[0.2:0.001:0.5]));

    %for QQQ=1:4
    %    Normailzed_image=(mean(stiched_volume(:,:,((QQQ-1)*4+1):((QQQ)*4)),3)-cmin)/cmax;
    %    Normailzed_image(Normailzed_image>1)=1;
    %    %stiched_images(:,:,NN)=mean(stiched_volume,3);
    %    imwrite(Normailzed_image,[cd,'\divide\',sprintf('stiched_image_%d.png',(NN-1)*4+QQQ),'.png']);
    %end
    %for QQQ=1:16
    %    Normailzed_image=(mean(stiched_volume(:,:,QQQ),3)-cmin)/cmax;
    %    Normailzed_image(Normailzed_image>1)=1;
        %stiched_images(:,:,NN)=mean(stiched_volume,3);
    %    imwrite(Normailzed_image,[cd,'\divide3\',sprintf('stiched_image_%d.png',(NN-1)*16+QQQ),'.png']);
    %end
    %stiched_volume=zeros(frame_width_eff*X_mosaic_number+X_overlapping,frame_height_eff*Y_mosaic_number+Y_overlapping,num_of_frame_per_division);
toc
%
disable=1;
if disable==0
frame_starting_index=85;
total_frame_number=20;

%stiched_image=mean(stiched_volume(:,:,frame_starting_index:(frame_starting_index+total_frame_number-1)),3);


%
NNN=1;
stiched_image=stiched_images(:,:,NNN);

imagesc(stiched_image);
colormap('gray');
caxis([cmin cmax]);
axis equal
xlim([0 size(stiched_image,2)]);
ylim([0 size(stiched_image,1)]);

Normailzed_image=(stiched_image-cmin)/cmax;
Normailzed_image(Normailzed_image>1)=1;

imagesc(Normailzed_image);
colormap('gray');
caxis([0 1]);
axis equal
xlim([0 size(stiched_image,2)]);
ylim([0 size(stiched_image,1)]);

%% Stacking
%Frame_spacing=5;
%Stacking_starting_frame=1;
%Total_stacked_frame_number=38;

%for qq=1:Total_stacked_frame_number
%    findex=Stacking_starting_frame+(qq-1)*Frame_spacing;
%    stiched_image=mean(stiched_volume(:,:,findex:(findex+total_frame_number-1)),3);
%    Normailzed_image=(stiched_image-cmin)/cmax;
%    Normailzed_image(Normailzed_image>1)=1;
%    imwrite(Normailzed_image,[cd,'\divide\',sprintf('averaged_frame_%d.png',qq),'.png']);
%end

end
 fclose all