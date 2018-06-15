clear all
tic
% Data reading & stacking

%parts=regexp(folder_path_without_index, '\','split');
%Parent_folder=folder_path_without_index(1:(length(folder_path_without_index)-length(parts{end})));

frame_width=648;
frame_height=488;

X_overlapping=25;
Y_overlapping=21;

cmin=10;
cmax_set=80;%20;


frame_width_eff=frame_width-X_overlapping;
frame_height_eff=frame_height-Y_overlapping;


num_of_frame_per_division=16;

normalization_factor=50;

X_mosaic_starting_index=0;
Y_mosaic_starting_index=0;


X_mosaic_number=12;%12;%12;
Y_mosaic_number=16;%16;%16;


starting_frame_index=25;   %0.2 micron per frame, start from 1
frame_average_factor=20;    %averaged to stack
total_stack_number=1;

%

total_FOV_number=X_mosaic_number*Y_mosaic_number;
k=0;


%stiched_volume=zeros(frame_width_eff*X_mosaic_number+X_overlapping,frame_height_eff*Y_mosaic_number+Y_overlapping,num_of_frame_per_division);
% correction image generation

% correction 1 (edge summation correction)

correction_A_Horizontal_L=ones(frame_width,frame_height);
correction_A_Horizontal_R=ones(frame_width,frame_height);
correction_A_Vertical_B=ones(frame_width,frame_height);
correction_A_Vertical_T=ones(frame_width,frame_height);
%correction_A=ones(frame_width,frame_height);
% left&right bound

for tt=1:X_overlapping
    correction_A_Horizontal_L(tt,:)=correction_A_Horizontal_L(tt,:)*((tt-1)/((X_overlapping-1)));
    correction_A_Horizontal_R(frame_width-tt+1,:)=correction_A_Horizontal_R(frame_width-tt+1,:)*((tt-1)/((X_overlapping-1)));
    %correction_A(tt,:)=correction_A(tt,:)*((tt-1)/((X_overlapping-1)));
    %correction_A(frame_width-tt+1,:)=correction_A(frame_width-tt+1,:)*((tt-1)/((X_overlapping-1)));
end
for tt=1:Y_overlapping
    correction_A_Vertical_B(:,tt)=correction_A_Vertical_B(:,tt)*((tt-1)/(Y_overlapping-1));
    correction_A_Vertical_T(:,frame_height-tt+1)=correction_A_Vertical_T(:,frame_height-tt+1)*((tt-1)/((Y_overlapping-1)));
    %correction_A(:,tt)=correction_A(:,tt)*((tt-1)/(Y_overlapping-1));
    %correction_A(:,frame_height-tt+1)=correction_A(:,frame_height-tt+1)*((tt-1)/((Y_overlapping-1)));

end
correction_A=correction_A_Horizontal_L.*correction_A_Horizontal_R.*correction_A_Vertical_B.*correction_A_Vertical_T;
correction_A_B=correction_A_Horizontal_L.*correction_A_Horizontal_R.*correction_A_Vertical_T;
correction_A_T=correction_A_Horizontal_L.*correction_A_Horizontal_R.*correction_A_Vertical_B;
correction_A_L=correction_A_Horizontal_R.*correction_A_Vertical_B.*correction_A_Vertical_T;
correction_A_R=correction_A_Horizontal_L.*correction_A_Vertical_B.*correction_A_Vertical_T;

correction_A_BL=correction_A_Horizontal_R.*correction_A_Vertical_T;
correction_A_BR=correction_A_Horizontal_L.*correction_A_Vertical_T;
correction_A_TL=correction_A_Horizontal_R.*correction_A_Vertical_B;
correction_A_TR=correction_A_Horizontal_L.*correction_A_Vertical_B;


imagesc(correction_A);
colormap('gray');
axis equal
xlim([0 size(correction_A,2)]);
ylim([0 size(correction_A,1)]);

G_Ratio=0.7;
Gaussian_X_width=400*G_Ratio;
Gaussian_Y_width=300*G_Ratio;
X_offset=10;
Y_offset=45;

correction_B_X=ones(frame_width,frame_height);
correction_B_Y=ones(frame_width,frame_height);

for tt=1:frame_height
    correction_B_X(:,tt)=gaussmf((1:frame_width),[Gaussian_X_width frame_width/2+X_offset]);
end
for tt=1:frame_width
    correction_B_Y(tt,:)=gaussmf((1:frame_height),[Gaussian_Y_width frame_height/2+Y_offset]);
end
correction_B=1./(correction_B_X.*correction_B_Y);

correction_image=correction_A;%.*correction_B;

dlmwrite('correction_B.txt',correction_B,'delimiter','\t','newline','pc');

%correction_image(:)=1;
imagesc(correction_B);
colormap('gray');
axis equal
xlim([0 size(correction_B,2)]);
ylim([0 size(correction_B,1)]);

    stiched_image=zeros(frame_width_eff*X_mosaic_number+X_overlapping,frame_height_eff*Y_mosaic_number+Y_overlapping,total_stack_number);
    for N=0:(total_FOV_number-1)
        X_number=rem(N,X_mosaic_number)+X_mosaic_starting_index; %0~2
        Y_number=floor(N/X_mosaic_number)+Y_mosaic_starting_index; %0~2
        X_FOV_number=X_mosaic_number-1-X_number+X_mosaic_starting_index;
        Y_FOV_number=Y_mosaic_number-1-Y_number+Y_mosaic_starting_index;
        
        Averaged_frame=ones(frame_width,frame_height);

        %Averaged_frame=Averaged_frame.*correction_B;  
        
        if (X_FOV_number==0) && (Y_FOV_number==0)
            Averaged_frame=Averaged_frame.*correction_A_BL;
        elseif (X_FOV_number==0) && (Y_FOV_number==(Y_mosaic_number-1))
            Averaged_frame=Averaged_frame.*correction_A_TL;
        elseif (X_FOV_number==(X_mosaic_number-1)) && (Y_FOV_number==0)
            Averaged_frame=Averaged_frame.*correction_A_BR;
        elseif (X_FOV_number==(X_mosaic_number-1)) && (Y_FOV_number==(Y_mosaic_number-1))
            Averaged_frame=Averaged_frame.*correction_A_TR;
        elseif (X_FOV_number==0)
            Averaged_frame=Averaged_frame.*correction_A_L;
        elseif (X_FOV_number==(X_mosaic_number-1))
            Averaged_frame=Averaged_frame.*correction_A_R;
        elseif (Y_FOV_number==0)
            Averaged_frame=Averaged_frame.*correction_A_B;
        elseif (Y_FOV_number==(Y_mosaic_number-1))
            Averaged_frame=Averaged_frame.*correction_A_T;
        else
            Averaged_frame=Averaged_frame.*correction_A;
        end

        
        
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
imagesc(stiched_image);
colormap(gray);
