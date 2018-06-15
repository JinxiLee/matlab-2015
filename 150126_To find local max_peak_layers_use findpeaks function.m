clear all

Frame_index=1:240;
unit_length=0.19;
Number_of_Layers_Min=5;
Threshold=30;
Distance_min=3;
ROI=[400 560];

cd('D:\Users\TuanShu\131126_skin');
Volume(1:length(Frame_index),1:ROI(1),1:ROI(2))=0;
for p=1:length(Frame_index)
    Image_Temp=imread(sprintf('%08i.tif',Frame_index(p)));
    Volume(p,:,:)=Image_Temp((size(Image_Temp,1)/2-ROI(1)/2):(size(Image_Temp,1)/2+ROI(1)/2)-1,(size(Image_Temp,2)/2-ROI(2)/2):(size(Image_Temp,2)/2+ROI(2)/2)-1,1);
    disp(p);
end
%%

imagesc(Volume(:,:,N)-./);
%%
Separation_Index_1=470;
Separation_Index_2=520;


for w=1:size(Volume,3)
    Frame=Volume(:,:,w);
    for p=1:size(Frame,2)    
        Signal=Frame(:,p);
        [max_temp max_index_temp]=findpeaks(Signal,'MINPEAKDISTANC',Distance_min,'MINPEAKHEIGHT',Threshold); %findpeaks如果在那條線找不到peak會跳error, 應該沒關係
        stored_max(1:length(max_index_temp),p,w)=max_index_temp;
    end
    disp(w);
end
stored_max(stored_max==0)=1;


%%

N=262;
Frame_display(:,:)=Volume(:,N,:);
imagesc(Frame_display,'ydata',[0:unit_length:size(Frame_display,1)*unit_length]);
hold on
stored_max_2d(:,:)=stored_max(:,N,:);
plot(stored_max_2d'*unit_length,'.');
hold off

%%
clear Number_of_Layers Total_Thickness
stored_max_Position=stored_max*unit_length;
stored_max_Position(stored_max_Position==unit_length)=0;
Thickness_Volume=diff(stored_max_Position,[],1);
Thickness_Volume(Thickness_Volume<0)=0;
Thickness_Volume(size(Thickness_Volume,1)+1,:,:)=0;
%算總厚度與layer數的2D map
for w=1:size(stored_max_Position,3)
    for p=1:size(stored_max_Position,2)
        linetemp=stored_max_Position(:,p,w);
        max_position_temp=linetemp(linetemp~=0);
        if length(max_position_temp)>Number_of_Layers_Min
            Number_of_Layers(p,w)=length(max_position_temp)-1;
            Total_Thickness(p,w)=max(max_position_temp)-min(max_position_temp);
        else            
            Number_of_Layers(p,w)=-1;
            Total_Thickness(p,w)=-1;
        end
    end
    disp(w);
end

%%
p_range=[1 400];
Number_of_Layers_Patial=Number_of_Layers(p_range(1):p_range(2),:);
Total_Thickness_Patial=Total_Thickness(p_range(1):p_range(2),:);
Thickness_Volume_Patial=Thickness_Volume(:,p_range(1):p_range(2),:);

Number_of_Layers_ALL=Number_of_Layers_Patial(:);
Total_Thickness_ALL=Total_Thickness_Patial(:);
Thickness_Volume_ALL=Thickness_Volume_Patial(:);

Number_of_Layers_ALL=Number_of_Layers_ALL(Number_of_Layers_ALL~=-1);
Total_Thickness_ALL=Total_Thickness_ALL(Total_Thickness_ALL~=-1);
Thickness_Volume_ALL=Thickness_Volume_ALL(Thickness_Volume_ALL~=0);

%%

Number_of_Layers_mean=mean(Number_of_Layers_ALL(Number_of_Layers_ALL~=-1));
Number_of_Layers_SD=std(Number_of_Layers_ALL(Number_of_Layers_ALL~=-1));
Number_of_Layers_max=max(Number_of_Layers_ALL(Number_of_Layers_ALL~=-1));
Number_of_Layers_min=min(Number_of_Layers_ALL(Number_of_Layers_ALL~=-1));

Total_Thickness_mean=mean(Total_Thickness_ALL(Total_Thickness_ALL~=-1));
Total_Thickness_SD=std(Total_Thickness_ALL(Total_Thickness_ALL~=-1));
Total_Thickness_max=max(Total_Thickness_ALL(Total_Thickness_ALL~=-1));
Total_Thickness_min=min(Total_Thickness_ALL(Total_Thickness_ALL~=-1));


Thickness_Volume_mean=mean(Thickness_Volume_ALL);
Thickness_Volume_SD=std(Thickness_Volume_ALL);
Thickness_Volume_max=max(Thickness_Volume_ALL);
Thickness_Volume_min=min(Thickness_Volume_ALL);
%%

fprintf('Number of Layers (mean): %f \n',Number_of_Layers_mean);
fprintf('Number of Layers (SD):   %f \n',Number_of_Layers_SD);
fprintf('Number of Layers (max):  %f \n',Number_of_Layers_max);
fprintf('Total Thickness (mean):  %f micron\n',Total_Thickness_mean);
fprintf('Total Thickness (SD):    %f micron\n',Total_Thickness_SD);
fprintf('Total Thickness (max):   %f micron\n',Total_Thickness_max);
fprintf('Total Thickness (min):   %f micron\n',Total_Thickness_min);


fprintf('Layer Thickness (mean):  %f micron\n',Thickness_Volume_mean);
fprintf('Layer Thickness (SD):    %f micron\n',Thickness_Volume_SD);
fprintf('Layer Thickness (max):   %f micron\n',Thickness_Volume_max);
fprintf('Layer Thickness (min):   %f micron\n',Thickness_Volume_min);

Output=[Number_of_Layers_mean Number_of_Layers_SD Number_of_Layers_max Total_Thickness_mean Total_Thickness_SD Total_Thickness_max Total_Thickness_min Thickness_Volume_mean Thickness_Volume_SD Thickness_Volume_max Thickness_Volume_min]';

%Number of Layers
%Ave Thickness
%Min Thickness
%Max Thikcness
%SD of Thickness