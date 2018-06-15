clear all

%x_size=7501;
x_size=7446;
%y_size=7493;
y_size=7508;
cd('D:\OCT data\150720\2015_0720_150720_pork_2nd(low p 50 micron)_11_15\divide');
file_list=cellstr(['stiched_image_2.0 micron']);%; 'Nevi_ 50'; 'Nevi_100'; 'Navi_140']);
file_name=char(file_list{1});
fin=fopen([cd,'\',file_name]);
ImageToClassify=fread(fin,[x_size,y_size],'float32','b');
ImageForTraining=ImageToClassify;

imagesc(ImageToClassify);
caxis([10 30]);
colormap('gray');
axis equal
xlim([0 size(ImageToClassify,2)]);
ylim([0 size(ImageToClassify,1)]);
%% Sample Define

Sample_Size=249;  %size of the matrix, must be odd number

%Weight_mask=ones(Sample_Size);
Weight_mask=fspecial('gaussian', Sample_Size,Sample_Size/3);    %gaussian mask

%% Sampling Parameter

Downsample_Factor_Sampling=50;

X_Fixed_Sampling=repmat([Sample_Size:Downsample_Factor_Sampling:(size(ImageToClassify,1)-Sample_Size)]',[1 length(Sample_Size:Downsample_Factor_Sampling:(size(ImageToClassify,2))-Sample_Size)]);
Y_Fixed_Sampling=repmat([Sample_Size:Downsample_Factor_Sampling:(size(ImageToClassify,2)-Sample_Size)],[length(Sample_Size:Downsample_Factor_Sampling:(size(ImageToClassify,1))-Sample_Size) 1]);
X_Fixed_Sampling_Array=Sample_Size:Downsample_Factor_Sampling:(size(ImageForTraining,1)-Sample_Size);
Y_Fixed_Sampling_Array=Sample_Size:Downsample_Factor_Sampling:(size(ImageForTraining,2)-Sample_Size);


    
%% Training

Number_of_Class=3;
Number_of_feature=3;
Number_of_Cycle_for_each_Class=1;

Feature_Map=zeros(size(X_Fixed_Sampling,1),size(X_Fixed_Sampling,2),Number_of_feature);

%Centers_X=[3709 5281 664]';
%Centers_X=[5957 1357 4357 2157 530 2730]';     %5957
Centers_X=[5957  4457  530 ]';
%Centers_Y=[3318 4821 5514]';
%Centers_Y=[852 1352 5452 6452 4850 1050]';  %852
Centers_Y=[852  5552  4850 ]';
%Radiuses=[450 300 550]';
%Radiuses=[800 800 1300 800 250 150]';   %800
Radiuses=[200  200  200 ]';
Downsample_Factor_Training=10;     %only every Downsample_Factor point take one point

imagesc(ImageForTraining);
caxis([10 30]);
colormap('gray');
axis equal
xlim([0 size(ImageForTraining,2)]);
ylim([0 size(ImageForTraining,1)]);

hold on

viscircles([Centers_Y,Centers_X],Radiuses,'EdgeColor','b');     %X和Y在顯示上會是反的 所以只好改顯示方向!

hold off

%% Generating training sample list
Training_Sample_Set=cell(Number_of_Class,2);    %2 for x and y coordinate
X_Fixed_Training=repmat([Sample_Size:Downsample_Factor_Training:(size(ImageForTraining,1)-Sample_Size)]',[1 length(Sample_Size:Downsample_Factor_Training:(size(ImageForTraining,2))-Sample_Size)]);
Y_Fixed_Training=repmat([Sample_Size:Downsample_Factor_Training:(size(ImageForTraining,2)-Sample_Size)],[length(Sample_Size:Downsample_Factor_Training:(size(ImageForTraining,1))-Sample_Size) 1]); 

for p=1:Number_of_Class
    X_mask=X_Fixed_Training;
    X=X_Fixed_Training;
    Y=Y_Fixed_Training;
    for q=1:Number_of_Cycle_for_each_Class
        X_mask(((X_Fixed_Training-Centers_X((p-1)*Number_of_Cycle_for_each_Class+q)).^2+(Y_Fixed_Training-Centers_Y((p-1)*Number_of_Cycle_for_each_Class+q)).^2).^0.5<Radiuses((p-1)*Number_of_Cycle_for_each_Class+q))=0;
    end
    %Y(((X_Fixed-Centers_X(p)).^2+(Y_Fixed-Centers_Y(p)).^2).^0.5>Radiuses(p))=0;
    Y(X_mask~=0)=[];
    X(X_mask~=0)=[];
    Training_Sample_Set{p,1}=X(:);
    Training_Sample_Set{p,2}=Y(:);
end
    
%% Extract features in trianing set
Features=cell(Number_of_Class,Number_of_feature);   %we use 3 different features (weighted mean, std, weighted absolute gradient)
for p=1:Number_of_Class
    for u=1:Number_of_feature
        Features{p,u}=zeros([length(Training_Sample_Set{p,1}) 1]);
    end

    for q=1:length(Training_Sample_Set{p,1})
        Temp=ImageForTraining((Training_Sample_Set{p,1}(q)-((Sample_Size-1)/2)):(Training_Sample_Set{p,1}(q)+((Sample_Size-1)/2)),(Training_Sample_Set{p,2}(q)-((Sample_Size-1)/2)):(Training_Sample_Set{p,2}(q)+((Sample_Size-1)/2)));
    %Feature: weighted mean (注意這裡是要用sum)
        Features{p,1}(q)=sum(Temp(:).*Weight_mask(:));
    %Feature: wighted std (normalized)
        Features{p,2}(q)=(sum(((Temp(:)-Features{p,1}(q)).^2).*Weight_mask(:)).^0.5)/Features{p,1}(q);
    %Feature: Kurtosis
        Features{p,3}(q)=sum(((Temp(:)-Features{p,1}(q)).^4).*Weight_mask(:))./Features{p,2}(q)^2;

    %Feature: weighted absolute gradient
        %grad_X=abs(diff(Temp,1,1));
        %grad_X(2:(end+1),:)=grad_X;
        %grad_Y=abs(diff(Temp,1,2));
        %grad_Y(:,2:(end+1))=grad_Y;
        %Features{p,3}(q)=0.5*sum(grad_X(:).*Weight_mask(:))+0.5*sum(grad_Y(:).*Weight_mask(:));
    end
end

%% To calculate the neccessary parameters (mean, covaraince matrix)
Mean_of_Features=cell(Number_of_Class);   %考慮到後面, 這裡用array structure比較好
Covariance_Matrix=cell(Number_of_Class);    %每個cell裡有一個3*3 matrix
for p=1:Number_of_Class
    Mean_of_Features{q}=zeros(Number_of_feature,1);
    Covariance_Matrix{p}=zeros(Number_of_feature);
end
for p=1:Number_of_Class
    for q=1:Number_of_feature
        Mean_of_Features{p}(q)=mean(Features{p,q});
    end
    for u=1:Number_of_feature
        for v=1:Number_of_feature
            Covariance_Matrix{p}(u,v)=sum((Features{p,u}-Mean_of_Features{p}(u)).*(Features{p,v}-Mean_of_Features{p}(v)))/(length(Features{p,u}-1));  %先不用預設function
        end
    end
end

%% Extract features in the large image
                 %這裡用3D matrix看看, 意思有點像3個channel的color image
for p=1:size(Feature_Map,1)
    for q=1:size(Feature_Map,2)
        Temp=ImageToClassify((X_Fixed_Sampling_Array(p)-((Sample_Size-1)/2)):(X_Fixed_Sampling_Array(p)+((Sample_Size-1)/2)),(Y_Fixed_Sampling_Array(q)-((Sample_Size-1)/2)):(Y_Fixed_Sampling_Array(q)+((Sample_Size-1)/2)));
        %Feature: weighted mean (注意這裡是要用sum)
        Feature_Map(p,q,1)=sum(Temp(:).*Weight_mask(:));
        %Feature: wighted std
        Feature_Map(p,q,2)=(sum(((Temp(:)-Feature_Map(p,q,1)).^2).*Weight_mask(:)).^0.5)/Feature_Map(p,q,1);
        %Feature: Kurtosis
        Feature_Map(p,q,3)=sum(((Temp(:)-Feature_Map(p,q,1)).^4).*Weight_mask(:))./Feature_Map(p,q,2)^2;

        %Feature: weighted absolute gradient
        %grad_X=abs(diff(Temp,1,1));
        %grad_X(2:(end+1),:)=grad_X;
        %grad_Y=abs(diff(Temp,1,2));
        %grad_Y(:,2:(end+1))=grad_Y;
        %Feature_Map(p,q,3)=0.5*sum(grad_X(:).*Weight_mask(:))+0.5*sum(grad_Y(:).*Weight_mask(:));
        
    end
    disp(p);
end

%% Classification
QDA= @ (Cov_Matrix,FeatureValueArray,MeanFeatureValueArray) -0.5*log(det(Cov_Matrix))-0.5*(FeatureValueArray-MeanFeatureValueArray')'*inv(Cov_Matrix)*(FeatureValueArray-MeanFeatureValueArray'); 
QDA_classifier=zeros(Number_of_Class,1);
QDA_map=zeros(size(Feature_Map,1),size(Feature_Map,2));
QDA_classifier_3D=zeros(size(Feature_Map,1),size(Feature_Map,2),Number_of_Class);
for p=1:size(Feature_Map,1)
    for q=1:size(Feature_Map,2)
        for u=1:Number_of_Class
            Temp_Feature_Array=Feature_Map(p,q,:);
            Temp_Feature_Array=Temp_Feature_Array(:);
            QDA_classifier(u)=QDA(Covariance_Matrix{u},Temp_Feature_Array,Mean_of_Features{u});
            QDA_classifier_3D(p,q,u)=QDA(Covariance_Matrix{u},Temp_Feature_Array,Mean_of_Features{u});
        end
        [maxvalue maxindex]=max(QDA_classifier);
        QDA_map(p,q)=maxindex;
    end
    disp(p);
end
%mkdir('Adipose');
    %imwrite(Normailzed_Image,[cd,'\Adipose\',file_name,'_normalized_small FOV_brightness reduced_0.15.png'],'png');

imagesc(QDA_map);
caxis([1 3]);
colormap('gray');
axis equal
xlim([0 size(ImageForTraining,2)/Downsample_Factor_Sampling]);
ylim([0 size(ImageForTraining,1)/Downsample_Factor_Sampling]);

hold on

viscircles([Centers_Y/Downsample_Factor_Sampling,Centers_X/Downsample_Factor_Sampling],Radiuses/Downsample_Factor_Sampling,'EdgeColor','b');

hold off
%%

imagesc(QDA_classifier_3D(:,:,1));
caxis([-100 0]);


%% Generate color map
Offset=100;
Color_map=QDA_classifier_3D+Offset;
Sum_map=zeros(size(Color_map,1),size(Color_map,2),3);
for p=1:3
    Sum_map(:,:,p)=sum(QDA_classifier_3D,3);
end
Color_map=Color_map./Sum_map;

%% Generate large QDA map


X_large=repmat([1:(size(ImageToClassify,1))]',[1 length(1:(size(ImageToClassify,2)))]);
Y_large=repmat([1:(size(ImageToClassify,2))],[length(1:(size(ImageToClassify,1))) 1]);
QDA_classifier_3D_Large=zeros(size(ImageToClassify,1),size(ImageToClassify,2),Number_of_Class);

for p=1:Number_of_Class
    QDA_classifier_3D_Large(:,:,p)=interp2(Y_Fixed_Sampling,X_Fixed_Sampling,QDA_classifier_3D(:,:,p),Y_large,X_large,'cubic',1);
end
