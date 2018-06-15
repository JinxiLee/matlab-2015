clear all

x_size=7446;
y_size=7508;
cd('D:\OCT data\Raw data\');
file_list=cellstr(['150715_4th']);
%file_list=cellstr(['150715_1st  ';'150715_2nd  ';'150715_3rd  ';'150715_4th  ';'150715_5th  ';'150715_6th  ';'150716_1st  ';'150716_2nd_1';'150716_2nd_2']);
cmin_std=0.05;
cmax_std=0.1;

x_Start=6000;
y_Start=6000;

Size=1200;


for p=1:length(file_list)

    file_name=char(file_list{p});
    fin=fopen([cd,'\',file_name]);
    TheLargeImage=fread(fin,[x_size,y_size],'float32','b');
    TheSmallImage=TheLargeImage(x_Start:(x_Start+Size-1),y_Start:(y_Start+Size-1));

    %% Normalize first assume min is zero

    %TheLargeImage=TheLargeImage/max(TheLargeImage(:));

    %%

    %Histogram_TheLargeImage=hist(TheLargeImage(:),1000);
    %hist(TheLargeImage(:),1000);

    %%
    cmin=10*mean(mean(TheSmallImage))/14.3428;
    cmax=20*mean(mean(TheSmallImage))/14.3428;
    
    
    Normailzed_Image=(TheSmallImage-cmin)/cmax;
    Normailzed_Image(Normailzed_Image>1)=1;
    Normailzed_Image(Normailzed_Image<0)=0;

    imagesc(Normailzed_Image);
    caxis([0 1]);
    colormap('gray');
    axis equal
    xlim([0 size(Normailzed_Image,2)]);
    ylim([0 size(Normailzed_Image,1)]);
    %%
    mkdir('dynamic_normalized_small FOV10');
    imwrite(Normailzed_Image,[cd,'\dynamic_normalized_small FOV10\',file_name,'_normalized_small FOV.png'],'png');
end
