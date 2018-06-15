clear all

x_size=7501;
y_size=7493;
cd('D:\\OCT data\\Raw data');
file_list=cellstr(['Nevi_ 25'; 'Nevi_ 50'; 'Nevi_100'; 'Navi_140']);
%file_list=cellstr(['Dermis_10'; 'Dermis_25'; 'Dermis_50']); ; 
%file_list=cellstr(['Adipose_10'; 'Adipose_15'; 'Adipose_20'; 'Adipose_25'; 'Adipose_50']);
%file_list=cellstr(['150715_1st  ';'150715_2nd  ';'150715_3rd  ';'150715_4th  ';'150715_5th  ';'150715_6th  ';'150716_1st  ';'150716_2nd_1';'150716_2nd_2']);

x_Start=1;
y_Start=1;

Size=7493;


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
    cmin=0;
    cmax=0.15;
    
    
    Normailzed_Image=(TheSmallImage-cmin)/cmax;
    Normailzed_Image(Normailzed_Image>1)=1;
    Normailzed_Image(Normailzed_Image<0)=0;

    imagesc(Normailzed_Image);
    caxis([0 2]);
    colormap('gray');
    axis equal
    xlim([0 size(Normailzed_Image,2)]);
    ylim([0 size(Normailzed_Image,1)]);
    %%
    mkdir('Adipose');
    imwrite(Normailzed_Image,[cd,'\Adipose\',file_name,'_normalized_small FOV_brightness reduced_0.15.png'],'png');
end
