clear all

%path='D:\Users\TuanShu\140626_RBC\1002.bin';
folder_path=sprintf('D:\\OCT data\\2015_0610_1251_Mirau 6_glass_it 200\\');
cd(folder_path);
mkdir('divide');

N_frame=19;
width=648;
height=488;
num_of_division=16;
k=0;
normalization_factor=50;

for slice_num=1:length(N_frame)
    file_path=[folder_path sprintf('%08d',N_frame(slice_num))];
    fin=fopen(file_path);
    A=fread(fin,[width,height*num_of_division],'float32','b');
    if fin ==-1
         k=k+1;
         fclose('all');
    else
        for i=1:num_of_division
            slice=A(:,(height*(i-1)+1):height*i);
            %     dlmwrite([cd,'\divide\',sprintf('%d',(slice_mun-1)*10+i),'.txt'],slice);
            imwrite(slice/normalization_factor,[cd,'\divide\',sprintf('%d',(slice_num-1)*num_of_division+i),'.png']);
        end
         %imwrite(A/65535,[cd,'\divide\',sprintf('%d',slice_num-k),'.png'],'BitDepth',16);
         fclose('all');
    end
end

