%load data and sound
load('ECE417_MP4_AV_Data.mat');
test_sound = audioread('test.wav');



img_len = size(av_validate.visual,2);


for i =1:480
    
 w = av_validate.visual(1,i);
 h1 = av_validate.visual(2,i);
 h2 = av_validate.visual(3,i);
 
    %generate jpeg based on the order and naming
 
     if i<10

        [warped_image] = warpedimage(w,h1,h2);    
        imwrite(warped_image,strcat('/Users/jiahong/Desktop/mp4picture/test_000',num2str(i),'.jpeg'));
     elseif i<100
        [warped_image] = warpedimage(w,h1,h2);    
        imwrite(warped_image,strcat('/Users/jiahong/Desktop/mp4picture/test_00',num2str(i),'.jpeg'));
     else    
         [warped_image] = warpedimage(w,h1,h2);    
        imwrite(warped_image,strcat('/Users/jiahong/Desktop/mp4picture/test_0',num2str(i),'.jpeg'));
     end


end