%%%% Team members: Jiahong Chen (jchen171) & Arvind Kamal ()
%%%% ECE 417 Machine Problem 1
%%%% Due 09/25/18



[wavesound,fs] = audioread('s5.wav');
%wavesound1=10000* wavesound;
wavesound = wavesound(1:21504);
%wavesound = wavesound+aver;

%number of samples per frame
Tframe = fs*0.025;

%how many sample to skip
Tskip = fs*0.010;

%calculate how many total frames are there
columns = 1+(floor((length(wavesound)-Tframe)/Tskip)); %find how many columns are needed

%1024 FFT
N=1024; 

%w[n]
win = hamming(Tframe);

%fft of hamming window of length
winFFT = fft(win,N); 

%s[n] matrix initialization
frame = zeros(Tframe,columns);

%sw[n] = s[n]w[n] matrix initialization
framewin = zeros(Tframe,columns);

%Sw(W) matrix initialization
framewinFFT= zeros(N,columns); %Sw(W)


for i = 1:columns
    
    %generate frames of 200 samples with Tskip of 80 samples
    frame(1:Tframe,i) = wavesound(1+(Tskip*(i-1)):Tframe+(Tskip*(i-1)),1);
    
    %multiply the frames with hamming window
    framewin(1:Tframe,i) = frame(1:Tframe,i).*win;
    
    %Sw(W) which takes the fft of 1024 points of sn[n]
    framewinFFT(1:N,i) = fft(framewin(1:Tframe,i),N); 
end


%%%%%%autocorrelation calculation%%%%%%%

%pad the windown function with 399 zeros
win1=[zeros(399,1);win];

%pad the frames with 399x298 of zeros
frame1 = [zeros(399,columns);frame(:,:)]; 

%pad the window function with 200 zeros then the window fun. then 199 more zeros
win2 = [zeros(200,1);win;zeros(199,1)]; 

%pad the frame with 200x298 zeros then the frames then 199x298 more zeros
frame2 = [zeros(200,columns);frame(:,:);zeros(199,columns)];

%preset zeros for temp1 in the for loop
temp1 = zeros(200,1);

%autocorrelation matrix initialization
autocorr = zeros(399,columns); %preset zeros for autocorrelation 

for i = 1:columns
    for m = 0:398
        for n = 400:599
            
           %autocorrelation equation of w^2[n]s[n] (not summed up yet)          
           temp1(n-399,i) =  win1(n).^2 .* frame1(n,i).* win2(n-m).^2 .* frame2(n-m,i); 
        end
           %summing the autocorrelation arrays
           autocorr(m+1,i) = sum(temp1(:,i));
           
           %%%%% autocorrelaton index in matlab starts from 1 to 399 but
           %%%%% matlab index(1) is -199 and index(399) is 199
    end
end

 


%phi(kp) matrix initialization


%matrix initialization for p estimation from [20,90]
p_sum_phi_kp = zeros(71,1);

%arg max P matrix initialization
max_p0=zeros(1,columns);







for i = 1:columns
    for p = 20:90
        for k = 1:floor(199/p)
            
            phi_kp(k,i) =  autocorr(200-k*p,i);
            phi_kp(50+k,i) =   autocorr(200+k*p,i);
                    
            phi_kp(100,i) = autocorr(200,i);
                        
        end       
        p_sum_phi_kp(p,i) = p*sum(phi_kp(:,i));     
        
    end
    [maxp,tempIndex] = max(p_sum_phi_kp(:,i));
    max_p0(i) = tempIndex;
    
end



% w0 = 2*pi./max_p0;
% 
% 
% for i= 1:columns
%     am(1:floor(max_p0(i))-1,i) = ceil(((1:floor(max_p0(i))-1)-1/2).*w0(i).*N/(2*pi));
%     bm(1:floor(max_p0(i))-1,i) = floor(((1:floor(max_p0(i))-1)+1/2).*w0(i).*N/(2*pi));
% end
% 
% 
% for y= 1:columns
%     for j = 1:max_p0(y)-1
%         for x = am(j,y):bm(j,y)
%             
%             
%             top(1:bm(j,y)-am(j,y)+1,y) = framewinFFT(x,y)  .* conj(winFFT(x,1));
%             
%             bottom(1:bm(j,y)-am(j,y)+1,y) = abs(winFFT(x,1)).^2; 
%         end
%         
%     Am_voiced(j,y) = sum(top(:,y))/sum(bottom(:,y));   
%     end    
%     
% end















%to find the refined p, the original p is +- 2 with step of 0.2
refined_len = length(-2:0.2:2);

%p0 refined matrix initialization
po_refined= zeros(refined_len,columns);

%for loop to generate the different p values
for i = 1: refined_len
    for j = 1: columns
    po_refined(i,j) = max_p0(j)-2+(i-1)*(0.2);
    end
end


totalerror_refined = zeros(refined_len,columns);


%Use the Am_calculation function to use different p0 values 
for i = 1:refined_len
    [RealAm, totalerror,z1] = Am_calculation(po_refined(i,:),framewinFFT);
    totalerror_refined(i,:) = totalerror;   
end
% % 
% 
Index=zeros(1,columns);
Refined_p=zeros(1,columns);


%find the p0 with the smallest error
%Refined_p is the array with p0's that yield smallest error
for i = 1:columns

    [Minp,P] = min(totalerror_refined(:,i));
    Index(i) = P;
    Refined_p(i) = po_refined(Index(i),i);
    
end


%use the refined p0 and run the Am_calculation function again
[Refined_p_Am,Refined_p_error,z2] = Am_calculation(Refined_p(1,:),framewinFFT);






%%%% synethesis part 


for f= 1:columns;
    for n = 80*(f-1)+1 : 80*(f)
        for m= 1: floor(Refined_p(f))-1 
            
            if f == columns
               Am_syn_voiced(n) = real(Refined_p_Am(m,f));
               w0_syn_voiced(n) = 2*pi/Refined_p(f);
            else               
               Am_syn_voiced(n) =  (f-((n-1)/Tskip)) * real(Refined_p_Am(m,f)) + ((n-1)/Tskip-(f-1))* real(Refined_p_Am(m,f+1));
               w0_syn_voiced(n) =  (f-((n-1)/Tskip)) * 2*pi/Refined_p(f) + ((n-1)/Tskip-(f-1))* 2*pi/Refined_p(f+1);                
            end
            
            
            if f==1  && n==1              
                Theta_syn_voiced(n)= m*w0_syn_voiced(n);      
            else               
                Theta_syn_voiced(n) = Theta_syn_voiced(n-1)+m*w0_syn_voiced(n);
            end  
            
            Am_mult_Cos_voiced(m) = Am_syn_voiced(n).* cos(Theta_syn_voiced(n));
               
        end
        Syn_voiced(n,1) = sum(Am_mult_Cos_voiced(m));
    end
end



% 
% 
% % 
% % % %%%Guassian Voice
% % % var1 = normrnd(0, 0.01, 24000,1);
% % % var2 = normrnd(0, 0.05, 24000,1);
% % % var3 = normrnd(0, 0.10, 24000,1);
% % % 
% % noisewave1 = wavesound+var1;
% % noisewave2 = wavesound+var2;
% % noisewave3 = wavesound+var3;
% % 
% % for i = 1:columns
% %     
% %     
% %     frame_noise1(1:Tframe,i) = wavesound(1+(Tskip*(i-1)):Tframe+(Tskip*(i-1)),1);
% %     frame_noise2(1:Tframe,i) = wavesound(1+(Tskip*(i-1)):Tframe+(Tskip*(i-1)),1);
% %     frame_noise3(1:Tframe,i) = wavesound(1+(Tskip*(i-1)):Tframe+(Tskip*(i-1)),1);
% %    
% %     framewin_noise1(1:Tframe,i) = frame(1:Tframe,i).*win;
% %     framewin_noise2(1:Tframe,i) = frame(1:Tframe,i).*win;
% %     framewin_noise3(1:Tframe,i) = frame(1:Tframe,i).*win;
% %     
% %     framewinFFT_noise1(1:N,i) = fft(framewin(1:Tframe,i),N); 
% %     framewinFFT_noise2(1:N,i) = fft(framewin(1:Tframe,i),N); 
% %     framewinFFT_noise3(1:N,i) = fft(framewin(1:Tframe,i),N); 
% % end
% % 
% % for i = 1:columns
% %     for m = 0:398
% %         for n = 400:599
% %             
% %            %autocorrelation equation of w^2[n]s[n] (not summed up yet)          
% %            temp_noise1(n-399,1) =  win1(n).^2 .* frame_noise1(n,i).* win2(n-m).^2 .* frame_noise1(n-m,i); 
% %            temp_noise2(n-399,1) =  win1(n).^2 .* frame_noise2(n,i).* win2(n-m).^2 .* frame_noise2(n-m,i); 
% %            temp_noise3(n-399,1) =  win1(n).^2 .* frame_noise3(n,i).* win2(n-m).^2 .* frame_noise3(n-m,i); 
% %            
% %         end
% %            %summing the autocorrelation arrays
% %            autocorr_noise1(m+1,i) = sum(temp_noise1);
% %            autocorr_noise2(m+1,i) = sum(temp_noise2);
% %            autocorr_noise3(m+1,i) = sum(temp_noise3);
% %            
% %            %%%%% autocorrelaton index in matlab starts from 1 to 399 but
% %            %%%%% matlab index(1) is -199 and index(399) is 199
% %     end
% % end
% % 
% % for i = 1:columns
% %     for p = 20:90
% %         for k = 1:floor(398/p)
% %             
% %             %phi(k*p) which is just autocorrelation index of k*p
% %             phi_kp_noise1(k,1) = autocorr_noise1(k*p,i);
% %             phi_kp_noise2(k,1) = autocorr_noise2(k*p,i);
% %             phi_kp_noise3(k,1) = autocorr_noise3(k*p,i);
% %         end
% %         
% %         %sum all phi(k*p) up and multiply it by p
% %         p_sum_phi_kp_noise1(p-19,1) = p*sum(phi_kp_noise1);
% %         p_sum_phi_kp_noise2(p-19,1) = p*sum(phi_kp_noise2);
% %         p_sum_phi_kp_noise3(p-19,1) = p*sum(phi_kp_noise3);
% %     end
% %     
% %     %find the index of the largest P0 and +19 is for index correlation
% %     [maxp_noise1,tempIndex_noise1] = max(p_sum_phi_kp_noise1);
% %     [maxp_noise2,tempIndex_noise2] = max(p_sum_phi_kp_noise2);
% %     [maxp_noise3,tempIndex_noise3] = max(p_sum_phi_kp_noise3);
% %     
% %     
% %     max_p0_noise1(i) = tempIndex_noise1+19;
% %     max_p0_noise2(i) = tempIndex_noise2+19;
% %     max_p0_noise3(i) = tempIndex_noise3+19;
% %     
% % end
% % 
