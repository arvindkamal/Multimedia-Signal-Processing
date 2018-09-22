[wavesound,fs] = audioread('s5.wav');
wavesound1(:,1)=10000 .* wavesound(:,1);

% sound(wavesound);

Tframe = fs*0.025;
Tskip = fs*0.010;
columns = 1+(floor((length(wavesound)-Tframe)/Tskip)); %find how many columns are needed
N=1024; %1024 FFT

win = hamming(Tframe); %w[n]
frame = zeros(Tframe,columns); %s[n]
framewin = zeros(Tframe,columns); %sw[n] = s[n]w[n]
framewinFFT= zeros(N,columns); %Sw(W)

for i = 1:columns
   
    frame(1:Tframe,i) = wavesound1(1+(Tskip*(i-1)):Tframe+(Tskip*(i-1)),1);    
    framewin(1:Tframe,i) = frame(1:Tframe,i).*win;
    framewinFFT(1:N,i) = fft(framewin(1:Tframe,i),N);
end



%%%%%%autocorrelation calculation%%%%%%%
win1=[zeros(399,1);win]; %pad the windown function with 399 zeros
frame1 = [zeros(399,298);frame(:,:)]; %pad the frames with 399x298 of zeros
win2 = [zeros(200,1);win;zeros(199,1)]; %pad the window function with 200 zeros then the window fun. then 199 more zeros
frame2 = [zeros(200,298);frame(:,:);zeros(199,298)]; %pad the frame with 200x298 zeros then the frames then 199x298 more zeros

temp1 = zeros(200,1); %preset zeros for temp1 in the for loop
autocorr = zeros(399,298); %preset zeros for autocorrelation 

for i = 1:columns
    for m = 0:398
        for n = 400:599
            temp1(n-399,1) =  win1(n).^2 .* frame1(n,i).* win2(n-m).^2 .* frame2(n-m,i); 
        end
           autocorr(m+1,i) = sum(temp1);
           %%%%% autocorrelaton index in matlab starts from 1 to 399 but
           %%%%% matlab index(1) is -199 and index(399) is 199
    end
end

 

%spectrogram of original sound
%spectrogram(wavesound)

% test = autocorr(:,100);
% 
% %for p = 20:1:90
%     for k = 1:floor(398/20)
%         
%         result(k,1) = test(20*k) 
%         
%         
%     end 
%end
    
    
    







    