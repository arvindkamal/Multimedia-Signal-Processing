[wavesound,fs] = audioread('s5.wav');
wavesound1=10000* wavesound;

% sound(wavesound);

Tframe = fs*0.025;
Tskip = fs*0.010;
columns = 1+(floor((length(wavesound)-Tframe)/Tskip)); %find how many columns are needed
N=1024; %1024 FFT

win = hamming(Tframe); %w[n]

winFFT = fft(win,N); %fft of hamming window of length

frame = zeros(Tframe,columns); %s[n]
framewin = zeros(Tframe,columns); %sw[n] = s[n]w[n]
framewinFFT= zeros(N,columns); %Sw(W)


for i = 1:columns
   
    frame(1:Tframe,i) = wavesound1(1+(Tskip*(i-1)):Tframe+(Tskip*(i-1)),1);    
    framewin(1:Tframe,i) = frame(1:Tframe,i).*win;
    framewinFFT(1:N,i) = fft(framewin(1:Tframe,i),N); %Sw(w)
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

phi_kp = zeros(40,1);
p_sum_phi_kp = zeros(71,1);
max_p0=zeros(1,298);


for i = 1:columns
    for p = 20:90
        for k = 1:floor(398/p)
            phi_kp(k,1) = autocorr(k*p,i);
        end
        
        p_sum_phi_kp(p-19,1) = p*sum(phi_kp);
    end    
    
    [maxp,tempIndex] = max(p_sum_phi_kp);
    max_p0(i) = tempIndex+19;

end


w0 = 2*pi./max_p0;


for i= 1:columns
    am(1:max_p0(i)-1,i) = ceil(((1:max_p0(i)-1)-1/2).*w0(i).*N/(2*pi));
    bm(1:max_p0(i)-1,i) = floor(((1:max_p0(i)-1)+1/2).*w0(i).*N/(2*pi));
end

%winFFT = fft(win,N); %fft of hamming window of length
%framewinFFT(1:N,i) = fft(framewin(1:Tframe,i),N); %Sw(w)



for i= 1:columns
    
    for j = 1:max_p0(i)-1
        
        top = framewinFFT(am(j,i):bm(j,i),1)  .* conj(winFFT(am(j,i):bm(j,i),1));
        bottom = abs(winFFT(am(j,i):bm(j,i),1)).^2;   
        Am_voiced(j,i) = sum(top)./sum(bottom);
        
        %error_voiced(j,i) = framewinFFT(am(j,i):bm(j,i),1)  .* winFFT(am(j,i):bm(j,i),1);
        
        
        
        
    end    
end




winFFTunvoiced= ones(1024,1);

for i= 1:columns
    
    for j = 1:max_p0(i)-1
        
        top = framewinFFT(am(j,i):bm(j,i),i)  .* conj(winFFTunvoiced(am(j,i):bm(j,i),1));
        bottom = abs(winFFTunvoiced(am(j,i):bm(j,i),1)).^2;   
        Am_unvoiced(j,i) = sum(top)./sum(bottom);
        
    end    
end



    