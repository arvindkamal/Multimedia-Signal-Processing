%%%% Team members: Jiahong Chen (jchen171) & Arvind Kamal ()
%%%% ECE 417 Machine Problem 1
%%%% Due 09/25/18



[wavesound,fs] = audioread('s5.wav');
%wavesound1=10000* wavesound;

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
frame1 = [zeros(399,298);frame(:,:)]; 

%pad the window function with 200 zeros then the window fun. then 199 more zeros
win2 = [zeros(200,1);win;zeros(199,1)]; 

%pad the frame with 200x298 zeros then the frames then 199x298 more zeros
frame2 = [zeros(200,298);frame(:,:);zeros(199,298)];

%preset zeros for temp1 in the for loop
temp1 = zeros(200,1);

%autocorrelation matrix initialization
autocorr = zeros(399,298); %preset zeros for autocorrelation 

for i = 1:columns
    for m = 0:398
        for n = 400:599
            
           %autocorrelation equation of w^2[n]s[n] (not summed up yet)          
           temp1(n-399,1) =  win1(n).^2 .* frame1(n,i).* win2(n-m).^2 .* frame2(n-m,i); 
        end
           %summing the autocorrelation arrays
           autocorr(m+1,i) = sum(temp1);
           
           %%%%% autocorrelaton index in matlab starts from 1 to 399 but
           %%%%% matlab index(1) is -199 and index(399) is 199
    end
end

 


%phi(kp) matrix initialization
phi_kp = zeros(40,1);

%matrix initialization for p estimation from [20,90]
p_sum_phi_kp = zeros(71,1);

%arg max P matrix initialization
max_p0=zeros(1,298);



for i = 1:columns
    for p = 20:90
        for k = 1:floor(398/p)
            
            %phi(k*p) which is just autocorrelation index of k*p
            phi_kp(k,1) = autocorr(k*p,i);
        end
        
        %sum all phi(k*p) up and multiply it by p
        p_sum_phi_kp(p-19,1) = p*sum(phi_kp);
    end
    
    %find the index of the largest P0 and +19 is for index correlation
    [maxp,tempIndex] = max(p_sum_phi_kp);
    max_p0(i) = tempIndex+19;

end


% to find the refined p, the original p is +- 2 with step of 0.2
refined_len = length(-2:0.2:2);

%p0 refined matrix initialization
po_refined= zeros(refined_len,columns);

%%for loop to generate the different p values
for i = 1: refined_len
    for j = 1: columns
    po_refined(i,j) = max_p0(j)-2+(i-1)*(0.2);
    end
end


totalerror_refined = zeros(refined_len,columns);

for i = 1:refined_len
    [RealAm, totalerror] = Am_calculation(po_refined(i,:),framewinFFT);
    totalerror_refined(i,:) = totalerror; 
    
end

Index=zeros(1,columns);
Refined_p=zeros(1,columns);

for i = 1:columns

    [Minp,P] = min(totalerror_refined(:,i));
    Index(i) = P;
    Refined_p(i) = po_refined(Index(i),i);
    
end


[Refined_p_Am,Refined_p_error] = Am_calculation(Refined_p(1,:),framewinFFT);







    %[RealAm, totalerror] = Am_calculation(max_p0,framewinFFT);





%Tskip is 80 sample
%f= 0:columns-1;

%w0 = (f+1)-(n/Tskip)*(2*pi)/P(index) +(n/Tskip-f)*(2*pi)/P(index+1)








% w0 = 2*pi./max_p0;
% 
% 
% for i= 1:columns
%     am(1:max_p0(i)-1,i) = ceil(((1:max_p0(i)-1)-1/2).*w0(i).*N/(2*pi));
%     bm(1:max_p0(i)-1,i) = floor(((1:max_p0(i)-1)+1/2).*w0(i).*N/(2*pi));
% end
% 
% %winFFT = fft(win,N); %fft of hamming window of length
% %framewinFFT(1:N,i) = fft(framewin(1:Tframe,i),N); %Sw(w)
% 
% 
% 
% for i= 1:columns
%     
%     for j = 1:max_p0(i)-1
%         
%         top = framewinFFT(am(j,i):bm(j,i),1)  .* conj(winFFT(am(j,i):bm(j,i),1));
%         bottom = abs(winFFT(am(j,i):bm(j,i),1)).^2;   
%         Am_voiced(j,i) = sum(top)./sum(bottom);
%         
%         
%         %error1 = abs(framewinFFT(am(j,i):bm(j,i),1) - Am_voiced(j,i).* winFFT(am(j,i):bm(j,i),1)).^2;
%         error1 = abs(real(framewinFFT(am(j,i):bm(j,i),1)) - real(Am_voiced(j,i)).* real(winFFT(am(j,i):bm(j,i),1)).^2);
%         error_voiced(j,i) = 1/(2*pi) *sum(error1);
%         
%         
%     end    
% end
% 
% 
% 
% 
% winFFTunvoiced= ones(1024,1);
% 
% for i= 1:columns
%     
%     for j = 1:max_p0(i)-1
%         
%         top = framewinFFT(am(j,i):bm(j,i),i)  .* conj(winFFTunvoiced(am(j,i):bm(j,i),1));
%         bottom = abs(winFFTunvoiced(am(j,i):bm(j,i),1)).^2;   
%         Am_unvoiced(j,i) = sum(top)./sum(bottom);
%         
%         %error2 = abs(framewinFFT(am(j,i):bm(j,i),1) - Am_unvoiced(j,i).* winFFTunvoiced(am(j,i):bm(j,i),1)).^2;
%         error2 = abs(real(framewinFFT(am(j,i):bm(j,i),1)) - real(Am_unvoiced(j,i)).* real(winFFTunvoiced(am(j,i):bm(j,i),1)).^2);
%         error_unvoiced(j,i) = 1/(2*pi) *sum(error2);
%         
%     end    
% end
% 
% for i = 1:84
%     for j = 1:298
%         
%         if (error_voiced(i,j)>error_unvoiced(i,j))
%             
%             
%             RealAm(i,j) = Am_unvoiced(i,j);
%             Realerror(i,j) = error_unvoiced(i,j);
%         else
%             
%             RealAm(i,j) = Am_voiced(i,j);
%             Realerror(i,j) = error_voiced(i,j);
%             
%         end     
%     end
% end
% 
% 
% 
% for i = 1:298
%     
%     totalerror(i) = sum(Realerror(:,i));
%     
% end
% 
%     