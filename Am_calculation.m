function [RealAm, totalerror] = Am_calculation(max_p0,Sw)

w0 = 2*pi./max_p0;
win = hamming(200); %w[n]
winFFT = fft(win,1024);
N =1024;
columns = 298;

for i= 1:columns
    am(1:floor(max_p0(i))-1,i) = ceil(((1:floor(max_p0(i))-1)-1/2).*w0(i).*N/(2*pi));
    bm(1:floor(max_p0(i))-1,i) = floor(((1:floor(max_p0(i))-1)+1/2).*w0(i).*N/(2*pi));
end


for i= 1:columns
    for j = 1:max_p0(i)-1        
        top = Sw(am(j,i):bm(j,i),1)  .* conj(winFFT(am(j,i):bm(j,i),1));
        bottom = abs(winFFT(am(j,i):bm(j,i),1)).^2;   
        Am_voiced(j,i) = sum(top)./sum(bottom);
        
        
        %error1 = abs(framewinFFT(am(j,i):bm(j,i),1) - Am_voiced(j,i).* winFFT(am(j,i):bm(j,i),1)).^2;
        error1 = abs(real(Sw(am(j,i):bm(j,i),1)) - real(Am_voiced(j,i)).* real(winFFT(am(j,i):bm(j,i),1)).^2);
        error_voiced(j,i) = 1/(2*pi) *sum(error1);
        
    end    
end





winFFTunvoiced= ones(1024,1);

for i= 1:columns
    
    for j = 1:max_p0(i)-1
        
        top = Sw(am(j,i):bm(j,i),i)  .* conj(winFFTunvoiced(am(j,i):bm(j,i),1));
        bottom = abs(winFFTunvoiced(am(j,i):bm(j,i),1)).^2;   
        Am_unvoiced(j,i) = sum(top)./sum(bottom);
        
        
        %error2 = abs(framewinFFT(am(j,i):bm(j,i),1) - Am_unvoiced(j,i).* winFFTunvoiced(am(j,i):bm(j,i),1)).^2;
        error2 = abs(real(Sw(am(j,i):bm(j,i),1)) - real(Am_unvoiced(j,i)).* real(winFFTunvoiced(am(j,i):bm(j,i),1)).^2);
        error_unvoiced(j,i) = 1/(2*pi) *sum(error2);
        
    end    
end


for i = 1:length(error_voiced(:,1))
    for j = 1:298
        
        if (error_voiced(i,j)>error_unvoiced(i,j))
            
            
            RealAm(i,j) = Am_unvoiced(i,j);
            Realerror(i,j) = error_unvoiced(i,j);
            
            
            
        else
            
            RealAm(i,j) = Am_voiced(i,j);
            Realerror(i,j) = error_voiced(i,j);
            
        end     
    end
end



totalerror = zeros(1,298);

for i = 1:298
    
    totalerror(i) = sum(Realerror(:,i));
    
end



end

