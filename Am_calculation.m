function [RealAm, totalerror, voiced_unvoiced] = Am_calculation(max_p0,Sw)

w0 = 2*pi./max_p0;
win = hamming(200); %w[n]
winFFT = fft(win,1024);
N =1024;
columns = 267;

for i= 1:columns
    am(1:floor(max_p0(i))-1,i) = ceil(((1:floor(max_p0(i))-1)-1/2).*w0(i).*N/(2*pi));
    bm(1:floor(max_p0(i))-1,i) = floor(((1:floor(max_p0(i))-1)+1/2).*w0(i).*N/(2*pi));
end




for y= 1:columns
    for j = 1:max_p0(y)-1
        for x = am(j,y):bm(j,y)
            
            
            top(1:bm(j,y)-am(j,y)+1,y) = Sw(x,y)  .* conj(winFFT(x,1));
            
            bottom(1:bm(j,y)-am(j,y)+1,y) = abs(winFFT(x,1)).^2; 
        end
        
    Am_voiced(j,y) = sum(top(:,y))/sum(bottom(:,y));   
    end    
    
end




for i= 1:columns
    for j = 1:max_p0(i)-1    
        for x = am(j,i):bm(j,i)
        
        error1 = abs(Sw(x,1) - (Am_voiced(j,i)).*(winFFT(x,1))).^2;
        end
        
        error_voiced(j,i) = 1/(2*pi) *sum(error1);
    end
    
end




winFFTunvoiced= ones(1024,1);

for i= 1:columns
    
    for j = 1:max_p0(i)-1
        for x = am(j,i):bm(j,i)
        
            top = Sw(x,i)  .* conj(winFFTunvoiced(x,1));
            bottom = abs(winFFTunvoiced(x,1)).^2;   
        
        end
        Am_unvoiced(j,i) = sum(top)./sum(bottom);
    end  
    
end




for i= 1:columns
    
    for j = 1:max_p0(i)-1
        
        for x = am(j,i):bm(j,i)
        error2 = abs(real(Sw(x,1)) - real(Am_unvoiced(j,i)).* real(winFFTunvoiced(x,1)).^2);
        
        end
        
        error_unvoiced(j,i) = 1/(2*pi) *sum(error2);
    end
    
end







for j = 1:length(error_voiced(:,1))
    for i = 1:columns
        
        if (error_voiced(j,i)>error_unvoiced(j,i))
            
            
            RealAm(j,i) = Am_unvoiced(j,i);
            Realerror(j,i) = error_unvoiced(j,i);
            voiced_unvoiced(j,i)= 0;
    
            
        else
            
            RealAm(j,i) = Am_voiced(j,i);
            Realerror(j,i) = error_voiced(j,i);
            voiced_unvoiced(j,i)= 1;
            
        end     
    end
end



totalerror = zeros(1,columns);

for i = 1:columns
    
    totalerror(i) = sum(Realerror(:,i));
    
end



end

