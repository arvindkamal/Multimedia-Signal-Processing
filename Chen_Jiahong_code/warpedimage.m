function [warped_image] = warpedimage(w,h1,h2)

neutral_w = 0;
neutral_h1 = 0;
neutral_h2 = 0;
fScale = 1;
min_x = 20;
max_x = 110;
min_y = 30;
max_y = 70;
[mouth] = imread('mouth.jpg');
mouth = mouth';
dim_x = size(mouth,1);
dim_y = size(mouth,2);
[inVertX inVertY, inTriangles] = readmesh('mesh.txt');
dim_tri = size(inTriangles,1);


% for i = 1:size(inVertX)
%     mouth(inVertY(i), inVertX(i)) = 0;
% end

% figure(1)
% imshow(mouth);

[retVertX, retVertY] = interpVert(inVertX, inVertY, neutral_w, neutral_h1, neutral_h2, w, h1, h2, fScale);

X = ones(dim_tri,3,3);
U = ones(dim_tri,3,3);
X_Inverse = zeros(dim_tri,3,3);



%setting up the U, X, and X inverse matrices
for i = 1:dim_tri
    
    temp = inTriangles(i,:);
    
    U(i,1,: )= [inVertX(temp(1)),inVertX(temp(2)), inVertX(temp(3))]; 
    U(i,2,:) = [inVertY(temp(1)),inVertY(temp(2)), inVertY(temp(3))];
    
    X(i,1,:)= [retVertX(temp(1)),retVertX(temp(2)),retVertX(temp(3))];
    X(i,2,:)= [retVertY(temp(1)),retVertY(temp(2)),retVertY(temp(3))];
    
    temp2 = reshape(X(i,:,:),[3,3]);
    X_Inverse(i,:,:) = inv(temp2);
   
end


temp5 = mouth;

for i = 1:dim_x
    for j = 1:dim_y
        
        x = [i,j,1]';
        
        found_k = false;  
        for k  = 1:dim_tri
            lamba = reshape(X_Inverse(k,:,:),[3,3]) * x;
            
            if all( lamba >= 0) && all(lamba<=1)
                found_k = true;
                break;
            end
        end
        
        
        if not(found_k)
                 if and(i > min_x , and(i< max_x, and(j > min_y, j< max_y)))
                      mouth(i,j) = 0;
                 end
        else
            
            temp3 = reshape(U(k,:,:),[3,3]) * lamba;      
            u = temp3(1);
            v = temp3(2);

            m = floor(u); %no int function, so floor function is used
            n = floor(v);
            f = u-m;
            e = v-n;
            
            %bilinear Interpolation formula
            mouth(i,j) = (1-f)*(1-e)*temp5(m,n)+f*(1-e)*temp5(m+1,n)+e*(1-f)*temp5(m,n+1)+e*f*temp5(m+1,n+1);     
        end        
    end
end

warped_image= uint8(mouth');
 
%imshow(warped_image);
end

