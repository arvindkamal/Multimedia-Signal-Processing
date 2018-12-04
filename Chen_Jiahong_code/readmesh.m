function [inVertx, inVerty, inTriangles] = readmesh(meshfile)

mesh_file = csvread(meshfile);

index1 = mesh_file(1);


inVertXY = mesh_file(2:index1+1,1:2);

inVertx = inVertXY(:,1);
inVerty = inVertXY(:,2);

index2 = mesh_file(index1+2);


for i =  (index1+3) : (index1+index2+2)
    inTriangles(i-(index1+2),1) = mesh_file(i,1);
    inTriangles(i-(index1+2),2) = mesh_file(i,2);
    inTriangles(i-(index1+2),3) = mesh_file(i,3);
    
end

end

