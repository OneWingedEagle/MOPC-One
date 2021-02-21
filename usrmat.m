
a=1.17e-0;
nL=8;

Nx=2;
Ny=800;
Nz=2;
X1=-a/2;
Xm=a/2;

Y1=-nL*a/2;
Ym=nL*a/2;

tk=1e-0;

Z1=-tk/2;
Zm=tk/2;

nk=zeros(Nx,Ny,Nz,3);

for iz=1:Nz
 
for iy=1:Ny 
  
for ix=1:Nx
nk(ix,iy,iz)=1+abs(sin(2*pi*4*iy/Ny));
end
end
end


fid = fopen('usr_mat.txt','wt'); 
  
fprintf(fid,'%d\t %f\t %f\n',Nx,X1,Xm);
fprintf(fid,'%d\t %f\t %f\n',Ny,Y1,Ym);
fprintf(fid,'%d\t %f\t %f\n',Nz,Z1,Zm);


for iz=1:Nz
 
for iy=1:Ny 
  
for ix=1:Nx
fprintf(fid,'%f\n',nk(ix,iy,iz));  

end
end
end

fclose(fid);