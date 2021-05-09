% in the name of allah
% in the name of allah
% MATLAB Program for Computation of Dispersion Diagrams of 2D PhC Waveguide

function PhC_Igor
%The program for the computation of the dispersion
%characteristic of the PhC waveguide based upon 2D PhC with square latice
clear all
clc
%%
% The variable a defines the period of the structure %

%mode=0  defines the TE mode and mode=1 defines TM mode %
mode=0
a=1e-6;
a=1.17e-6
%The variable r contains elements radius. Here it is defined as a part of the period.

r=0.4*a;
r=0.25e-6
%The variable eps1 contains the information about the relative background material permittivity.
%Eps1=2.72;
Eps1=5.6;

%The variable eps2 contains the information about permittivity of the elements composing PhC.
Eps2=1.0;

%Variable numElem defines the number of PhC elements
numElem=1;

%Number of holes missed
numDef=0;

%Setting transversal wave vector to a single value
kx=0;
%kx=0:(pi/a)/10:pi/a;

color='-b';
if(mode==1)
    color='ok';
end


%PWE(mode, a, r, Eps1, Eps2, numElem, numDef, kx,color, 2);
%numElem=1;
%numDef=0;

mode=1;
color='-r';

PWE(mode, a, r, Eps1, Eps2, numElem, numDef, kx,color, 2);
%Setting transversal wave vector to a range of values within Brillouin zone
kx=0:(pi/a)/20:pi/a;

%PWE(mode, a, r, Eps1, Eps2, numElem, numDef, kx, 'r', 2);

function PWE(mode, a, r, Eps1, Eps2, numElem,numDefect, kx_param, color, linewidth)

%The variable precis defines the number of k-vector points between high symmetry points
precis=15;
%The variable nG defines the number of plane waves.
%total number of plane waves may be determined as (nG*2-1)^2
nG=3;

%kapth segments 1: gama to X, 2: Gama to M, 3: whole path
nKpath=3;

%discretization mesh elements.
precisStruct=60;

%The following loop carries out the definition of the unit cell.
for n=0:numElem-1
nx=1;
for countX=-a/2:a/precisStruct:a/2
ny=1;
for countY=-a/2:a/precisStruct:a/2
%Following condition allows to define the circle with of radius r
if((sqrt(countX^2+countY^2)<r)&&(n>=numDefect))
%Setting the value of the inversed dielectric function to the mesh node
struct(nx+precisStruct*n,ny)=1/Eps2;
%Saving the node coordinate
xSet(nx+precisStruct*n)=countX+a*n;
ySet(ny)=countY;
else
struct(nx+precisStruct*n,ny)=1/Eps1;
xSet(nx+precisStruct*n)=countX+a*n;
ySet(ny)=countY;
end
ny=ny+1;
end
nx=nx+1;
end
end



%The computation of the area of the mesh cell.
dS=(a/precisStruct)^2;
%Forming 2D arrays of nods coordinates
Lx=length(xSet)-1;
Ly=length(ySet)-1;

[xMesh,yMesh]=meshgrid(xSet(1:Lx),ySet(1:Ly));
structMesh=(struct(1:Lx,1:Ly)*dS/((max(xSet)-min(xSet))*(max(ySet)-min(ySet))))';
b1=[2*pi/a 0]./numElem;
b2=[0 2*pi/a];

%Defining the k-path within the Brilluoin zone.

kxm=pi/a;
dkx=kxm/precis;

Lk=nKpath*precis+1;
if(nKpath==3) 
    Lk=Lk-1;
end
kv=zeros(Lk,2);
for i=2:Lk
    if(i<=precis+1)
        kv(i,2)=kv(i-1,2)+dkx;
    elseif(i<=2*precis+1)
         kv(i,1)=kv(i-1,1)+dkx;
          kv(i,2)=kv(i-1,2);
    else
        kv(i,1)=kv(i-1,1)-dkx;
        kv(i,2)= kv(i-1,2)-dkx;
    end
end


Nk=size(kv,1);
if(nKpath==3)
    Nk=Nk+1;
end

numG=1;
for Gx=-nG*numElem:nG*numElem
for Gy=-nG:nG
G(numG,1)=Gx*b1(1)+Gy*b2(1);
G(numG,2)=Gx*b1(2)+Gy*b2(2);
numG=numG+1;
end
end

%The next loop computes the Fourier expansion coefficients
for countG=1:numG-1
for countG1=1:numG-1
CN2D_N(countG,countG1)=sum(sum(structMesh.*exp(1i*((G(countG,1)-G(countG1,1))*xMesh+(G(countG,2)-G(countG1,2))*yMesh))));
end
end


%The computation is carried out for each of earlier defined
for kx_count=1:length(kx_param)

for countG=1:numG-1
for countG1=1:numG-1
for countK=1:size(kv,1)
    kx=kv(countK,1);
     ky=kv(countK,2);
     
    w1=0;
   
    if(mode==0)
          w1=((kx+G(countG,1))*(kx+G(countG1,1))+(ky+G(countG,2))*(ky+G(countG1,2))); %TE
    else
     w1=  (sqrt((kx+G(countG,1))^2+(ky+G(countG,2))^2)*sqrt((kx+G(countG1,1))^2+(ky+G(countG1,2))^2));%TM
    end
  
   M(countK,countG,countG1)=CN2D_N(countG,countG1)*w1;
 
end
end
end

%The computation of eigen-states is also carried
for countK=1:size(kv,1)
    
%Taking the matrix differential operator for current wave vector.
MM(:,:)=M(countK,:,:);

%Computing the eigen-vectors and eigen-states of the matrix
   V=eig(MM);

%Transforming matrix eigen-states to the form of normalized frequency
dispe(:,countK)=sort(sqrt(V)*a/2/pi);
end
Nk=size(dispe,2);
nCurve=15;
result=zeros(Nk,nCurve+1);
%% plot %%
%%%%%%%%%%
figure(1);
hold on;
for u=1:nCurve
  yy=abs(dispe(u,:));
        if(nKpath==3)
            yy=[yy yy(1)];
        end
        
     for m=1:Nk   
       
       result(m,u+1)=yy(m);
     end
        plot(yy,color,'LineWidth',linewidth);
        
        
end

 
      for m=1:Nk   
       result(m,1)=m;
     end


save('bandsTM.txt', 'result', '-ascii')

ylabel('Frequency \omegaa/2\pic','FontSize',18);
xlabel('Propagation constant \beta, m^{-1}','FontSize',18);
title('The dispersion diagram of the PhC waveguide','FontSize',10)
axis([1,Nk,0,0.6]);

grid on;
str={'G','X','M','G'};
     set(gca,...
          'xtick',[1,precis+1,2*precis+1,3*precis+1],...
          'xticklabel',str);

drawnow;
end



