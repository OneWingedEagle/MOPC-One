function PhC_waveguide
%The program for the computation of the dispersion
%characteristic of the PhC waveguide based upon 2D PhC with square latice
clear all
clc
%%

% The variable mode defines mode of the wave.

% The variable a defines the period of the structure %
a=1e-6;

%The variable r contains elements radius. Here it is defined as a part of the period.
r=0.504*a;

%The variable epsb and Mub contains the information about permittivity of the background.
Epsb=2.72;
gyrb=.0;
%The variables epsa and Mua contains the information about hole.
Epsa=1;

gyra=0.00;

%The variables epsd and Mud contains the information about the defect.

Epsd=1;
gyrd=.0;

gyrDir1=[1 0 0];
gyrDir=[0 0 0];
if(norm(gyrDir1)>0)
    gyrDir=gyrDir1/norm(gyrDir1);
end

%Variable numElem defines the number of PhC elements
numElem=1;

%Number of holes missed
numDef=0;

%The variable precis defines the number of k-vector points between high symmetry points
precis=15;

%The variable nG defines the number of plane waves.
%total number of plane waves may be determined as (nG*2-1)^2
nG=2;

%Setting transversal wave vector to a single value
kx=0;

%kapth segments 1: gama to X, 2: Gama to M, 3: whole path
nKpath=3;

color='r'

PWE(a, r, Epsb, Epsa, Epsd,gyra,gyrb,gyrd,gyrDir, numElem, numDef,nG,nKpath,precis, kx,color, 1);

%Setting transversal wave vector to a range of values within Brillouin zone
kx=0:(pi/a)/2:pi/a;

gyra=0;
gyrb=0;
gyrd=0;

%PWE(a, r, Epsb, Epsa, Epsd,gyra,gyrb,gyrd,gyrDir, numElem, numDef,nG,nKpath,precis, kx,'k', 1);




    function PWE(a, r,Epsb, Epsa, Epsd,gyra,gyrb,gyrd,gyrDir, numElem,numDefect,nG,nKpath,precis, kx_param, color , linewidth)
        
        nGx=nG*numElem;
        nGy=nG;
        % nGx and nGy define the number of plane waves in x and y direction
        % respectively
        
        cmp=3;
        %discretization mesh elements. each cell is discretized to a grid of
        %ngrid*ngrid square elememts for numerical compuation of Fourier series.
        
        ngrid=60;
        
        epsTensInv=zeros(cmp,cmp,ngrid,ngrid);
        
        %The following loop carries out the definition of the unit cell.
        for n=0:numElem-1
            nx=1;
            for countX=-a/2:a/ngrid:a/2
                ny=1;
                for countY=-a/2:a/ngrid:a/2
                    %Following condition allows to define the circle with of radius r
                    if(countX^2+countY^2<r^2)
                        if(n<numDefect)
                            eps=Epsd;
                            gyr=gyrd*gyrDir;
                            
                        else
                            eps=Epsa;
                            gyr=gyra*gyrDir;
                        end
                    else
                        eps=Epsb;
                        gyr=gyrb*gyrDir;
                    end
                    
                    epsTens=zeros(cmp,cmp);
                    
                    
                    epsTens(1,1)=eps;
                    epsTens(2,2)=eps;
                    epsTens(3,3)=eps;
                    
                    epsTens(1,2)=-1i*gyr(3)*eps;
                    epsTens(2,1)=1i*gyr(3)*eps;
                    epsTens(1,3)=1i*gyr(2)*eps;
                    epsTens(3,1)=-1i*gyr(2)*eps;
                    epsTens(2,3)=-1i*gyr(1)*eps;
                    epsTens(3,1)=1i*gyr(1)*eps;
                    
                    
                    epsTensInv(:,:,nx+ngrid*n,ny)=inv(epsTens);
                    
                    xSet(nx+ngrid*n)=countX+a*n;
                    ySet(ny)=countY;
                    ny=ny+1;
                end
                nx=nx+1;
            end
        end
        
        
        
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
        
         
        MtNt=(length(xSet)-1)*(length(ySet)-1);
        
        bx=2*pi/a/numElem;
        by=2*pi/a;
        
        LEFT=zeros(cmp,cmp,4*nGx+1,4*nGy+1);
        %OFFDIAG=zeros(4*nGx+1,4*nGy+1);
        %The next loop computes the Fourier expansion coefficients
        for dGx=-2*nGx:2*nGx
            for dGy=-2*nGy:2*nGy
                dGxp=dGx+1+2*nGx;
                dGyp=dGy+1+2*nGy;
                
                for nx=1:length(xSet)-1
                    for ny=1:length(ySet)-1
                        x=xSet(nx);
                        y=ySet(ny);
                        tt=dGx*bx*x+dGy*by*y;
                        
                        LEFT(:,:,dGxp,dGyp)=LEFT(:,:,dGxp,dGyp)+epsTensInv(:,:,nx,ny)*exp(-1i*(tt));
                        
                    end
                end
                
                
            end
        end
        
        
        LEFT= LEFT/MtNt;
        
        
        numG=0;
        for Gx=-nGx:nGx
            for Gy=-nGy:nGy
                numG=numG+1;
            end
        end
        
        M=  zeros(size(kv,1),cmp*numG,cmp*numG);
        
        countG=0;
        
        for Gx=-nGx:nGx
            for Gy=-nGy:nGy
                countG=countG+1;
                countG3=cmp*(countG-1);
                countG1=0;
                for Gx1=-nGx:nGx
                    
                    dGxp=2*nGx+1+Gx-Gx1;
                    
                    for Gy1=-nGy:nGy
                        
                        dGyp=2*nGy+1+Gy-Gy1;
                        
                        countG1=countG1+1;
                        
                        countG4=cmp*(countG1-1);
                        
                        for countK=1:size(kv,1)
                            kx=kv(countK,1);
                            ky=kv(countK,2);
                            
                            Akg=zeros(3,3);
                            
                            Akg(1,1)=(ky+Gy1*by)*(ky+Gy*by);
                            Akg(1,2)=-(kx+Gx1*bx)*(ky+Gy*by);
                            Akg(2,1)=-(kx+Gx*bx)*(ky+Gy1*by);
                            Akg(2,2)=(kx+Gx1*bx)*(kx+Gx*bx);
                            
                            Akg(3,3)=Akg(1,1)+ Akg(2,2);
                            
                            KA=LEFT(:,:, dGxp, dGyp)*Akg;
                            
                            
                            for ii=1:cmp
                                for jj=1:cmp
                                    M(countK,countG3+ii,countG4+jj)= KA(ii,jj);
                                    
                                end
                            end
                            
                            
                        end
                    end
                    
                    
                end
            end
            
        end
        
        
        disp('Solving eigenvalue problem ....')
        t1=cputime;
    
        %The computation of eigen-states is also carried
        disp('size kv')
        size(kv,1);
        for countK=1:size(kv,1)
          countK
            %Taking the matrix differential operator for current wave vector.
            MM(:,:)=M(countK,:,:);
            
            %Computing the eigen-vectors and eigen-states of the matrix
            [Q D]=eig(MM);
            
            V=diag(D);

            [dispe(:,countK), index]=sort(abs(sqrt(V)*a/2/pi));
        
 
     end;
         
        t2=cputime-t1;
        disp('-----');
        sprintf('%s  %8f', 'Eigenvalue calculation completed within', t2, '  seconds.')
        
        %% plot %%
        %%%%%%%%%%
        % close all;
        nz=0;
        
        if(Nk>1)
            for i=1:size(V)
                if(dispe(i,2)<1e-4)
                    nz=nz+1;
                end;
            end
        end

nCurve=15;


figure(1);
hold on;


    x=0;
     if(nKpath==3)
         x=1;
     end;
    
    for u=nz+1:nz+nCurve
  yy(1:Nk)=zeros(Nk,1);
     yy(1)=dispe(u,1);
          
 
     for m=2:Nk-x
        
           yy(m)=dispe(u,m);
      
     end
         

    if(nKpath==3)
        yy(Nk)=yy(1);
    end
    
    kindex=1:length(yy);
    plot(kindex,yy,color,'LineWidth',linewidth);
    
    end
 

ylabel('Frequency \omegaa/2\pic','FontSize',18);
xlabel('Propagation constant \beta, m^{-1}','FontSize',18);
title('The dispersion diagram of the PhC waveguide','FontSize',10)
axis([1,Nk,0,0.6]);

grid on;
str={'Gamma','X','M','Gamma'};

set(gca,...
    'xtick',[1,precis+1,2*precis+1,3*precis+1],...
    'xticklabel',str);
drawnow;