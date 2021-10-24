clear 
clc
ind=0;
disp('This is a course project for the solution of flow field in a cavity for Reynolds numbers 1 and 100 using the Finite Volume approach')
disp('*******************************************************')
Re=input('Enter Reynolds number : ');
U0=1;
L=1;
% imax=input('Enter number of nodes in X direction : ');
% jmax=input('Enter number of nodes in Y direction : ');
imax=input('Enter number of nodes in each direction : ');
jmax=imax;
URFm=0.8;
%% Grid
dx=L/(imax-1);
dy=L/(jmax-1);
for i=1:imax
    x(i)=(i-1)*dx;
end
for j=1:jmax
   y(j)=(j-1)*dy;
end
%% Boundary conditions & initialization
Vstar=zeros(jmax,imax);
V=zeros(jmax,imax);
Ustar=zeros(jmax,imax);
U=zeros(jmax,imax);
Ustar(jmax,:)=1;
U(jmax,:)=1;
Pstar=zeros(jmax,imax);
P=zeros(jmax,imax);
Pc=zeros(jmax,imax);
Tu=zeros(jmax,imax);
Tv=zeros(jmax,imax);
%% LOOP
ss=1;
while ss>0
%% X-momentum
for i=3:imax-1
    for j=2:jmax-1
            Fe=(Ustar(j,i)+Ustar(j,i+1))/2*dy;     %(UP+UE)/2
            Fw=(Ustar(j,i)+Ustar(j,i-1))/2*dy;         %(UP+UW)/2
            Fn=(Vstar(j+1,i)+Vstar(j+1,i-1))/2*dx;       %(UP+UN)/2
            Fs=(Vstar(j,i)+Vstar(j,i-1))/2*dx;                %(UP+US)/2
            De=dy/Re/dx;  
            Dw=dy/Re/dx;
            Dn=dx/Re/dy;
            Ds=dx/Re/dy;
            if   j==2 
                Ds=2*Ds;
            end
            if   j==jmax-1
                Dn=2*Dn;
            end
            aE=De-Fe/2;
            aW=Dw+Fw/2;
            aN=Dn-Fn/2;
            aS=Ds+Fs/2;
            aP=aE+aW+aN+aS+(Fe+Fn-Fw-Fs);
            Tu(j,i)=dy/aP;  %for pressure correction equation
            SU=(Pstar(j,i-1)-Pstar(j,i))*dy;
            Sum=aE*U(j,i+1)+aW*U(j,i-1)+aN*U(j+1,i)+aS*U(j-1,i)+SU;
            U(j,i)=URFm*Sum/aP+(1-URFm)*Ustar(j,i);
    end
end
%% Y-momentum
  for j=3:jmax-1
        for i=2:imax-1
            Fe=(Ustar(j,i+1)+Ustar(j-1,i+1))/2*dy;
            Fw=(Ustar(j,i)+Ustar(j-1,i))/2*dy;
            Fn=(Vstar(j,i)+Vstar(j+1,i))/2*dx;
            Fs=(Vstar(j,i)+Vstar(j-1,i))/2*dx;
            De=dy/Re/dx;
            Dw=dy/Re/dx;
            Dn=dx/Re/dy;
            Ds=dx/Re/dy;
            if i==2
                Dw=2*Dw;
            end
            if i==imax-1
                De=2*De;
            end
            aE=De-Fe/2;
            aW=Dw+Fw/2;
            aN=Dn-Fn/2;
            aS=Ds+Fs/2;
            aP=aE+aW+aN+aS+(Fe+Fn-Fw-Fs);
            Tv(j,i)=dx/aP;    %for pressure correction equation
            SU=(Pstar(j-1,i)-Pstar(j,i))*dx;
            Sum=aE*V(j,i+1)+aW*V(j,i-1)+aN*V(j+1,i)+aS*V(j-1,i)+SU;
            V(j,i)=URFm*Sum/aP+(1-URFm)*Vstar(j,i);
        end
  end
 %% Continuety
     for i=2:imax-1
        for j=2:jmax-1
            aE=Tu(j,i+1)*dy;
            aW=Tu(j,i)*dy;
            aN=Tv(j+1,i)*dx;
            aS=Tv(j,i)*dx;
            aP=aE+aW+aN+aS;
            SU=(U(j,i)-U(j,i+1))*dy+(V(j,i)-V(j+1,i))*dx;
            sum=aE*Pc(j,i+1)+aW*Pc(j,i-1)+aN*Pc(j+1,i)+aS*Pc(j-1,i);
            Pc(j,i)=(sum+SU)/aP;
        end
     end
     %% Pressure correction
    for i=2:imax-1
        for j=2:jmax-1
            P(j,i)=Pstar(j,i)+Pc(j,i);
        end
    end
    %% Velocity correction
 for i=3:imax-1
     for j=2:jmax-1
         U(j,i)=U(j,i)+Tu(j,i)*(Pc(j,i-1)-Pc(j,i));
     end
 end 
for i=2:imax-1
     for j=3:jmax-1
            V(j,i)=V(j,i)+Tv(j,i)*(Pc(j-1,i)-Pc(j,i));
     end
end
%% Convergence  
ss=0;
sp=0;
su=0;
sv=0;
for i=2:imax-1
    for j=2:jmax-1
        if abs((U(j,i+1)-U(j,i))*dy+(V(j+1,i)-V(j,i))*dx)<=10^(-6)
            s=0;
            ss=ss+s;
        else
            s=1;
            ss=ss+s;
        end
            su=su+(U(j,i)-Ustar(j,i))^2;
            sv=sv+(V(j,i)-Vstar(j,i))^2;
            sp=sp+(P(j,i)-Pstar(j,i))^2;
    end
end
Pc=zeros(jmax,imax);
Pstar=P;
Vstar=V;
Ustar=U;
ind=ind+1;
su=sqrt(su);
sv=sqrt(sv);
sp=sqrt(sp);
subplot(1,3,1)
plot(ind,su,'o')
title('X-momentum Residual')
hold on
pause(0.001)
subplot(1,3,2)
plot(ind,sv,'o')
title('Y-momentum Residual')
hold on
pause(0.001)
subplot(1,3,3)
plot(ind,sp,'o')
title('Pressure Residual')
hold on
pause(0.001)
end
%% CFD post
   %% Calculation of boundary pressure
   for i=2:imax-1
       P(1,i)=P(2,i)-(1/(Re*dy))*(2*U(1,i)-5*U(2,i)+4*U(3,i)-U(4,i));
       P(jmax,i)=P(jmax-1,i)+(1/(Re*dy))*(2*U(jmax,i)-5*U(jmax-1,i)+4*U(jmax-2,i)-U(jmax-3,i));
   end
   for j=2:jmax-1
       P(j,1)=P(j,2)-(1/(Re*dx))*(2*V(j,1)-5*V(j,2)+4*V(j,3)-V(j,4))*dx;
       P(j,imax)=P(j,imax-1)+((Re*dx))*(2*V(j,imax)-5*V(j,imax-1)+4*V(j,imax-2)-V(j,imax-3))*dx;
   end
   P(1,1)=0.5*(P(1,2)+P(2,1));
   P(1,imax)=0.5*(P(1,imax-1)+P(2,imax));
   P(jmax,1)=0.5*(P(jmax,2)+P(jmax-1,1));
   P(jmax,imax)=0.5*(P(jmax,imax-1)+P(jmax-1,imax));
 %% Calculation of stream function
 psi=zeros(jmax,imax);
 opsi=psi;
 cc=1;
 while cc>0
     cc=0;
     for i=2:imax-1
         for j=jmax-1:-1:2
                 psi(j,i)=psi(j+1,i)-U(j,i)*dy;
         end
     end
  for i=2:imax-1
      for j=2:jmax-1
          if abs(opsi(j,i)-psi(j,i))<=(10^-7)
              c=0;
              cc=cc+c;
          else
              c=1;
              cc=cc+c;
          end
      end
  end
  opsi=psi;
 end
%  %% Tecplot
% tt=imax*jmax;
% tec=zeros(tt,3);
%  ii=1;
% for i=1:imax
%     for j=1:jmax
%         tec(ii,1)=x(i);
%         tec(ii,2)=y(j);
%         tec(ii,3)=psi(j,i);
%         tec(ii,4)=U(j,i);
%         tec(ii,5)=V(j,i);
%         tec(ii,6)=P(j,i);
%         ii=ii+1;
%     end
% end
% xlswrite('tecplotR1',tec)