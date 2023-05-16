clear
clc
lx=241;%x�����ܳ���
ly=81;%y�����ܳ���
Re=40;%��ŵ��
uo=0.05;%���ٶ�
mstep=20000;%���㲽��

%�ϰ�������
obst_r=10;%�ϰ���뾶
obst_x1=2*obst_r+obst_r/2;%Բ��x����
obst_y1=(ly+1)/2;%Բ��y����
nu=uo*2*obst_r/Re;%�˶�ճ��
tau=1/2+3*nu;%�ɳ�ʱ��
omega=1/tau;%�ɳ�Ƶ��
w  = [4/9, 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36]; 
cx = [  0,   1,  0, -1,  0,    1,  -1,  -1,   1]; 
cy = [  0,   0,  1,  0, -1,    1,   1,  -1,  -1]; 
opp = [ 1,   4,  5,  2,  3,    8,   9,   6,   7];%������ķ��� 

in=1;out=lx;%����������������
col=2:(ly-1);
[y,x]=meshgrid(1:ly,1:lx);

%��ʼ��
ux=zeros(1,lx,ly);
uy=zeros(1,lx,ly);
rho=ones(1,lx,ly);
feq=zeros(9,lx,ly);%ƽ�⺯����ʼ
f=zeros(9,lx,ly);%������ʼ
for k=1:9
    feq(k,:,:)=w(k);
    f=feq;
end

%�����ϰ���
obst1=(x-obst_x1).^2 + (y-obst_y1).^2 <= obst_r.^2;%Բ�ı�ʾ����
obst(:,[1,ly]) = 1;    % ����������ұ߽߱�
bbRegion = find(obst1); %���������õ��ϰ��� 
%��ѭ��
tic
for kk=1:mstep
    %������ײ
   for i=1:9 
       cu = 3*(cx(i)*ux+cy(i)*uy); 
       feq(i,:,:)  = rho .* w(i) .* ... 
                       ( 1 + cu + 1/2*(cu.*cu)  - 3/2*(ux.^2+uy.^2) ); 
       f(i,:,:) = f(i,:,:) - omega .* (f(i,:,:)-feq(i,:,:)); 
    end     
    
     %����Ǩ��
    for i=1:9 
       f(i,:,:) = circshift(f(i,:,:), [0,cx(i),cy(i)]); 
    end
     %�������  ��zou-heģ�ͣ�
       ux(1,in,col)=uo;
       uy(1,in,col)=0;
       rho(1,in,col) = 1 ./ (1-ux(1,in,col)) .* (f(5,in,col)+f(1,in,col)+f(3,in,col)+2*(f(4,in,col)+f(7,in,col)+f(8,in,col)));  
       f(2,in,col)=f(4,in,col)+2/3*rho(1,in,col).*ux(1,in,col);
       f(6,in,col)=f(8,in,col)+0.5*(f(5,in,col)-f(3,in,col))+1/6*rho(1,in,col).*ux(1,in,col);
       f(9,in,col)=f(7,in,col)+0.5*(f(3,in,col)-f(5,in,col))+1/6*rho(1,in,col).*ux(1,in,col);
     %��������   ѹ���߽�
    rho(1,out,col) = 1; 
    ux(1,out,col) = -1 + 1  / (rho(1,out,col)) .* (f(5,out,col)+f(1,out,col)+f(3,out,col)+2*(f(2,out,col)+f(6,out,col)+f(9,out,col))); 
    uy(1,out,col) = 0; 
    f(4,out,col) = f(2,out,col) - 2/3*rho(1,out,col).*ux(1,out,col);  
    f(8,out,col) = f(6,out,col) + 1/2*(f(3,out,col)-f(5,out,col))- 1/6*rho(1,out,col).*ux(1,out,col);                                                                            
    f(7,out,col) = f(9,out,col) + 1/2*(f(5,out,col)-f(3,out,col)) - 1/6*rho(1,out,col).*ux(1,out,col);   
   %�ϰ��������ұ߽練�� 
    for i=1:9 
         f(i,bbRegion) = f(opp(i),bbRegion); 
    end 
     %��۱���
    rho = sum(f); 
    ux  = reshape ( (cx * reshape(f,9,lx*ly)), 1,lx,ly) ./rho; 
    uy  = reshape ( (cy * reshape(f,9,lx*ly)), 1,lx,ly) ./rho; 
    %���ӻ�
    if mod(kk,10)==0
        u = reshape(sqrt(ux.^2+uy.^2),lx,ly); 
        u(bbRegion) = nan; 
         imagesc(u'); 
         colorbar('northoutside')
         colormap('jet')
         axis off
         pbaspect([3 1 1])
        axis ([1,241,1,81]); drawnow 
    end 
end
