function [ xn,yn,vap,liq,fail] = calflash( P,T,Pc,Tc,zn,w,ap,bp,ip )

k=1;
kmax=20;
Ki=zeros(1,2);
for i = 1:2
    Ki(i) = exp(log(Pc(i)/P)+5.37*(1+w(i))*(1-Tc(i)/T));
end
theta=zeros(kmax,1);
vap=0;
liq=0;
fail=0;

for k=1:kmax
    
    P0=zn(1)*(Ki(1)-1)+zn(2)*(Ki(2)-1);
    P1=1-((zn(1)/Ki(1))+(zn(2)/Ki(2)));
    
    if P0*P1<=0
        fun=@(vf) zn(1)*(Ki(1)-1)/(1+vf*(Ki(1)-1))+zn(2)*(Ki(2)-1)/(1+vf*(Ki(2)-1));
        [v]=nr(fun,0.5,0,1,0.0000001,Ki,zn);
        np=2;
    else
        if P0<=0
            v=0;
            np=1;
            %fprintf('Only liquid phase\n');
        else
            if P1<=0
                error('Error in the system\n');
            else
                v=1;
                np=1;
                %fprintf('Only vapor phase\n');
            end
        end
    end
    
    %xi represents liquid fraction of single components, yi vapor fraction
    xi=zn./(1+v.*(Ki-1));
    yi=Ki.*xi;
    
    %normalizing
    xn=xi./sum(xi);
    yn=yi./sum(yi);
    
    %calculate bv, bl, av, al (v is the first element in the matrix)
    bvl=zeros(1,2);
    avl=zeros(1,2);
    
    bvl(1)=calb(yn,bp);
    bvl(2)=calb(xn,bp);
    
    %aij is a 2x2 matrix calculated with interaction coefficient
    aij=zeros(2,2);
    [avl(1),aij]=cala(yn,ap,ip);
    [avl(2),aij]=cala(xn,ap,ip);
    
    %solve the EOS
    Avl=zeros(1,2);
    Bvl=zeros(1,2);
    Zvl=zeros(1,2);
    Vvl=zeros(1,2);
    
    [rv,Avl(1),Bvl(1)]=solvroot(avl(1),bvl(1),P,T);
    [rl,Avl(2),Bvl(2)]=solvroot(avl(2),bvl(2),P,T);
    
    Zvl(1)=max(rv);
    Zvl(2)=min(rl);
    
    %Vvl(1)=Zvl(1)*((8.314*T)/P);
    %Vvl(2)=Zvl(2)*((8.314*T)/P);
    
    %Calculate fugacity [fug1v, fug1l; fug2v, fug2l]
    fug=zeros(2,2);
    fug(1,1)=calfug(yn,1,P,bp(1),bvl(1),avl(1),Zvl(1),Bvl(1),Avl(1),aij(1,:));
    fug(1,2)=calfug(xn,1,P,bp(1),bvl(2),avl(2),Zvl(2),Bvl(2),Avl(2),aij(1,:));
    fug(2,1)=calfug(yn,2,P,bp(2),bvl(1),avl(1),Zvl(1),Bvl(1),Avl(1),aij(2,:));
    fug(2,2)=calfug(xn,2,P,bp(2),bvl(2),avl(2),Zvl(2),Bvl(2),Avl(2),aij(2,:));
    
    %calculate fugacity coefficient
    for i=1:2
        theta(k)=theta(k)+yn(i)*log(fug(i,1)/fug(i,2));
    end
    
    if abs(theta(k))<=1E-7
        if xn==zn
            liq=1;
            break
        end
        %fprintf('Flash converged, mole fraction of carbon dioxide in vapor is %d, in liquid is %d.\n',yn(2),xn(2));
        break
    else
        if k~=1
            if abs(theta(k)-theta(k-1))<=1e-8
                if theta(k)>0 && v==0
                    %fprintf('Feed is in liquid.\n');
                    liq=1;
                    break
                else if theta(k)<0 && v==1
                        %fprintf('Feed is in vapor.\n');
                        vap=1;
                        break
                    end
                end
            end
        end
        
        Ki(1)=Ki(1)*(fug(1,2)/fug(1,1));
        Ki(2)=Ki(2)*(fug(2,2)/fug(2,1));
        
        if np~=2
            Ki=Ki.*exp(theta(k));
        end
        
        if k>=kmax
            %fprintf('Fail\n');
            fail=1;
            break
        end
    end
end


end



