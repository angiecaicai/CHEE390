function [ fug ] = calfug( xy,j,P,bp,bvl,avl,Zvl,Bvl,Avl,aij )

sum=0;
for k=1:2
    sum=sum+xy(k)*aij(k);
end
    lnfop=(bp/bvl)*(Zvl-1)-log(Zvl-Bvl)-(Avl/(2*sqrt(2)*Bvl))*((2*sum)/avl-bp/bvl)*log((Zvl+2.414*Bvl)/(Zvl-0.414*Bvl));
    fug=exp(lnfop)*xy(j)*P;
    
end

