function [ avl,aij] = cala( xy,ap,ip )

avl=0.0000;
aij=zeros(2,2);

for i=1:2
    for j=1:2
        aij(i,j)=sqrt(ap(i))*sqrt(ap(j))*(1-ip(i,j));
        avl=avl+xy(i)*xy(j)*aij(i,j);
    end
end


end

