function [ bvl ] = calb( xy,bp )

bvl=0.00000;
for i=1:2
    bvl=bvl+xy(i)*bp(i);
end
  


end

