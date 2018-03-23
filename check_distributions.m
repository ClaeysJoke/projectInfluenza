kappa = textread('kappa.txt','%f')
histogram(kappa,100)

%[PI, PT] = textread('truncated_normal.txt','%f %f')
%hist3([PI PT],[20 20])