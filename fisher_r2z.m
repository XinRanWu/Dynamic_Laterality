function [z] = fisher_r2z(r)
z = 0.5.* log((1+r)./(1-r)) ;