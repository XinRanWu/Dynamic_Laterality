function [gs] = roi2gs(tc,roi_size)
gs = sum(repmat(roi_size',size(tc,1),1).*tc,2);