function out = more_than(in,prc)
in  = in(:);
in  = in(~isnan(in));
out = in(in>=prctile(in,prc));
end