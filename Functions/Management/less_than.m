function out = less_than(in,prc)
in  = in(:);
in  = in(~isnan(in));
out = in(in<prctile(in,prc));
end