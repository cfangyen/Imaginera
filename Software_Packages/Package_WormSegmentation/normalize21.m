function v=normalize21(v,m)
if nargin<2,
    m=min(v(:));    
else
    assert(numel(m)==1,'size of input must be 1');    
end
m=double(m);
v=double(v);
v=v-m;
v=v./max(v(:));

