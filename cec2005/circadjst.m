function fit = circadjst(x,v)
if (x<=0)
  x=v+x;
else
    if(x>v)
        x=x-v;
    end
end
fit=x;
