
NG = 9;



gs = (cos([NG-1:-1:0]/(NG-1)*pi)+1)/2;

  NT = 2^3;
  ts = [0:NT-1]/NT;


  NG = size(gs,2);
  NT = size(ts,2);
  
  tmp = zeros(NG,NT);
  for a=1:NG
    gud = [1:a-1 a+1:NG];
    for b=1:NT
      cur = (ts(b)-gs(gud))./(gs(a)-gs(gud));
      tmp(a,b) = prod(cur);
    end
  end
  res = tmp;
  

res
