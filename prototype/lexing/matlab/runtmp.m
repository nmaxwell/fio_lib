close all;

N = 4096;
fi = 0;
EPS = 9;
EXT = 0;
EL = 3;
SL = log2(N)-EL; 

if(1)
  switch fi
   case 0
    fun = @fun0;
   case 1
    fun = @fun1;
   case 2
    fun = @fun2;
  end
  
  [mats,dir] = bfio_prep(EL,EPS);
  
  if(1)
    f = randn(N,N) + i*randn(N,N);  %mid = [N/4:3*N/4];  f(mid,mid) = 0;
    binstr = sprintf('f_%d.bin', N);
    fid = fopen(binstr,'w');
    string = {'CpxNumMat'};
    serialize(fid, f, string);
  end
  if(0)
    binstr = sprintf('f_%d.bin', N);
    fid = fopen(binstr,'r');
    string = {'CpxNumMat'};
    f = deserialize(fid, string);
  end
	quit()
  t0 = cputime;
  
  %u = bfio_eval(N,SL,EL,EXT,EPS,fun,f,mats,dir); %LEXING
  %profile on;
  bfio_eval;
  %profile report;
  
  te = cputime-t0;
  
  NC = 128;
  t0 = cputime;
  relerr = bfio_check(N,fun,f,u,NC);
  tc = (cputime-t0)*N*N/NC;
  rt = tc/te;
  
  fprintf(1,'N %d\n', N);
  fprintf(1,'EPS %d\n', EPS);
  fprintf(1,'relerr %d\n', relerr);
  fprintf(1,'eval time %d\n',te);
  fprintf(1,'check time %d\n',tc);
  fprintf(1,'ratio %d\n', rt);
end



