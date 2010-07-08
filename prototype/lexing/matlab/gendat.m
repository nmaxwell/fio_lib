if(0)
  N = 128;
  f = randn(N,N) + i*randn(N,N);  %mid = [N/4:3*N/4];  f(mid,mid) = 0;
  binstr = sprintf('f_%d.bin', N);
  fid = fopen(binstr,'w');
  string = {'CpxNumMat'};
  serialize(fid, f, string)
  fclose(fid);
  
  N = 512;
  f = randn(N,N) + i*randn(N,N);  %mid = [N/4:3*N/4];  f(mid,mid) = 0;
  binstr = sprintf('f_%d.bin', N);
  fid = fopen(binstr,'w');
  string = {'CpxNumMat'};
  serialize(fid, f, string)
  fclose(fid);
  
  N = 1024;
  f = randn(N,N) + i*randn(N,N);  %mid = [N/4:3*N/4];  f(mid,mid) = 0;
  binstr = sprintf('f_%d.bin', N);
  fid = fopen(binstr,'w');
  string = {'CpxNumMat'};
  serialize(fid, f, string)
  fclose(fid);
  
  N = 2048;
  f = randn(N,N) + i*randn(N,N);  %mid = [N/4:3*N/4];  f(mid,mid) = 0;
  binstr = sprintf('f_%d.bin', N);
  fid = fopen(binstr,'w');
  string = {'CpxNumMat'};
  serialize(fid, f, string)
  fclose(fid);
  
  N = 4096;
  f = randn(N,N) + i*randn(N,N);  %mid = [N/4:3*N/4];  f(mid,mid) = 0;
  binstr = sprintf('f_%d.bin', N);
  fid = fopen(binstr,'w');
  string = {'CpxNumMat'};
  serialize(fid, f, string)
  fclose(fid);
end

if(0)
  EL = 3;
  
  Eall = [5 7 9 11 13 15];
  
  datall = cell(length(Eall), 2);
  for g=1:length(Eall)
    EPS = Eall(g);
    fprintf(1, '%d\n', EPS);
    grid = bfio_grid(EPS);
    [mats,dir] = bfio_prep(EL,EPS);
    datall{g,1} = EPS;
    datall{g,2} = {grid, mats, dir};
  end
  
  binstr = sprintf('bfio.bin');
  fid = fopen(binstr,'w');
  string = {'map' ...
            {'int'} ...
            {'tuple' ...
             {'DblNumVec'} ...
             {'NumVec' ...
              {'CpxNumMat'} ...
             } ...
             {'CpxNumMat'} ...
            } ...
           };
  serialize(fid, datall, string);
  fclose(fid);
  
  binstr = sprintf('bfio.bin');
  fid = fopen(binstr,'r');
  string = {'map' ...
            {'int'} ...
            {'tuple' ...
             {'DblNumVec'} ...
             {'NumVec' ...
              {'CpxNumMat'} ...
             } ...
             {'CpxNumMat'} ...
            } ...
           };
  newall = deserialize(fid, string);
  fclose(fid);
end

