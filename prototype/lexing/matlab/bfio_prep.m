function [mats,dir] = bfio_prep(EL,EPS)
  
%gs - grid points
%ts - test points

  gs = bfio_grid(EPS);
  
  mats = cell(2,1);
  ts = gs/2;
  mats{1} = bfio_prep_aux(gs,ts);
  ts = gs/2+1/2;
  mats{2} = bfio_prep_aux(gs,ts);
  
  NT = 2^EL;
  ts = [0:NT-1]/NT;
  dir = bfio_prep_aux(gs,ts);
  
