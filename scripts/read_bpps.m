function [nres, all_bpp ] = read_bpps( outpath );
%  [nres, all_bpp ] = read_bpps( outpath );
%  outpath    = path to directory with .bpp files.
%
% (C) R. Das, Stanford University, 2011

N = 2000;

for i = 1:N
  bpp_file = [outpath,'/seq',num2str(i-1),'.bpp' ];
  all_bpp{i} = [];
  if exist( bpp_file )
    all_bpp{i} = load( bpp_file );
    fprintf( 'Reading: %s\n', bpp_file );
    nres(i) = size( all_bpp{i}, 1 );
    end
end
