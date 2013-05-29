function [ align_to_ref, align_matrix, nres_ref, sequences ] = convert_alignment_to_matrix( align_lines, ref_idx );

numlines = length( align_lines );
nres_align = length( align_lines{1} );
align_matrix = zeros( nres_align, numlines );

for i = 1:length( align_lines )
  line = align_lines{i};
  count = 0;
  for j = 1:length( line )
    if line(j) ~= '.' & line(j)~='-'
      count = count+1;
      align_matrix(j,i) = count;
    end
  end
  sequences{i} = line( find(align_matrix(:,i)>0) );
end

nres_ref = max( align_matrix( :, ref_idx ) );

gp = find( align_matrix(:,ref_idx) > 0 );
align_to_ref = align_matrix( gp, : );
