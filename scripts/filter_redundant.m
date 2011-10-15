function gp_filter_redundant = filter_redundant( align_lines, num_lines, CUTOFF );

gp_filter_redundant = [];
 
if ~exist( 'CUTOFF'); CUTOFF = 0.99; end;

line_length = length( align_lines{1} );
for i = 1:num_lines

  match = 0;
  %fprintf( 'Checking redundancy: %d\n', i );
  
    
  for j = (i+1):num_lines
	
    seq_similarity = get_seq_similarity( align_lines{i}, align_lines{j} );
    
    if ( seq_similarity > CUTOFF ) 
      %fprintf( 'Match? %d %d\n', i, j );
      match = 1;
      break;
    end
        
  end
  
  if ~match  & length( align_lines{i} ) == line_length
    gp_filter_redundant =  [ gp_filter_redundant i ];
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function seq_similarity = get_seq_similarity( line1, line2 )

n = 0;
n_same = 0;

seq_similarity = 0.0;

if length( line1 ) ~= length( line2 )
  return;
end

for j = 1:length( line1 )
  if ( line1(j) ~= '.' |...
       line2(j) ~= '.' )
    n = n+1;
    n_same = n_same + ( line1(j) == line2(j) );
  end
end
seq_similarity = n_same/ n ;
     