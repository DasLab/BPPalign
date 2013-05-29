function [ids, align_lines] = read_align_file( align_file );

fid = fopen( align_file );
count = 0;
while  ~feof( fid )
  count = count+1;
  line = fgetl( fid );
  cols = split_string( line );
  align_lines{count} = cols{end};
  ids{count} = join_string( cols(1:end-1) );
end
  
