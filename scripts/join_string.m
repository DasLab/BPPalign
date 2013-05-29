function s = join_string( cols, delimiter );
if ~exist( 'delimiter') delimiter = ' '; end;
s = cols{1};
for i = [2:length(cols)]
  s = [s,delimiter,cols{i}];
end


