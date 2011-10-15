function show_AUG_and_polyU_script( outpath, ref )

% need to align to some template sequence.
align_file = [outpath,'/new_align.txt'];
[ids, align_lines ] = textread( align_file, '%s %s' );

% find the sequence in the alignment that corresponds to desired 'reference', e.g., V. vulnificus for
% adenine riboswitch as ID: AE016796.1
ref_idx = 1;
for i = 1:length( ids ); if ( ~isempty( strfind( ids{i}, ref ) ) ); ref_idx = i; break;end;end;
fprintf( 'Found reference at number: %d\n', ref_idx );

[ align_to_ref, align_matrix, nres_ref, sequences ] = convert_alignment_to_matrix( align_lines, ref_idx );

imagex = ones( length( align_lines ) , nres_ref , 3);

for i = 1:length( align_lines )

  goodres = find( align_to_ref(:,i) > 0 );
  goodres_mapped = align_to_ref(goodres,i);

  UUUUidx = strfind( sequences{i}( goodres_mapped ), 'UUUU' );
  UUUUidx = intersect( UUUUidx, 1:nres_ref);
  imagex( i,UUUUidx,[2 3] ) = 0.0;

  
  AUGidx = strfind( sequences{i}( goodres_mapped ), 'AUG' );
  AUGidx = intersect( AUGidx, 1:nres_ref );
  imagex( i,AUGidx,[1 2] ) = 0.0;

  if ( i == ref_idx)
    sequences{i}
    UUUUidx
    AUGidx
  end

end

clf


plot( -100,-100,'s','markerfacecolor','b','color','b');
hold on
plot( -100,-100,'s','markerfacecolor','r','color','r');

gp = 1:length( align_lines );
%gp = filter_redundant( align_lines, length( align_lines ) );
length(gp)
image( imagex( gp, :, : ) );

axis( [ 0 nres_ref  0 length(gp)] );
make_lines_horizontal( [ref_idx-1, ref_idx], 'k',0.5 );

legend( 'AUG','UUUU',2 );
set( gca,'fontsize',12,'fontweight','bold');h=title( outpath ); set(h,'interpreter','none');
xlabel( ['seqpos [aligned to ',ref,']'] )
ylabel( 'sequence number in rfam alignment' )
