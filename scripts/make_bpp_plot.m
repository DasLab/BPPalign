function bpp_mean = make_bpp_plot( nres, all_bpp, outpath, ref, structures , box_bounds, whichseq, REDUNDANCY_CUTOFF);
% bpp_mean = make_bpp_plot( nres, all_bpp, outpath, ref, structures , box_bounds, whichseq, REDUNDANCY_CUTOFF);
%
%  nres       = number of residues
%  all_bpp    = cell containing all bpp matrices (use read_bpps.m)
%  outpath    = path to directory with new_align.txt 
%  structures = secondary structures to show on matrix plot (in dot-bracket notation)
%  box_bounds = number of nucleotides upstream and downstream of domain that were pulled in from GenBank
%
%  whichseq   = [optional] which bpps to average over. 
%  REDUNDANCY_CUTOFF = [default 0.99] throw away sequences that are similar at this level
%
% (C) R. Das, Stanford University, 2011

% need to align to some template sequence.
align_file = [outpath,'/new_align.txt'];
[ids, align_lines ] = read_align_file( align_file );

if exist( 'whichseq','var' ) & ~isempty( whichseq); ids = ids( whichseq); align_lines = align_lines( whichseq ); nres = nres( whichseq); all_bpp = all_bpp( whichseq); end;
if ~exist( 'REDUNDANCY_CUTOFF'); REDUNDANCY_CUTOFF = 0.99; end;

% find the sequence in the alignment that corresponds to desired 'reference', e.g., V. vulnificus for
% adenine riboswitch as ID: AE016796.1
ref_idx = 1;
for i = 1:length( ids ); if ( ~isempty( strfind( ids{i}, ref ) ) ); ref_idx = i; break;end;end;
fprintf( 'Found reference at number: %d\n', ref_idx );

[ align_to_ref, align_matrix, nres_ref ] = convert_alignment_to_matrix( align_lines, ref_idx );

% filter out which rows to show.
gp_length_filter = find( nres > 0 & nres < 1000 );
gp = gp_length_filter;

gp = intersect( gp, [1:length(all_bpp)] );

% a surprising number of lines are redundant; RFAM is not filtered well for redundancy.
gp_filter_redundant = filter_redundant( align_lines, max(gp), REDUNDANCY_CUTOFF );
gp = intersect( gp_length_filter, gp_filter_redundant );

n_found_data = length( gp );
fprintf( 'Calculating bpp from %d non-redundant sequences.\n', n_found_data );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bpp_mean = zeros( nres_ref, nres_ref );
for i = gp
  goodres = find( align_to_ref(:,i) > 0 );
  goodres_mapped = align_to_ref(goodres,i);
  if max( goodres_mapped ) <= length( all_bpp{i} )  
    bpp_contribution = all_bpp{i}( goodres_mapped, goodres_mapped );
    bpp_mean( goodres, goodres ) = bpp_mean( goodres, goodres ) + bpp_contribution;
    
    %if length( bpp_contribution ) >192 &  max( bpp_contribution( 98, 192) ) > 0.1; fprintf( 'feature1 --> %d\n', i);end;
    %if length( bpp_contribution )> 312 & max( bpp_contribution( 195, 308:312) ) > 0.1; fprintf( 'feature2 --> %d\n', i);end;
  end
end
bpp_mean = bpp_mean / n_found_data;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clf;
%image( (smooth2d(bpp_mean,2)-0.005) * 1000 );
image( (bpp_mean-0.005) * 1000 );
hold on
plot( [1 nres_ref], [1 nres_ref],'k' );
colormap( 1 - gray(100) );
set(gca,'xgrid','on','ygrid','on');
h=title( outpath ); set(h,'interpreter','none')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist( 'structures' )
  if ~iscell( structures ); structures = { structures } ; end;
  
  colorcode = [1 0 0; 0 0 1];
  for m = 1:length( structures)
    structure = structures{m};
    bps = convert_structure_to_bps( structure );
    for i = 1:size( bps, 1 )
      rectangle( 'Position', [bps(i,1)-0.5 bps(i,2)-0.5 1 1 ],'EdgeColor',colorcode(m,:) );
      %rectangle( 'Position', [bps(i,2)-0.5 bps(i,1)-0.5 1 1 ],'EdgeColor','r' );
    end
  end
end

if exist( 'box_bounds' ) & ~isempty( box_bounds )
  rectangle( 'Position', [ (box_bounds(1) + 0.5)  ( box_bounds(1) + 0.5 ) ...
		    ( nres_ref - box_bounds(2) - box_bounds(1) ) ...
		    ( nres_ref - box_bounds(2) - box_bounds(1) ) ] );
end

hold off
