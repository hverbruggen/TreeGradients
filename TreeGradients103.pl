#!/usr/bin/perl

#
#	TreeGradients.pl     version 1.02
#
#   heroen.verbruggen@gmail.com
#
# 
#   Version history:
#   1.03  bugfixes for continuous character plotting
#         added linear gradient approximating the Ocean Colors sea surface temperature gradient
#   1.02  multistate drawing added for the three state situation (triangle gradients)
#   1.01  initial release

use warnings;
use strict;
use Bio::TreeIO;
use Sort::Array qw(Sort_Table Discard_Duplicates);
use SVG;
use Graphics::ColorUtils;
use XML::DOM;
use Math::Trig;
use Data::Dumper;

message("\nTreeGradients.pl\n\n");

###################################################################################################################################
#   global variables
###################################################################################################################################

	my (
		# IO files
			$infile,				# input tree file
			$outfile,				# file to which the output is written (in SVG format)

		# tree drawing specifications
			$tree_shape,			# square | circle
			$colorinfo,				# color | bootstrap | value_name  --  name of the value to be used for calculating color information
									#    if 'color', the color field of the node definitions will be used
									#    if 'bootstrap1', 'bootstrap2', 'posterior1', 'posterior2' the support values in the tree will be used and the node definition file becomes non-mandatory
									#    otherwise, the colors will be calculated using the specified variable from the values field of the node definitions
			$variable_type,			# type of variable (discrete, double_discrete, continuous, multistate) -- this is set by program
			$colormode,				# off | plain | gradient  --  method of plotting color on branches
			$gradient_type,			# type of gradient
			$gradient_res,			# resolution of gradient (how many steps)
			$legend_orientation,	# h | v  --  orientation to use for the gradient legend (horizontal | vertical)
			$round_edges,			# 0 | 1  --  specifies whether or not branch edges need to be rounded
			$branch_spacing,		# branch spacing (in mm)
			$branch_width,			# branch width (in mm)
			$tree_hw_ratio,			# height-width ration of the tree
			$node_info,				# value | bootstrap | posterior | nodeID | none  --  determines the values that will be plotted at nodes
			$node_decimals,			# number of decimals to print for node information
			$print_taxon_names,		# 0 | 1  --  print taxon names?
			$l_gap,					# open space in center of circle tree
			
		# data containers
			$tree,					# tree that is read from file
			$taxorder,				# order of taxa for drawing
			$values,				# contains some statistics about the values of the variable that needs to be plotted on the tree
			$svg,					# contains the SVG hash structure
			$multistates,			# sorted array containing the states that the multistate character under consideration can take

		# other
			$branch_radius,			# half of the branch width
			$scale_factor,			# factor with which the branch lengths need to be multiplied to yield the wanted height-width ratio
			$maxtreelength,			# maximum root-to-tip distance (in branch length units)
	);

my $defaults = {
	'tree_shape' => 'square',
	'colormode' => 'gradient',
	'gradient_type' => 'royg',
	'legend_orientation' => 'h',
	'round_edges' => 1,
	'branch_spacing' => 7,
	'branch_width' => 2,
	'tree_hw_ratio' => 1.5,
	'node_info' => 'none',
	'node_decimals' => 2,
	'print_taxon_names' => 1,
	'gradient_res' => 10,
};

my $recognized = {
	'vartype' => {
		'all' => {	'discrete' => 1,
					'continuous' => 1,
					'double_discrete' => 1,
					'multistate' => 1,
		},
	},
	'legend_orientation' => {
		'h' => 1,
		'v' => 1
	},
	'node_info' => {
		'nodeID' => 1,
		'value' => 1,
		'bootstrap' => 1,
		'posterior' => 1,
		'none' => 1,
	},
	'colormode' => {
		'off' => 1,
		'gradient' => 1,
		'plain' => 1,
	},
	'tree_shape' => {
		'square' => 1,
#		'circle' => 1,
	}
};

my $pi = 3.14159265359;

###################################################################################################################################
#   command line parsing
###################################################################################################################################

	unless (	($ARGV[0]) and (substr($ARGV[0],0,1) eq "-") and
				defined($ARGV[1]) and
				($ARGV[2]) and (substr($ARGV[2],0,1) eq "-") and
				defined($ARGV[3]) )
		{ usage(); }
	
	for (my $i=0; $i<scalar(@ARGV); $i+=2) {
		if ($ARGV[$i] eq "-t") { $infile=$ARGV[$i+1]; }
		elsif ($ARGV[$i] eq "-o") { $outfile=$ARGV[$i+1]; }
		elsif ($ARGV[$i] eq "-ts") { $tree_shape=$ARGV[$i+1]; }
		elsif ($ARGV[$i] eq "-vn") { $colorinfo=$ARGV[$i+1]; }
		elsif ($ARGV[$i] eq "-cm") { $colormode=$ARGV[$i+1]; }
		elsif ($ARGV[$i] eq "-gt") { $gradient_type=$ARGV[$i+1]; }
		elsif ($ARGV[$i] eq "-gr") { $gradient_res=$ARGV[$i+1]; }
		elsif ($ARGV[$i] eq "-lo") { $legend_orientation=$ARGV[$i+1]; }
		elsif ($ARGV[$i] eq "-re") { $round_edges=$ARGV[$i+1]; }
		elsif ($ARGV[$i] eq "-bs") { $branch_spacing=$ARGV[$i+1]; }
		elsif ($ARGV[$i] eq "-bw") { $branch_width=$ARGV[$i+1]; }
		elsif ($ARGV[$i] eq "-hw") { $tree_hw_ratio=$ARGV[$i+1]; }
		elsif ($ARGV[$i] eq "-ni") { $node_info=$ARGV[$i+1]; }
		elsif ($ARGV[$i] eq "-nd") { $node_decimals=$ARGV[$i+1]; }
		elsif ($ARGV[$i] eq "-pt") { $print_taxon_names=$ARGV[$i+1]; }
		elsif ($ARGV[$i] eq "-sc") { $l_gap=$ARGV[$i+1]; }
		else { usage(); }
	}
	($infile && $outfile) || usage();
	unless (defined($tree_shape)) {$tree_shape = $defaults->{'tree_shape'};}
	unless ($recognized->{'tree_shape'}->{$tree_shape}) {error("error in specifications: tree shape (-ts switch) cannot take value $tree_shape\n"); usage();}unless (defined($colormode)) {$colormode = $defaults->{'colormode'};}
	unless (defined($colormode)) {$colormode = $defaults->{'colormode'};}
	unless ($recognized->{'colormode'}->{$colormode}) {error("error in specifications: color mode (-cm switch) cannot take value $colormode\n"); usage();}
	unless ($colormode eq 'off') {unless (defined($colorinfo)) {error("\nerror in specifications: if color mode is set to ".$defaults->{'colormode'}.", you have to specify the name of the variable to be plotted on the tree using the -vn option\n"); usage();}}
	unless (defined($gradient_type)) {$gradient_type = $defaults->{'gradient_type'};}
	unless (defined($gradient_res)) {$gradient_res = $defaults->{'gradient_res'};}
	unless (is_valid_gradient_res($gradient_res)) {error("error in specifications: gradient resolution (-gr switch) cannot take value $gradient_res; use an integer value (not 0) instead\n"); usage();}
	unless (defined($legend_orientation)) {$legend_orientation = $defaults->{'legend_orientation'};}
	unless ($recognized->{'legend_orientation'}->{$legend_orientation}) {error("error in specifications: legend orientation (-lo switch) cannot take value $legend_orientation\n"); usage();}
	unless (defined($round_edges)) {$round_edges = $defaults->{'round_edges'};}
	unless (is_binary($round_edges)) {error("error in specifications: round edges (-re switch) cannot take value $round_edges; use a binary value (0 or 1) instead\n"); usage();}
	unless (defined($branch_spacing)) {$branch_spacing = $defaults->{'branch_spacing'};}
	unless (is_numeric($branch_spacing) && ($branch_spacing > 0)) {error("error in specifications: branch spacing (-bs switch) cannot take value $branch_spacing\n"); usage();}
	unless (defined($branch_width)) {$branch_width = $defaults->{'branch_width'};}
	unless (is_numeric($branch_width) && ($branch_width > 0)) {error("error in specifications: branch width (-bw switch) cannot take value $branch_width\n"); usage();}
	unless (defined($tree_hw_ratio)) {$tree_hw_ratio = $defaults->{'tree_hw_ratio'};}
	unless (is_numeric($tree_hw_ratio) && ($tree_hw_ratio > 0)) {error("error in specifications: tree height-width ratio (-hw switch) cannot take value $tree_hw_ratio\n"); usage();}
	unless (defined($node_info)) {$node_info = $defaults->{'node_info'};}
	unless ($recognized->{'node_info'}->{$node_info}) {error("error in specifications: print node information (-ni switch) cannot take value $node_info\n"); usage();}
	unless (defined($node_decimals)) {$node_decimals = $defaults->{'node_decimals'};}
	unless (is_numeric($node_decimals)) {error("error in specifications: node value decimals (-nd switch) cannot take value $node_decimals\n"); usage();}
	unless (defined($print_taxon_names)) {$print_taxon_names = $defaults->{'print_taxon_names'};}
	unless (is_binary($print_taxon_names)) {error("error in specifications: print taxon names (-pt switch) cannot take value $print_taxon_names\n"); usage();}
	$branch_radius = $branch_width / 2;
	
	sub usage {			# prints warning about incorrect command line usage
		message( "\nusage:\n" ); 
		message( "\nmandatory parameters\n" );
		message( "   -t    tree file (in phyloXML format)\n" );
		message( "   -o    name for the output file\n" );
		message( "\noptional parameters -- general\n" );
#		message( "   -ts   tree shape (square|circle -- default: ".$defaults->{'tree_shape'}.")\n" );
		message( "   -vn   name of variable to be plotted on the tree\n" );
		message( "           special values:\n" );
		message( "             bootstrap1    uses the bootstrap values from the tree\n" );
		message( "             bootstrap2    ditto, but values < 50 are given the color for 50\n" );
		message( "             posterior1    uses the posterior probability values from the tree\n" );
		message( "             posterior2    ditto, but values < 0.50 are given the color for 0.50\n" );
		message( "   -pt   print taxon names (0|1 -- default:".$defaults->{'print_taxon_names'}.")\n" );
		message( "   -ni   type of information to print at nodes (default:".$defaults->{'node_info'}.")\n" );
		message( "            value       prints the value of the variable that is plotted on the tree at the node in question\n" );
		message( "            bootstrap   prints bootstrap values\n" );
		message( "            posterior   prints posterior probabilities\n" );
		message( "            nodeID      prints node IDs\n" );
		message( "            none        doesn't print anything at nodes\n" );
		message( "   -nd   number of decimals to print for node information (default:".$defaults->{'node_decimals'}.")\n" );
		message( "   -cm   color mode (off|plain|gradient -- default: ".$defaults->{'colormode'}.")\n" );
		message( "   -gt   gradient type (default: ".$defaults->{'gradient_type'}.")\n" );
		message( "           (1) linear gradients for continuous characters or probabilities of binary characters\n");
		message( "                 bw       black-white\n" );
		message( "                 g[x-y]   shades of gray, from darkness x to darkness y (in the range 0-1)\n" );
		message( "                 royg     red-orange-yellow-green\n" );
		message( "                 rainbow  red-orange-yellow-green-blue-indigo-violet\n" );
		message( "                 bly      blue-yellow\n" );
		message( "                 blw      blue-white\n" );
		message( "                 blg      blue-green\n" );
		message( "                 ocsst    approximation of Ocean Colors sea surface temperature gradient\n" );
		message( "           (2) triangular gradients for probabilities of three-state characters\n");
		message( "                 rgb      red-green-blue\n" );
		message( "                 gby      green-blue-yellow\n" );
		message( "                 gyr      green-yellow-red\n" );
		message( "                 byr      blue-yellow-red\n" );
		message( "   -gr   gradient resolution (default: ".$defaults->{'gradient_res'}.")\n" );
		message( "            number of steps along each gradient (high numbers imply large output files)\n" );
		message( "   -lo   legend orientation (h|v (horizontal or vertical) -- default: ".$defaults->{'legend_orientation'}.")\n" );
		message( "\noptional parameters -- square trees\n" );
		message( "   -re   round edges (1|0 -- default: ".$defaults->{'round_edges'}.")\n" );
		message( "   -bs   branch spacing (in mm -- default: ".$defaults->{'branch_spacing'}.")\n" );
		message( "   -bw   branch width (in mm -- default: ".$defaults->{'branch_width'}.")\n" );
		message( "   -hw   height-width ratio of the tree (default: ".$defaults->{'tree_hw_ratio'}.")\n" );
#		message( "\noptional parameters -- circle trees\n" );
#		message( "   -sc   size of open space in center of tree, as proportion of maximal root to tip path length (default: 0.1)\n" );
		message( "\n" );
		exit;
	}



#################################################################################
#   opening tree file
#################################################################################
{
	message("opening tree file $infile\n");
	my $trees = phyloxml2trees($infile);
	$tree = $trees->[0];
	if ($tree) {
		message("   tree with ".scalar($tree->get_leaf_nodes)." taxa loaded\n");
	} else {
		fatalerror("phyloXML tree parsing failed\n");
	}
	@$taxorder = grep {$_->is_Leaf} @$taxorder;
#	foreach my $tax (@$taxorder) {print "  ",$tax->id,"\n"}; exit;
	$tree->get_root_node->id('root_node');
	$tree = add_node_number_tags($tree);
	if ($tree_shape =~ /^circle$/i) {
		if (defined($l_gap)) {
			$l_gap = maximum_root_to_tip_path_length($tree) * $l_gap;
		} else {
			$l_gap = maximum_root_to_tip_path_length($tree) * 0.1;
		}
	}
}

#################################################################################
#   calculating coordinates
#################################################################################

	message( "calculating coordinates for tree drawing\n");
	if ($tree_shape =~ /^square$/) {
		calculate_coordinates_square_tree();
	} elsif ($tree_shape =~ /^circle$/i) {
		calculate_coordinates_circle_tree();
	}

	
#################################################################################
#   processing node color information
#################################################################################

	message("processing node color information\n");
	if ($colormode eq 'off') {
		message("   colors are switched off -- plotting branches in black\n");
		foreach my $node ($tree->get_nodes) {
			$node->remove_tag('color');
			$node->add_tag_value('color','000000');
		}
	} else {
		message("   color mode set to $colormode\n");
		if ($colorinfo eq 'color') {
			message("   using colors specified in the nodes' color tags\n");
		} elsif ($colorinfo =~ /^bootstrap(\d)$/) {
			message("   converting bootstrap values to colors\n");
			if ($1 == 1) {
				($values->{'matmin'}, $values->{'matmax'}) = (0,100);
			} elsif ($1 == 2) {
				($values->{'matmin'}, $values->{'matmax'}) = (50,100);
			} else {
				fatalerror("double-check command-line input: cannot interpret -cm $colorinfo\n");
			}
			$values->{'ordmagn'} = order_of_magnitude($values->{'matmax'} - $values->{'matmin'});
			$values->{'gradmax'} = gradlimit('upper',$values->{'matmax'},$values->{'ordmagn'});
			$values->{'gradmin'} = gradlimit('lower',$values->{'matmin'},$values->{'ordmagn'});
			message('         gradient range: '.$values->{'gradmin'}.'-'.$values->{'gradmax'}."\n");
			foreach my $node ($tree->get_nodes) {
				if ($node->is_Leaf) {
					if($node->get_tag_values('color')) {$node->remove_tag('color');}
					$node->add_tag_value('color',calc_color(100));
				} else {
					if ($node->has_tag('support')) {
						my $post;
						foreach my $r ($node->get_tag_values('support')) {
							if (defined($r->{'bootstrap'})) {$post = $r->{'bootstrap'};}
						}
						if($node->get_tag_values('color')) {$node->remove_tag('color');}
						if (($colorinfo eq 'bootstrap2') && ($post < 50)) {$post = 50;}
						$node->add_tag_value('color',calc_color($post));
					} elsif ($node->id eq "root_node") {
						my $post = get_support_average_of_descendents($node,'bootstrap');
						if($node->get_tag_values('color')) {$node->remove_tag('color');}
						if (($colorinfo eq 'bootstrap2') && ($post < 50)) {$post = 50;}
						$node->add_tag_value('color',calc_color($post));
					} else {
						fatalerror("node support is not defined for all internal nodes");
					}
				}
			}
		} elsif ($colorinfo =~ /^posterior(\d)$/) {
			message("   converting posterior probabilities to colors\n");
			if ($1 == 1) {
				($values->{'matmin'}, $values->{'matmax'}) = (0,1);
			} elsif ($1 == 2) {
				($values->{'matmin'}, $values->{'matmax'}) = (.5,1);
			} else {
				fatalerror("double-check command-line input: cannot interpret -cm $colorinfo\n");
			}
			$values->{'ordmagn'} = order_of_magnitude($values->{'matmax'} - $values->{'matmin'});
			$values->{'gradmax'} = gradlimit('upper',$values->{'matmax'},$values->{'ordmagn'});
			$values->{'gradmin'} = gradlimit('lower',$values->{'matmin'},$values->{'ordmagn'});
			message('         gradient range: '.$values->{'gradmin'}.'-'.$values->{'gradmax'}."\n");
			foreach my $node ($tree->get_nodes) {
				if ($node->is_Leaf) {
					if($node->get_tag_values('color')) {$node->remove_tag('color');}
					$node->add_tag_value('color',calc_color(1));
				} else {
					if ($node->has_tag('support')) {
						my $post;
						foreach my $r ($node->get_tag_values('support')) {
							if (defined($r->{'posterior'})) {$post = $r->{'posterior'};}
						}
						if($node->get_tag_values('color')) {$node->remove_tag('color');}
						if (($colorinfo eq 'posterior2') && ($post < .50)) {$post = .50;}
						$node->add_tag_value('color',calc_color($post));
					} elsif ($node->id eq "root_node") {
						my $post = get_support_average_of_descendents($node,'posterior');
						if($node->get_tag_values('color')) {$node->remove_tag('color');}
						if (($colorinfo eq 'posterior2') && ($post < .50)) {$post = .50;}
						$node->add_tag_value('color',calc_color($post));
					} else {
						fatalerror("node support is not defined for all internal nodes");
					}
				}
			}
		} else {
			message("   converting $colorinfo values to colors\n");
			unless ($variable_type) {
				fatalerror("Variable type was not set during tree parsing process.\nAre you sure the character \"$colorinfo\" is present in $infile?\n");
			}
			if ($variable_type =~ /^continuous$/i) {
				message("     $colorinfo represents a continuous character\n");
				my $val;
				foreach my $node ($tree->get_nodes) {
					unless ($node->has_tag($colorinfo)) {fatalerror("character $colorinfo is not defined for all nodes\n");}
					my $tags; @$tags = $node->get_all_tags;
					my $c; @$c = $node->get_tag_values($colorinfo);
					push @$val,$c->[0]->{value};
					#print Dumper $c;
				}
				($values->{'matmin'}, $values->{'matmax'}) = minmax($val);
				message('         value range: '.$values->{'matmin'}.'-'.$values->{'matmax'}."\n");
				$values->{'ordmagn'} = order_of_magnitude($values->{'matmax'} - $values->{'matmin'});
				$values->{'gradmax'} = gradlimit('upper',$values->{'matmax'},$values->{'ordmagn'});
				$values->{'gradmin'} = gradlimit('lower',$values->{'matmin'},$values->{'ordmagn'});
				message('         gradient range: '.$values->{'gradmin'}.'-'.$values->{'gradmax'}."\n");
				foreach my $node ($tree->get_nodes) {
					my @c = $node->get_tag_values($colorinfo);
					if($node->get_tag_values('color')) {$node->remove_tag('color');}
					$node->add_tag_value('color',calc_color($c[0]->{value}));
				}
			} elsif ($variable_type =~ /^double_discrete$/i) {
				fatalerror("plotting double_discrete characters is not implemented yet\n");
				message("     $colorinfo represents a set of two discrete characters\n");
				message("     using the default tetrahedron-type gradient\n");
				foreach my $node ($tree->get_nodes) {
					unless ($node->has_tag($colorinfo)) {fatalerror("character $colorinfo is not defined for all nodes\n");}
					my $c; @$c = $node->get_tag_values($colorinfo);
					my $matrix = double_discrete_states_to_matrix($c->[0]->{'state'});
					if($node->get_tag_values('color')) {$node->remove_tag('color');}
					my $color = calc_color_square_discrete($matrix);
					$node->add_tag_value('color',$color);
				}
			} elsif ($variable_type =~ /^discrete$/i) {
				message("     $colorinfo represents a discrete character\n");
				($values->{'matmin'},$values->{'matmax'}) = (0,1);
				$values->{'ordmagn'} = order_of_magnitude($values->{'matmax'} - $values->{'matmin'});
				$values->{'gradmax'} = gradlimit('upper',$values->{'matmax'},$values->{'ordmagn'});
				$values->{'gradmin'} = gradlimit('lower',$values->{'matmin'},$values->{'ordmagn'});
				message('         gradient range: '.$values->{'gradmin'}.'-'.$values->{'gradmax'}."\n");
				foreach my $node ($tree->get_nodes) {
					my $prob;
					unless ($node->has_tag($colorinfo)) {fatalerror("character $colorinfo is not defined for all nodes\n");}
					my @c = $node->get_tag_values($colorinfo);
					foreach my $s (@{$c[0]->{'state'}}) {
						if ($s->{'state'} == 1) {$prob = $s->{'probability'};}
					}
					if($node->get_tag_values('color')) {$node->remove_tag('color');}
					$node->add_tag_value('color',calc_color($prob));
				}
			} elsif ($variable_type =~ /^multistate$/i) {
				message("     $colorinfo represents a multistate character\n");
				foreach my $node ($tree->get_nodes) {
					unless ($node->has_tag($colorinfo)) {fatalerror("character $colorinfo is not defined for all nodes\n");}
					my $c; @$c = $node->get_tag_values($colorinfo);
					if($node->get_tag_values('color')) {$node->remove_tag('color');}
					unless (defined $multistates) {   # set multistates variable
						foreach my $s (@{$c->[0]->{'state'}}) {push @$multistates, $s->{'state'};}
						@$multistates = sort @$multistates;
					}
					my $color = calc_color_three_state_triangle($c->[0]->{'state'});
					$node->add_tag_value('color',$color);
				}
			}
		}
	}
	message("   color processing done\n");

	
#################################################################################
#   making SVG
#################################################################################

	message("generating svg graph and saving to $outfile\n");
	if ($tree_shape =~ /^square$/) {
		output_square_tree();
	} elsif ($tree_shape =~ /^circle$/i) {
		output_circle_tree();
	}
	message("   done\n");

	
	

############################################################################################
#
#  ALL SUBROUTINES FROM THIS POINT ONWARDS
#
############################################################################################
	
###################################################################################################################################
#   general drawing subroutines
###################################################################################################################################


sub generate_triangle_gradient_legend {
	my $svg = shift;
	my $leg = $svg->group('id' => 'legend');
	my $sts = $branch_width; # small triangle size
	my $dimension = 20; # number of small triangles in the longest row of the big triangle
	my $triangles;
	# calculate coordinates for all small triangles and the distances from the base and left edge of the big triangle
	for (my $column = 0; $column < $dimension; ++$column) {
		for (my $row = 0; $row < $dimension - $column; ++$row) {
			my $triangle;
			$triangle->{xbase} = ($column + (0.5 * $row)) * $sts;
			$triangle->{ybase} = $row * $sts * sqrt(0.75);
			$triangle->{xright} = $triangle->{xbase} + $sts;
			$triangle->{yright} = $triangle->{ybase};
			$triangle->{xtop} = $triangle->{xbase} + (0.5 * $sts);
			$triangle->{ytop} = $triangle->{ybase} + (sqrt(0.75) * $sts);
			$triangle->{xcenter} = $triangle->{xtop};
			$triangle->{ycenter} = $triangle->{ybase} + ($sts * 0.3660254);
			$triangle->{distbase} = $triangle->{ycenter};
			my $x3 = $dimension * $sts * 0.5;
			my $y3 = $dimension * $sts * sqrt(0.75);
			my $xleft = $triangle->{ycenter} * $x3 / $y3;
			$triangle->{distleft} = 0.8660254 * abs($xleft - $triangle->{xcenter});
			push @$triangles,$triangle;
			# inverse triangle if $row > 0
			if ($row > 0) {
				my $triangle;
				$triangle->{xbase} = ($column + (0.5 * $row)) * $sts;
				$triangle->{ybase} = $row * $sts * sqrt(0.75);
				$triangle->{xright} = $triangle->{xbase} + $sts;
				$triangle->{yright} = $triangle->{ybase};
				$triangle->{xtop} = $triangle->{xbase} + (0.5 * $sts);
				$triangle->{ytop} = $triangle->{ybase} - (sqrt(0.75) * $sts);
				$triangle->{xcenter} = $triangle->{xtop};
				$triangle->{ycenter} = $triangle->{ybase} - ($sts * 0.3660254);
				$triangle->{distbase} = $triangle->{ycenter};
				my $x3 = $dimension * $sts * 0.5;
				my $y3 = $dimension * $sts * sqrt(0.75);
				my $xleft = $triangle->{ycenter} * $x3 / $y3;
				$triangle->{distleft} = 0.8660254 * abs($xleft - $triangle->{xcenter});
				push @$triangles,$triangle;
			}
		}
	}
	# write the triangles to svg
	my $counter; $counter = 0;
	foreach my $t (@$triangles) {
		++$counter;
		my $d = $sts * $dimension * sqrt(0.75);
		my $p2 = $t->{distbase} / $d;
		my $p3 = $t->{distleft} / $d;
		my $p1 = 1 - $p2 - $p3;
		my $array = [
			{'state' => $multistates->[0],'probability' => $p1},
			{'state' => $multistates->[1],'probability' => $p2},
			{'state' => $multistates->[2],'probability' => $p3},
		];
		my $xv = [$t->{xbase},$t->{xright},$t->{xtop}];
		my $yv = [
			$t->{ybase} + (scalar($tree->get_leaf_nodes) * $branch_spacing) + (5 * $branch_spacing),
			$t->{yright} + (scalar($tree->get_leaf_nodes) * $branch_spacing) + (5 * $branch_spacing),
			$t->{ytop} + (scalar($tree->get_leaf_nodes) * $branch_spacing) + (5 * $branch_spacing),
		];
		my $points = $leg->get_path(x => $xv, y => $yv, -type => 'polygon');
		$leg->polygon(
			%$points,
			id => 'triangle'.$counter,
			'style' => 'fill:'.calc_color_three_state_triangle($array).';stroke:'.calc_color_three_state_triangle($array).';stroke-width:0.05',
		);
	}
	# probabilities near node with lowest state
	$leg->text(
		x		=> -5 * $sts,
		y		=> scalar($tree->get_leaf_nodes) * $branch_spacing + (3.5 * $branch_spacing),
		-cdata	=> 'Pr(state '.$multistates->[0].') = 1',
		'font-size'	=> 4,
	);
	$leg->text(
		x		=> -5 * $sts,
		y		=> scalar($tree->get_leaf_nodes) * $branch_spacing + (4.2 * $branch_spacing),
		-cdata	=> 'Pr(state '.$multistates->[1].') = 0',
		'font-size'	=> 4,
	);
	$leg->text(
		x		=> -5 * $sts,
		y		=> scalar($tree->get_leaf_nodes) * $branch_spacing + (4.9 * $branch_spacing),
		-cdata	=> 'Pr(state '.$multistates->[2].') = 0',
		'font-size'	=> 4,
	);
	# probabilities near node with highest state
	$leg->text(
		x		=> ($dimension - 5) * $sts,
		y		=> scalar($tree->get_leaf_nodes) * $branch_spacing + (3.5 * $branch_spacing),
		-cdata	=> 'Pr(state '.$multistates->[0].') = 0',
		'font-size'	=> 4,
	);
	$leg->text(
		x		=> ($dimension - 5) * $sts,
		y		=> scalar($tree->get_leaf_nodes) * $branch_spacing + (4.2 * $branch_spacing),
		-cdata	=> 'Pr(state '.$multistates->[1].') = 0',
		'font-size'	=> 4,
	);
	$leg->text(
		x		=> ($dimension - 5) * $sts,
		y		=> scalar($tree->get_leaf_nodes) * $branch_spacing + (4.9 * $branch_spacing),
		-cdata	=> 'Pr(state '.$multistates->[2].') = 1',
		'font-size'	=> 4,
	);
	# probabilities near node with intermediate state
	$leg->text(
		x		=> (($dimension / 2) - 5) * $sts,
		y		=> scalar($tree->get_leaf_nodes) * $branch_spacing + ((3.5 + 2.1) * $branch_spacing) + ($dimension * $sts * sqrt(0.75)),
		-cdata	=> 'Pr(state '.$multistates->[0].') = 0',
		'font-size'	=> 4,
	);
	$leg->text(
		x		=> (($dimension / 2) - 5) * $sts,
		y		=> scalar($tree->get_leaf_nodes) * $branch_spacing + ((4.2 + 2.1) * $branch_spacing) + ($dimension * $sts * sqrt(0.75)),
		-cdata	=> 'Pr(state '.$multistates->[1].') = 1',
		'font-size'	=> 4,
	);
	$leg->text(
		x		=> (($dimension / 2) - 5) * $sts,
		y		=> scalar($tree->get_leaf_nodes) * $branch_spacing + ((4.9 + 2.1) * $branch_spacing) + ($dimension * $sts * sqrt(0.75)),
		-cdata	=> 'Pr(state '.$multistates->[2].') = 0',
		'font-size'	=> 4,
	);
}

sub generate_gradient_legend {
	my $svg = shift;
	my $orientation = shift;
	my $leg = $svg->group('id' => 'legend');
	if ($orientation eq 'v') {
		# draw the gradient
		my $gradient = $leg->gradient(
			'-type' => "linear",
			'id' => 'leggrad',
			'x1' => "0%",
			'y1' => "0%",
			'x2' => "0%",
			'y2' => "100%",
		);
		my $step = 100 / $gradient_res;
		for (my $i=0; $i<=100; $i += $step) {
			my $color = calc_color($values->{'gradmin'} + ( ($i/100) * ($values->{'gradmax'} - $values->{'gradmin'}) ));
			$gradient->stop(
				'offset' => $i."%",
				'style' => 'stop-color:#'.$color.';stop-opacity:1',
			);
		}
		$leg->rect(
			x		=>	-2 * $branch_spacing,
			y		=>	scalar($tree->get_leaf_nodes) * $branch_spacing - (5 * $branch_spacing),
			width	=>	$branch_spacing,
			height	=>	$branch_spacing * 5,
			'style' => 'fill:url(#leggrad)',
			'stroke-width'	=>	0,
		);
		if ($colorinfo =~ /((bootstrap|posterior))\d/) {
			my $val = $1; if ($val =~ /bootstrap/) {$val = 'Bootstrap value';} else {$val = 'Posterior probability';}
			$leg->text(
				x		=> -2 * $branch_spacing,
				y		=> scalar($tree->get_leaf_nodes) * $branch_spacing - (5.5 * $branch_spacing),
				-cdata	=> $val,
				'font-size'	=> 4,
			);
		} else {
			if ($variable_type =~ /^discrete$/i) {
				$leg->text(
					x		=> -2 * $branch_spacing,
					y		=> scalar($tree->get_leaf_nodes) * $branch_spacing - (5.5 * $branch_spacing),
					-cdata	=> 'Pr (character state 1)',
					'font-size'	=> 4,
				);
			} elsif ($variable_type =~ /^continuous$/i) {
				$leg->text(
					x		=> -2 * $branch_spacing,
					y		=> scalar($tree->get_leaf_nodes) * $branch_spacing - (5.5 * $branch_spacing),
					-cdata	=> 'Character state',
					'font-size'	=> 4,
				);
			}
		}
		# draw some labels along the gradient
		for (my $i=0; $i <=100; $i+=20) {
			my $x = - $branch_spacing;
			my $y = scalar($tree->get_leaf_nodes) * $branch_spacing - (($i+1) * $branch_spacing / 20);
			$leg->rect(
				x		=>	$x,
				y		=>	$y + 0.75 * $branch_spacing / 20,
				width	=>  1,
				height	=>	0.25 * $branch_spacing / 10,
				'stroke-width'	=>	0,
			);
			$leg->text(
				x		=>	$x + 2,
				y		=>	$y + 1,
				-cdata	=>	($values->{'gradmin'} + ( ($i/100) * ($values->{'gradmax'} - $values->{'gradmin'}) )),
				'font-size'	=>	2,
			);
		}
	} elsif ($orientation eq 'h') {
		# draw the gradient
		my $gradient = $leg->gradient(
			'-type' => "linear",
			'id' => 'leggrad',
			'x1' => "0%",
			'y1' => "0%",
			'x2' => "100%",
			'y2' => "0%",
		);
		my $step = 100 / $gradient_res;
		for (my $i=0; $i<=100; $i += $step) {
			my $color = calc_color($values->{'gradmin'} + ( ($i/100) * ($values->{'gradmax'} - $values->{'gradmin'}) ));
			$gradient->stop(
				'offset' => $i."%",
				'style' => 'stop-color:#'.$color.';stop-opacity:1',
			);
		}
		$leg->rect(
			x		=>	0,
			y		=>	(scalar($tree->get_leaf_nodes) + 3) * $branch_spacing,
			width	=>	$branch_spacing * 5,
			height	=>	$branch_spacing,
			'style' => 'fill:url(#leggrad)',
			'stroke-width'	=>	0,
		);
		if ($colorinfo =~ /((bootstrap|posterior))\d/) {
			my $val = $1; if ($val =~ /bootstrap/) {$val = 'Bootstrap value';} else {$val = 'Posterior probability';}
			$leg->text(
				x		=> 0,
				y		=> (scalar($tree->get_leaf_nodes) + 2.2) * $branch_spacing,
				-cdata	=> $val,
				'font-size'	=> 4,
			);
		} else {
			if ($variable_type =~ /^discrete$/i) {
				$leg->text(
					x		=> 0,
					y		=> (scalar($tree->get_leaf_nodes) + 2.2) * $branch_spacing,
					-cdata	=> 'Pr (character state 1)',
					'font-size'	=> 4,
				);
			} elsif ($variable_type =~ /^continuous$/i) {
				$leg->text(
					x		=> 0,
					y		=> (scalar($tree->get_leaf_nodes) + 2.2) * $branch_spacing,
					-cdata	=> 'Character state',
					'font-size'	=> 4,
				);
			}
		}
		# draw some labels along the gradient
		for (my $i=0; $i <=100; $i+=20) {
			my $x = ($i * $branch_spacing / 20) - .0125 * $branch_spacing;
			my $y = (scalar($tree->get_leaf_nodes) + 3) * $branch_spacing;
			$leg->rect(
				x		=>	$x,
				y		=>	$y - 1,
				width	=>  0.25 * $branch_spacing / 10,
				height	=>	1,
				'stroke-width'	=>	0,
			);
			my $label = ($values->{'gradmin'} + ( ($i/100) * ($values->{'gradmax'} - $values->{'gradmin'}) ));
			$leg->text(
				x		=>	$x - (length($label)/2),
				y		=>	$y - 2,
				-cdata	=>	$label,
				'font-size'	=>	2,
			);
		}
	}
}

###################################################################################################################################
#   circle tree subroutines
###################################################################################################################################

sub calculate_coordinates_circle_tree {
	# calculating scale factor for branch lengths
		{
		my $max; $max = 0;
		my $node1 = $tree->get_root_node;
		foreach my $term ($tree->get_leaf_nodes) {
			my $nodes; @$nodes = ($node1,$term);
			my $dist = $tree->distance('-nodes' => $nodes);
			if ($dist > $max) {$max = $dist;}
		}
		my $maxdraw = scalar($tree->get_leaf_nodes) * $branch_spacing / $tree_hw_ratio;
		$scale_factor = $maxdraw / $max;
		$maxtreelength = $max;
		}
	# calculating distances
	#    (from the root node outwards to the terminals)
		{
		my $start_node = $tree->get_root_node;
		calc_distance_circle($start_node,0);
		message( "   distances calculated for all nodes\n");
		}
	# angles of terminal nodes are pretty straightforward
	{
		my $counter; $counter = 0;
		my $n = scalar($tree->get_leaf_nodes);
		foreach my $node ($tree->get_leaf_nodes) {
			++$counter;
			$node->add_tag_value('_angle', (($pi/2)-(($counter-1)*(11*$pi)/(6*$n))) );
		}
		message( "   angles calculated for terminal nodes\n" );
	# angles of internal nodes need to be resolved from the terminal nodes inwards
		my $nr_descendents_of_nodes;
		foreach my $node (grep {!($_->is_Leaf)} $tree->get_nodes) {
			my $term = calc_nr_terminals($node);
			push @{$nr_descendents_of_nodes->{$term}},$node;
		}
		my $keys; @$keys = keys(%$nr_descendents_of_nodes);
		@$keys = Discard_Duplicates('sorting' => 'ascending', 'data' => $keys);
		foreach my $level (@$keys) {
			foreach my $node (@{$nr_descendents_of_nodes->{$level}}) {
				my $descs; @$descs = $node->each_Descendent;
				my $int = calc_internal_angle_circle($descs);
				$node->add_tag_value('_angle',$int);
			}
		}
		message( "   angles calculated for internal nodes\n" );
	}
	# convert everything to cartesian coordinates
		foreach my $node ($tree->get_nodes) {
			my @angle = $node->get_tag_values('_angle');
			my @dist = $node->get_tag_values('_dist');
			my ($x,$y) = polar2cartesian($dist[0],$angle[0]);
			$node->add_tag_value('_x',$x);
			$node->add_tag_value('_y',$y);
		}
		message( "   polar coordinates converted to cartesian coordinates\n" );
}

sub calc_internal_angle_circle {
	my $descs = shift;
	my $out; $out = 0;
	foreach my $desc (@$descs) {
		my @y = $desc->get_tag_values('_angle');
		$out += $y[0];
	}
	$out /= scalar(@$descs);
	return $out;
}

sub calc_distance_circle {
	my $node = shift;
	my $offset = shift;
	my $dist; $dist = $offset;
	if ($offset == 0) {$dist += $l_gap;}
	if ($node->branch_length) {$dist += $node->branch_length * $scale_factor;}
	$node->add_tag_value('_dist',$dist);
	if ($node->each_Descendent) {
		foreach my $desc ($node->each_Descendent) {
			calc_distance_circle($desc,$dist);
		}
	}
}

sub output_circle_tree {
	$svg = SVG->new('width' => 1000, 'height' => 1000);
	my $y = $svg->group('id' => 'branches');
	my $gr = $svg->defs('id' => 'gradients');
	my $start_node = $tree->get_root_node;
	generate_svg_tree_circle($svg,$y,$gr,$start_node);
	add_node_labels_to_svg_circle($svg,$tree);
	generate_scalebar_circle($svg,$maxtreelength);
	unless ($colormode =~ /^off$/i) { generate_gradient_legend($svg,$legend_orientation); }
	open FH,">$outfile";
	print FH $svg->xmlify;
	close FH;
}

sub generate_svg_tree_circle {
	my $svg = shift;
	my $y = shift;
	my $gr = shift;
	my $node = shift;
	my @nr = $node->get_tag_values('node_nr');
	my @xp = $node->get_tag_values('_x');
	my @yp = $node->get_tag_values('_y');
	my @ap = $node->get_tag_values('_angle');
	my @dp = $node->get_tag_values('_dist');
	my @ch = $node->get_tag_values('color');
	my @cp = $node->get_tag_values('values');
	my $vp = $cp[0]->{$colorinfo};
	if ($colormode eq 'gradient') {
		# continue here
	} else {
=add this circle for complete thing
		$y->circle('cx' => $xp[0], 'cy' => $yp[0], 'r' => $branch_radius, 'style' => 'fill:'.$ch[0]);    # plot circle at node
=cut
		if ($node->each_Descendent) {
			foreach my $desc ($node->each_Descendent) {
				my @ad = $desc->get_tag_values('_angle');
				my @dd = $desc->get_tag_values('_dist');
				my @xd = $desc->get_tag_values('_x');
				my @yd = $desc->get_tag_values('_y');
				my @chd = $desc->get_tag_values('color');
				my ($ix,$iy) = polar2cartesian($dp[0],$ad[0]);
				my $angle = $ad[0] * 180 / $pi;
=complete thing -- no lines but filled elements
				# plot circle at edge of curve
					$y->circle('cx' => $ix, 'cy' => $iy, 'r' => $branch_radius, 'style' => 'fill:'.$chd[0]);
				# plot line towards daughter node
					{
					my ($x1,$y1) = polar2cartesian($dp[0],$ad[0]-atan($branch_radius/$dp[0]));
					my ($x4,$y4) = polar2cartesian($dp[0],$ad[0]+atan($branch_radius/$dp[0]));
					my ($x2,$y2) = polar2cartesian($dd[0],$ad[0]-atan($branch_radius/$dd[0]));
					my ($x3,$y3) = polar2cartesian($dd[0],$ad[0]+atan($branch_radius/$dd[0]));
					my ($dx2,$dy2) = ($x2-$x1,$y2-$y1);
					my ($dx3,$dy3) = ($x3-$x2,$y3-$y2);
					my ($dx4,$dy4) = ($x4-$x3,$y4-$y3);
					$y->path('d' => "M$x1,$y1 l$dx2,$dy2 l$dx3,$dy3 l$dx4,$dy4 z", 'style' => 'fill:'.$chd[0], 'stroke-width' => "0");
					}
				# plot curve
					{
					my $inner_radius = $dp[0]-$branch_radius;
					my $outer_radius = $dp[0]+$branch_radius;
					my ($x1,$y1) = polar2cartesian($inner_radius,$ad[0]);
					my ($x4,$y4) = polar2cartesian($inner_radius,$ap[0]);
					my ($x2,$y2) = polar2cartesian($outer_radius,$ad[0]);
					my ($x3,$y3) = polar2cartesian($outer_radius,$ap[0]);
					my ($dx2,$dy2) = ($x2-$x1,$y2-$y1);
					my ($dx3,$dy3) = ($x3-$x2,$y3-$y2);
					my ($dx4,$dy4) = ($x4-$x3,$y4-$y3);
					my ($dx1,$dy1) = ($x1-$x4,$y1-$y4);
					my ($fl1,$fl2,$fl3); $fl1 = 0; if ($ad[0] - $ap[0] > 0) {$fl2 = 1; $fl3 = 0;} else {$fl2 = 0; $fl3 = 1}
					$y->path(
						'd' => 
							"M$x2,$y2 ".
							"a".$outer_radius.",".$outer_radius." 0 $fl1,$fl2 $dx3,$dy3 ".
							"l$dx4,$dy4 ".
							"a".$inner_radius.",".$inner_radius." 0 $fl1,$fl3 $dx1,$dy1 ".
							"z",
						'style' => 'fill:'.$chd[0],
						'stroke-width' => "1"
					);
					}
=cut
#=this is an alternative using lines instead of closed paths
				# do calculations for curve
				my ($rx,$ry) = ($dp[0],$dp[0]);
				my ($xdiff,$ydiff) = ($xp[0] - $ix, $yp[0] - $iy);
				my ($fl1,$fl2); $fl1 = 0;
				if ($ad[0] - $ap[0] > 0) {$fl2 = 1;} else {$fl2 = 0;}
				# draw curve
				$y->path('d' => "M$ix,$iy a$rx,$ry 0 $fl1,$fl2 $xdiff,$ydiff", 'fill' => "none", 'stroke' => "black", 'stroke-width' => $branch_radius*2);
				# draw line
				($ix,$iy) = polar2cartesian($dp[0],$ad[0]);
				($xdiff,$ydiff) = ($xd[0] - $ix, $yd[0] - $iy);
				$y->path('d' => "M$ix,$iy l$xdiff,$ydiff", 'fill' => "none", 'stroke' => "black", 'stroke-width' => $branch_radius*2);
#=cut
				# draw descendents
				generate_svg_tree_circle($svg,$y,$gr,$desc);
			}
		}
	}
}

sub add_node_labels_to_svg_circle {
	my $svg = shift;
	my $tree = shift;
	my $z = $svg->group('id' => 'node_labels');
	my $y = $svg->group('id' => 'taxon_labels');
	foreach my $node ($tree->get_leaf_nodes) {
		# add code to insert character values of terminal nodes
		my @a = $node->get_tag_values('_angle');
		my @d = $node->get_tag_values('_dist');
		my $amp = maximum_root_to_tip_path_length($tree) + $l_gap;
		my $ntaxa = scalar($tree->get_leaf_nodes);
		my $angle = -180 * $a[0] / $pi;
		my ($x1,$y1) = polar2cartesian($d[0]+($amp/100),$a[0]-(($pi/2)/$ntaxa));
		if ($print_taxon_names) {
			$y->text('x' => $x1, 'y' => $y1, 'transform' => 'rotate('.$angle.",$x1,$y1".')', 'font-size' => 2 * $amp / $ntaxa )->cdata($node->id);
		}
	}
	# adding labels for internal nodes
	unless ($node_info =~ /^none$/i) {
		foreach my $node (grep {!($_->is_Leaf)} $tree->get_nodes) {
			my @a = $node->get_tag_values('_angle'); my $angle=$a[0];
			my @d = $node->get_tag_values('_dist'); my $dist=$d[0];
			my $val;
			if (($node_info =~ /^nodeID$/i) and $node->id) {
				$val = $node->id;
				if ($val =~ /^[E\.\d\-]+$/i) {$val = $b = sprintf("%.".$node_decimals."f", $val);}
			} elsif ($node_info =~ /^value$/i) {
				my @c = $node->get_tag_values('values');
				$val = $c[0]->{$colorinfo};
				$val = $b = sprintf("%.".$node_decimals."f", $val);
			}
			my ($x,$y) = polar2cartesian($dist - $branch_width - length($val), $angle - atan($branch_width/$dist));
			$z->text('x' => $x, 'y' => $y - $branch_width, 'font-size' => $branch_spacing * 0.4)->cdata($val);
		}
	}
}

sub generate_scalebar_circle {
}

###################################################################################################################################
#   square tree subroutines
###################################################################################################################################

sub calculate_coordinates_square_tree {
	# calculating scale factor for branch lengths
		{
		my $max; $max = 0;
		my $node1 = $tree->get_root_node;
		foreach my $term ($tree->get_leaf_nodes) {
			my $nodes; @$nodes = ($node1,$term);
			my $dist = $tree->distance('-nodes' => $nodes);
			if ($dist > $max) {$max = $dist;}
		}
		my $maxdraw = scalar($tree->get_leaf_nodes) * $branch_spacing / $tree_hw_ratio;
		$scale_factor = $maxdraw / $max;
		$maxtreelength = $max;
		}
	# calculating x coordinates
	#    (from the root node outwards to the terminals)
		{
		my $start_node = $tree->get_root_node;
		calc_x_coord_square($start_node,0);
		message( "   x coordinates processed for all nodes\n");
		}
	# y coordinates of terminal nodes are pretty straightforward
	{
		my $counter; $counter = 0;
		foreach my $node ($tree->get_leaf_nodes) {
			++$counter;
			$node->add_tag_value('y',($branch_spacing * $counter));
		}
	# y coordinates of internal nodes need to be resolved from the terminal nodes inwards
		my $nr_descendents_of_nodes;
		foreach my $node ($tree->get_nodes) {
			unless ($node->is_Leaf) {
				my $term = calc_nr_terminals($node);
				push @{$nr_descendents_of_nodes->{$term}},$node;
			}
		}
		my $keys; @$keys = keys(%$nr_descendents_of_nodes);
		@$keys = Discard_Duplicates('sorting' => 'ascending', 'data' => $keys);
		foreach my $level (@$keys) {
			foreach my $node (@{$nr_descendents_of_nodes->{$level}}) {
				my $descs; @$descs = $node->each_Descendent;
				$node->add_tag_value('y',calc_avg_y_square($descs));
			}
		}
		message( "   y coordinates processed for all nodes\n" );
	}
}

sub output_square_tree {
	$svg = SVG->new('width' => (scalar($tree->get_leaf_nodes) * $branch_spacing * 1.5 / $tree_hw_ratio), 'height' => scalar($tree->get_leaf_nodes) * $branch_spacing * 1.5);
	my $y = $svg->group('id' => 'branches');
	my $gr = $svg->defs('id' => 'gradients');
	my $start_node = $tree->get_root_node;
	generate_svg_tree_square($svg,$y,$gr,$start_node);
	add_node_labels_to_svg_square($svg,$tree);
	generate_scalebar_square($svg,$maxtreelength);
	unless ($colormode =~ /^off$/i) {
		if ($colorinfo =~ /(bootstrap|posterior)\d/) {
			generate_gradient_legend($svg,$legend_orientation);
		} else {
			if ($variable_type =~ /^continuous$/i) {
				generate_gradient_legend($svg,$legend_orientation);
			} elsif ($variable_type =~ /^double_discrete$/i) {
				message("  gradient legend for double_discrete not implemented yet\n");
			} elsif ($variable_type =~ /^discrete$/i) {
				generate_gradient_legend($svg,$legend_orientation);
			} elsif ($variable_type =~ /^multistate$/i) {
				generate_triangle_gradient_legend($svg);
			}
		}
	}
	open FH,">$outfile";
	print FH $svg->xmlify;
	close FH;
}

sub calc_avg_y_square {
	my $descs = shift;
	my $out; $out = 0;
	foreach my $desc (@$descs) {
		my @y = $desc->get_tag_values('y');
		$out += $y[0];
	}
	$out /= scalar(@$descs);
	return $out;
}

sub calc_x_coord_square {
	my $node = shift;
	my $offset = shift;
	my $x; $x = $offset;
	if ($node->branch_length) {$x += $node->branch_length * $scale_factor;}
	$node->add_tag_value('x',$x);
	if ($node->each_Descendent) {
		foreach my $desc ($node->each_Descendent) {
			calc_x_coord_square($desc,$x);
		}
	}
}

sub generate_svg_tree_square {
	my $svg = shift;
	my $y = shift;
	my $gr = shift;
	my $node = shift;
	my @nr = $node->get_tag_values('node_nr');
	my @xp = $node->get_tag_values('x');
	my @yp = $node->get_tag_values('y');
	my @ch = $node->get_tag_values('color');
	my @cp = $node->get_tag_values($colorinfo);
	if ($colormode eq 'gradient') {
		if ($colorinfo =~ /((bootstrap|posterior))\d/) {
			if ($round_edges) {
				$y->circle('cx' => $xp[0], 'cy' => $yp[0], 'r' => $branch_radius, 'style' => 'fill:#'.$ch[0]);
			}
			my $support_type = $1;
			my $vp;
			if ($node->id && ($node->id eq 'root_node')) {
				$vp = get_support_average_of_descendents($node,$support_type);
			} elsif ($node->is_Leaf) {
				if ($support_type eq 'bootstrap') {
					$vp = 100;
				} elsif($support_type eq 'posterior') {
					$vp = 1;
				}
			} else {
				my $a; @$a = $node->get_tag_values('support');
				foreach my $b (@$a) {
					if (defined($b->{$support_type})) {
						$vp = $b->{$support_type};
					}
				}
			}
			if (($colorinfo eq 'bootstrap2') && ($vp < 50)) {$vp = 50;}
			if (($colorinfo eq 'posterior2') && ($vp < .50)) {$vp = .50;}
			if ($node->each_Descendent) {
				foreach my $desc ($node->each_Descendent) {
					my @nrd = $desc->get_tag_values('node_nr');
					my @xd = $desc->get_tag_values('x');
					my @yd = $desc->get_tag_values('y');
					my @chd = $desc->get_tag_values('color');
					my @cd = $desc->get_tag_values($colorinfo);
					my $vd;
					if ($desc->is_Leaf) {
						if ($support_type eq 'bootstrap') {
							$vd = 100;
						} elsif($support_type eq 'posterior') {
							$vd = 1;
						}
					} else {
						my $a; @$a = $desc->get_tag_values('support');
						foreach my $b (@$a) {
							if (defined($b->{$support_type})) {
								$vd = $b->{$support_type};
							}
						}
					}
					if (($colorinfo eq 'bootstrap2') && ($vd < 50)) {$vd = 50;}
					if (($colorinfo eq 'posterior2') && ($vd < .50)) {$vd = .50;}
					my $gradientname = 'gradient_between_nodes_'.$nr[0].'_and_'.$nrd[0];
					my $gradient = $gr->gradient(
						'-type' => "linear",
						'id' => $gradientname,
						'x1' => "0%",
						'y1' => "0%",
						'x2' => "100%",
						'y2' => "0%",
					);
					my $step = 100 / $gradient_res;
					for (my $p=0; $p <= 100; $p += $step) {
						my $val = (($vd-$vp) * $p / 100) + $vp;
						$gradient->stop(
							'offset' => $p."%",
							'style' => 'stop-color:#'.calc_color($val).';stop-opacity:1',
						);
					}
					if ($round_edges) {
						$y->circle('cx' => $xp[0], 'cy' => $yd[0], 'r' => $branch_radius, 'style' => 'fill:#'.$ch[0]);
						$y->rect('x' => $xp[0] - $branch_radius, 'y' => min2($yp[0],$yd[0]), 'width' => $branch_width, 'height' => abs($yd[0] - $yp[0]), 'style' => 'fill:#'.$ch[0]);
						$y->rect('x' => $xp[0], 'y' => $yd[0] - $branch_radius, 'width' => $xd[0] - $xp[0], 'height' => $branch_width, 'style' => 'fill:url(#'.$gradientname.')');
					}
					else {
						$y->rect('x' => $xp[0] - $branch_radius, 'y' => min2($yp[0],$yd[0]), 'width' => $branch_width, 'height' => abs($yd[0] - $yp[0]), 'style' => 'fill:#'.$ch[0]);
						$y->rect('x' => $xp[0] - $branch_radius, 'y' => $yd[0] - $branch_radius, 'width' => $xd[0] - $xp[0] + $branch_radius, 'height' => $branch_width, 'style' => 'fill:url(#'.$gradientname.')');
					}
					generate_svg_tree_square($svg,$y,$gr,$desc);
				}
			}
		} else {
			if ($variable_type =~ /^continuous$/i) {
				my $vp = $cp[0]->{'value'};
				if ($round_edges) {
					$y->circle('cx' => $xp[0], 'cy' => $yp[0], 'r' => $branch_radius, 'style' => 'fill:#'.$ch[0]);
				}
				if ($node->each_Descendent) {
					foreach my $desc ($node->each_Descendent) {
						my @nrd = $desc->get_tag_values('node_nr');
						my @xd = $desc->get_tag_values('x');
						my @yd = $desc->get_tag_values('y');
						my @chd = $desc->get_tag_values('color');
						my @cd = $desc->get_tag_values($colorinfo);
						my $vd = $cd[0]->{'value'};
						my $gradientname = 'gradient_between_nodes_'.$nr[0].'_and_'.$nrd[0];
						my $gradient = $gr->gradient(
							'-type' => "linear",
							'id' => $gradientname,
							'x1' => "0%",
							'y1' => "0%",
							'x2' => "100%",
							'y2' => "0%",
						);
						my $step = 100 / $gradient_res;
						for (my $p=0; $p <= 100; $p += $step) {
							my $val = (($vd-$vp) * $p / 100) + $vp;
							$gradient->stop(
								'offset' => $p."%",
								'style' => 'stop-color:#'.calc_color($val).';stop-opacity:1',
							);
						}
						if ($round_edges) {
							$y->circle('cx' => $xp[0], 'cy' => $yd[0], 'r' => $branch_radius, 'style' => 'fill:#'.$ch[0]);
							$y->rect('x' => $xp[0] - $branch_radius, 'y' => min2($yp[0],$yd[0]), 'width' => $branch_width, 'height' => abs($yd[0] - $yp[0]), 'style' => 'fill:#'.$ch[0]);
							$y->rect('x' => $xp[0], 'y' => $yd[0] - $branch_radius, 'width' => $xd[0] - $xp[0], 'height' => $branch_width, 'style' => 'fill:url(#'.$gradientname.')');
						}
						else {
							$y->rect('x' => $xp[0] - $branch_radius, 'y' => min2($yp[0],$yd[0]), 'width' => $branch_width, 'height' => abs($yd[0] - $yp[0]), 'style' => 'fill:#'.$ch[0]);
							$y->rect('x' => $xp[0] - $branch_radius, 'y' => $yd[0] - $branch_radius, 'width' => $xd[0] - $xp[0] + $branch_radius, 'height' => $branch_width, 'style' => 'fill:url(#'.$gradientname.')');
						}
						generate_svg_tree_square($svg,$y,$gr,$desc);
					}
				}
			} elsif ($variable_type =~ /^double_discrete$/i) {
				my $vp = $cp[0]->{'state'};
				if ($round_edges) {
					$y->circle('cx' => $xp[0], 'cy' => $yp[0], 'r' => $branch_radius, 'style' => 'fill:#'.$ch[0]);
				}
				if ($node->each_Descendent) {
					foreach my $desc ($node->each_Descendent) {
						my @nrd = $desc->get_tag_values('node_nr');
						my @xd = $desc->get_tag_values('x');
						my @yd = $desc->get_tag_values('y');
						my @chd = $desc->get_tag_values('color');
						my @cd = $desc->get_tag_values($colorinfo);
						my $vd = $cd[0]->{'state'};
						my $gradientname = 'gradient_between_nodes_'.$nr[0].'_and_'.$nrd[0];
						my $gradient = $gr->gradient(
							'-type' => "linear",
							'id' => $gradientname,
							'x1' => "0%",
							'y1' => "0%",
							'x2' => "100%",
							'y2' => "0%",
						);
						my $step = 100 / $gradient_res;
						for (my $p=0; $p <= 100; $p += $step) {
							my $pm = double_discrete_states_to_matrix($vp);
							my $dm = double_discrete_states_to_matrix($vd);
							my $matrix = [
								((($dm->[0]-$pm->[0]) * $p / 100) + $pm->[0]),
								((($dm->[1]-$pm->[1]) * $p / 100) + $pm->[1]),
								((($dm->[2]-$pm->[2]) * $p / 100) + $pm->[2]),
								((($dm->[3]-$pm->[3]) * $p / 100) + $pm->[3]),
							];
							$gradient->stop(
								'offset' => $p."%",
								'style' => 'stop-color:'.calc_color_square_discrete($matrix).';stop-opacity:1',
							);
						}
						if ($round_edges) {
							$y->circle('cx' => $xp[0], 'cy' => $yd[0], 'r' => $branch_radius, 'style' => 'fill:#'.$ch[0]);
							$y->rect('x' => $xp[0] - $branch_radius, 'y' => min2($yp[0],$yd[0]), 'width' => $branch_width, 'height' => abs($yd[0] - $yp[0]), 'style' => 'fill:#'.$ch[0]);
							$y->rect('x' => $xp[0], 'y' => $yd[0] - $branch_radius, 'width' => $xd[0] - $xp[0], 'height' => $branch_width, 'style' => 'fill:url(#'.$gradientname.')');
						}
						else {
							$y->rect('x' => $xp[0] - $branch_radius, 'y' => min2($yp[0],$yd[0]), 'width' => $branch_width, 'height' => abs($yd[0] - $yp[0]), 'style' => 'fill:#'.$ch[0]);
							$y->rect('x' => $xp[0] - $branch_radius, 'y' => $yd[0] - $branch_radius, 'width' => $xd[0] - $xp[0] + $branch_radius, 'height' => $branch_width, 'style' => 'fill:url(#'.$gradientname.')');
						}
						generate_svg_tree_square($svg,$y,$gr,$desc);
					}
				}
			} elsif ($variable_type =~ /^discrete$/i) {
				my $vp = $cp[0]->{'state'};
				if ($round_edges) {
					$y->circle('cx' => $xp[0], 'cy' => $yp[0], 'r' => $branch_radius, 'style' => 'fill:#'.$ch[0]);
				}
				foreach my $s (@$vp) {if ($s->{'state'} == 1) {$vp = $s->{'probability'};}}
				if ($node->each_Descendent) {
					foreach my $desc ($node->each_Descendent) {
						my @nrd = $desc->get_tag_values('node_nr');
						my @xd = $desc->get_tag_values('x');
						my @yd = $desc->get_tag_values('y');
						my @chd = $desc->get_tag_values('color');
						my @cd = $desc->get_tag_values($colorinfo);
						my $vd = $cd[0]->{'state'};
						foreach my $s (@$vd) {if ($s->{'state'} == 1) {$vd = $s->{'probability'};}}
						my $gradientname = 'gradient_between_nodes_'.$nr[0].'_and_'.$nrd[0];
						my $gradient = $gr->gradient(
							'-type' => "linear",
							'id' => $gradientname,
							'x1' => "0%",
							'y1' => "0%",
							'x2' => "100%",
							'y2' => "0%",
						);
						my $step = 100 / $gradient_res;
						for (my $p=0; $p <= 100; $p += $step) {
							my $val = (($vd-$vp) * $p / 100) + $vp;
							$gradient->stop(
								'offset' => $p."%",
								'style' => 'stop-color:#'.calc_color($val).';stop-opacity:1',
							);
						}
						if ($round_edges) {
							$y->circle('cx' => $xp[0], 'cy' => $yd[0], 'r' => $branch_radius, 'style' => 'fill:#'.$ch[0]);
							$y->rect('x' => $xp[0] - $branch_radius, 'y' => min2($yp[0],$yd[0]), 'width' => $branch_width, 'height' => abs($yd[0] - $yp[0]), 'style' => 'fill:#'.$ch[0]);
							$y->rect('x' => $xp[0], 'y' => $yd[0] - $branch_radius, 'width' => $xd[0] - $xp[0], 'height' => $branch_width, 'style' => 'fill:url(#'.$gradientname.')');
						}
						else {
							$y->rect('x' => $xp[0] - $branch_radius, 'y' => min2($yp[0],$yd[0]), 'width' => $branch_width, 'height' => abs($yd[0] - $yp[0]), 'style' => 'fill:#'.$ch[0]);
							$y->rect('x' => $xp[0] - $branch_radius, 'y' => $yd[0] - $branch_radius, 'width' => $xd[0] - $xp[0] + $branch_radius, 'height' => $branch_width, 'style' => 'fill:url(#'.$gradientname.')');
						}
						generate_svg_tree_square($svg,$y,$gr,$desc);
					}
				}
			} elsif ($variable_type =~ /^multistate$/i) {
				my $vp = $cp[0]->{'state'};
				if ($round_edges) {
					$y->circle('cx' => $xp[0], 'cy' => $yp[0], 'r' => $branch_radius, 'style' => 'fill:'.$ch[0]);
				}
				if ($node->each_Descendent) {
					foreach my $desc ($node->each_Descendent) {
						my @nrd = $desc->get_tag_values('node_nr');
						my @xd = $desc->get_tag_values('x');
						my @yd = $desc->get_tag_values('y');
						my @chd = $desc->get_tag_values('color');
						my @cd = $desc->get_tag_values($colorinfo);
						my $vd = $cd[0]->{'state'};
						my $gradientname = 'gradient_between_nodes_'.$nr[0].'_and_'.$nrd[0];
						my $gradient = $gr->gradient(
							'-type' => "linear",
							'id' => $gradientname,
							'x1' => "0%",
							'y1' => "0%",
							'x2' => "100%",
							'y2' => "0%",
						);
						my $step = 100 / $gradient_res;
						for (my $p=0; $p <= 100; $p += $step) {
							my $array;
							foreach my $pe (@$vp) {
								foreach my $de (@$vd) {
									if ($pe->{'state'} eq $de->{'state'}) {
										push @$array,{'state' => $pe->{'state'}, 'probability' => ( ((100-$p) * 0.01 * $pe->{'probability'}) + ($p * 0.01 * $de->{'probability'} ))};
									}
								}
							}
							$gradient->stop(
								'offset' => $p."%",
								'style' => 'stop-color:'.calc_color_three_state_triangle($array).';stop-opacity:1',
							);
						}
						if ($round_edges) {
							$y->circle('cx' => $xp[0], 'cy' => $yd[0], 'r' => $branch_radius, 'style' => 'fill:'.$ch[0]);
							$y->rect('x' => $xp[0] - $branch_radius, 'y' => min2($yp[0],$yd[0]), 'width' => $branch_width, 'height' => abs($yd[0] - $yp[0]), 'style' => 'fill:'.$ch[0]);
							$y->rect('x' => $xp[0], 'y' => $yd[0] - $branch_radius, 'width' => $xd[0] - $xp[0], 'height' => $branch_width, 'style' => 'fill:url(#'.$gradientname.')');
						}
						else {
							$y->rect('x' => $xp[0] - $branch_radius, 'y' => min2($yp[0],$yd[0]), 'width' => $branch_width, 'height' => abs($yd[0] - $yp[0]), 'style' => 'fill:'.$ch[0]);
							$y->rect('x' => $xp[0] - $branch_radius, 'y' => $yd[0] - $branch_radius, 'width' => $xd[0] - $xp[0] + $branch_radius, 'height' => $branch_width, 'style' => 'fill:url(#'.$gradientname.')');
						}
						generate_svg_tree_square($svg,$y,$gr,$desc);
					}
				}
			}
		}
	} else {
		if ($round_edges) {
			$y->circle('cx' => $xp[0], 'cy' => $yp[0], 'r' => $branch_radius, 'style' => 'fill:#'.$ch[0]);
		}
		if ($node->each_Descendent) {
			foreach my $desc ($node->each_Descendent) {
				my @xd = $desc->get_tag_values('x');
				my @yd = $desc->get_tag_values('y');
				my @chd = $desc->get_tag_values('color');
				if ($round_edges) {
					$y->circle('cx' => $xp[0], 'cy' => $yd[0], 'r' => $branch_radius, 'style' => 'fill:#'.$chd[0]);
					$y->rect('x' => $xp[0] - $branch_radius, 'y' => min2($yp[0],$yd[0]), 'width' => $branch_width, 'height' => abs($yd[0] - $yp[0]), 'style' => 'fill:#'.$chd[0]);
					$y->rect('x' => $xp[0], 'y' => $yd[0] - $branch_radius, 'width' => $xd[0] - $xp[0], 'height' => $branch_width, 'style' => 'fill:#'.$chd[0]);
				}
				else {
					$y->rect('x' => $xp[0] - $branch_radius, 'y' => min2($yp[0],$yd[0]), 'width' => $branch_width, 'height' => abs($yd[0] - $yp[0]), 'style' => 'fill:#'.$chd[0]);
					$y->rect('x' => $xp[0] - $branch_radius, 'y' => $yd[0] - $branch_radius, 'width' => $xd[0] - $xp[0] + $branch_radius, 'height' => $branch_width, 'style' => 'fill:#'.$chd[0]);
				}
				generate_svg_tree_square($svg,$y,$gr,$desc);
			}
		}
	}
}

sub add_node_labels_to_svg_square {
	my $svg = shift;
	my $tree = shift;
    my $y = $svg->group('id' => 'taxon_labels');
    my $z = $svg->group('id' => 'node_labels');
	foreach my $node ($tree->get_nodes) {
		my @xp = $node->get_tag_values('x');
		my @yp = $node->get_tag_values('y');
		if ($node->is_Leaf) {
			if ($node_info =~ /^value$/i) {
				if ((defined $colorinfo) && !($colorinfo =~ /^(bootstrap|posterior)\d$/)) {
					my $val;
					if ($variable_type =~ /^continuous$/i) {
						my $c; @$c = $node->get_tag_values($colorinfo); $c = $c->[0];
						$val = $c->{'value'};
						$val = '('.sprintf("%.".$node_decimals."f", $val).')';
						$y->text('x' => $xp[0] + (1.5 * $branch_width), 'y' => $yp[0] + ($branch_spacing * 0.2) - .5, 'font-size' => $branch_spacing * 0.4)->cdata($val);
						#message("  node labels for continuous characters not implemented yet\n");
					} elsif ($variable_type =~ /^double_discrete$/i) {
						message("  node labels for double_discrete characters not implemented yet\n");
					} elsif ($variable_type =~ /^discrete$/i) {
						my $c; @$c = $node->get_tag_values($colorinfo); $c = $c->[0];
						foreach my $state (@{$c->{'state'}}) {
							if ($state->{'state'} eq 1) {$val = $state->{'probability'};}
						}
						$val = '('.sprintf("%.".$node_decimals."f", $val).')';
						$y->text('x' => $xp[0] + (1.5 * $branch_width), 'y' => $yp[0] + ($branch_spacing * 0.2) - .5, 'font-size' => $branch_spacing * 0.4)->cdata($val);
					} elsif ($variable_type =~ /^multistate$/i) {
						my $c; @$c = $node->get_tag_values($colorinfo); $c = $c->[0];
						foreach my $s (@{$c->{'state'}}) {if ($s->{'probability'} == 1) {$val = '(P'.$s->{'state'}.' = 1.00)';}}
						$y->text('x' => $xp[0] + (1.5 * $branch_width), 'y' => $yp[0] + ($branch_spacing * 0.2) - .5, 'font-size' => $branch_spacing * 0.4)->cdata($val);
					}
					if ($print_taxon_names) {
						$y->text('x' => $xp[0] + (2 * $branch_width) + (1.5 * length($val)), 'y' => $yp[0] + ($branch_spacing * 0.2), 'font-size' => $branch_spacing * 0.7)->cdata($node->id);
					}
				} else {
					unless ($colorinfo =~ /^(bootstrap|posterior)\d$/) {
						error("cannot plot character values at external node because no character was specified\n");
					}
				}
			} else {
				if ($print_taxon_names) {
					$y->text('x' => $xp[0] + (2 * $branch_width), 'y' => $yp[0] + ($branch_spacing * 0.2), 'font-size' => $branch_spacing * 0.7)->cdata($node->id);
				}
			}
		} else {
			unless ($node_info =~ /^none$/i) {
				my $val;
				if ($node_info =~ /^nodeID$/i) {
					if (defined($node->{'_id'})) {
						unless ($node->id eq 'root_node') {
							$val = $node->id;
							if ($val =~ /^[E\.\d\-]+$/i) {$val = $b = sprintf("%.".$node_decimals."f", $val);}
						}
					} else {
						error("cannot plot nodeID of internal node because it is not defined\n");
					}
				} elsif ($node_info =~ /^bootstrap$/i) {
					my $st = 'bootstrap';
					unless ($node->id && ($node->id eq 'root_node')) {
						unless ($node->get_tag_values('support')) {error("cannot plot requested support value because it is not defined\n");}
						my $c; @$c = $node->get_tag_values('support'); $c = $c->[0];
						if (defined($c->{$st})) {$val = $c->{$st};} else {error("cannot plot requested support value because it is not defined\n");}
					}
				} elsif ($node_info =~ /^posterior$/i) {
					my $st = 'posterior';
					unless ($node->id && ($node->id eq 'root_node')) {
						unless ($node->get_tag_values('support')) {error("cannot plot requested support value because it is not defined\n");}
						my $c; @$c = $node->get_tag_values('support'); $c = $c->[0];
						if (defined($c->{$st})) {$val = $c->{$st};} else {error("cannot plot requested support value because it is not defined\n");}
					}
				} elsif ($node_info =~ /^value$/i) {
					if (defined $colorinfo){
						if ($colorinfo =~ /((bootstrap|posterior))\d/) {
							my $st = $1;
							unless ($node->id && ($node->id eq 'root_node')) {
								unless ($node->get_tag_values('support')) {error("cannot plot requested support value because it is not defined\n");}
								my $c; @$c = $node->get_tag_values('support'); $c = $c->[0];
								if (defined($c->{$st})) {$val = $c->{$st};} else {error("cannot plot requested support value because it is not defined\n");}
							}
						} else {
							if ($variable_type =~ /^continuous$/i) {
								my $c; @$c = $node->get_tag_values($colorinfo); $c = $c->[0];
								$val = $c->{'value'};
								$val = sprintf("%.".$node_decimals."f", $val);
								#$y->text('x' => $xp[0] + (1.5 * $branch_width), 'y' => $yp[0] + ($branch_spacing * 0.2) - .5, 'font-size' => $branch_spacing * 0.4)->cdata($val);
								#message("  node labels for continuous characters not implemented yet\n");
							} elsif ($variable_type =~ /^double_discrete$/i) {
								message("  node labels for double_discrete characters not implemented yet\n");
							} elsif ($variable_type =~ /^discrete$/i) {
								my $c; @$c = $node->get_tag_values($colorinfo); $c = $c->[0];
								foreach my $state (@{$c->{'state'}}) {
									if ($state->{'state'} eq 1) {$val = $state->{'probability'};}
								}
								$val = sprintf("%.".$node_decimals."f", $val);
							} elsif ($variable_type =~ /^multistate$/i) {
								my $c; @$c = $node->get_tag_values($colorinfo); $c = $c->[0];
								my $datahash;
								foreach my $s (@{$c->{'state'}}) {
									$datahash->{'P'.$s->{'state'}} = $s->{'probability'};
								}
								$val = "";
								foreach my $key (sort keys %$datahash) {
									$val .= $key.' = '.sprintf("%.".$node_decimals."f", $datahash->{$key}).'\n';
								}
							}
						}
					} else {
						error("cannot plot character values at internal node because no character was specified\n");
					}
				}
				if (defined $val) {
					if ($val =~ /\\n/) {
						my $array; @$array = split /\\n/,$val;
						my $i; $i=0;
						foreach my $entry (@$array) { unless ($entry =~ /^\s*$/) {
							$z->text('x' => $xp[0] - $branch_width - length($entry) - 1, 'y' => $yp[0] - ((1 - (1.2*$i)) * $branch_width), 'font-size' => $branch_spacing * 0.4)->cdata($entry);
							++$i;
						}}
					} else {
						$z->text('x' => $xp[0] - $branch_width - length($val) - 1, 'y' => $yp[0] - $branch_width, 'font-size' => $branch_spacing * 0.4)->cdata($val);
					}
				}
			}
		}
	}
}

sub generate_scalebar_square {
	my $svg = shift;
	my $maxtreelength = shift;
	my $sc = $svg->group();
	my $magn = order_of_magnitude($maxtreelength);
	my $val = 10 ** $magn;
	my $len = $val * $scale_factor;
	$sc->line(
		x1	=>	0,
		y1	=>	(scalar($tree->get_leaf_nodes)+1) * $branch_spacing,
		x2	=>	$len,
		y2	=>	(scalar($tree->get_leaf_nodes)+1) * $branch_spacing,
		style => "stroke-width: 0.5; stroke: black",
	);
	$sc->line(
		x1	=>	0,
		y1	=>	(scalar($tree->get_leaf_nodes)+1) * $branch_spacing - 1,
		x2	=>	0,
		y2	=>	(scalar($tree->get_leaf_nodes)+1) * $branch_spacing + 1,
		style => "stroke-width: 0.5; stroke: black",
	);
	$sc->line(
		x1	=>	$len,
		y1	=>	(scalar($tree->get_leaf_nodes)+1) * $branch_spacing - 1,
		x2	=>	$len,
		y2	=>	(scalar($tree->get_leaf_nodes)+1) * $branch_spacing + 1,
		style => "stroke-width: 0.5; stroke: black",
	);
	$sc->text(
		x		=>	$len / 10,
		y		=>	(scalar($tree->get_leaf_nodes)+1) * $branch_spacing - 1,
		id		=>	'scale',
		'font-size' => $branch_spacing * 0.3,
		-cdata	=>	$val,
	);
}

###################################################################################################################################
#   color processing
###################################################################################################################################

sub calc_color_three_state_triangle {
	my $in = shift;
	my ($r,$g,$b);
	if ($gradient_type =~ /^rgb$/i) {
		my $data; foreach my $state (@$in) {$data->{$state->{'state'}} = 255 * $state->{'probability'};}
		my $keys; @$keys = sort keys %$data;
		($r,$g,$b) = ($data->{$keys->[0]},$data->{$keys->[1]},$data->{$keys->[2]});
	} elsif ($gradient_type =~ /^gby$/i) {
		my $data; foreach my $state (@$in) {$data->{$state->{'state'}} = $state->{'probability'};}
		my $keys; @$keys = sort keys %$data;
		my ($p1,$p2) = ($data->{$keys->[0]},$data->{$keys->[1]});
		$r = 242 - 180 * $p1 - 221 * $p2;
		$g = 252 -  88 * $p1 - 3.2 * $p2;
		$b =  12 + 203 * $p1 + 7.2 * $p2;
	} elsif ($gradient_type =~ /^gyr$/i) {
		my $data; foreach my $state (@$in) {$data->{$state->{'state'}} = $state->{'probability'};}
		my $keys; @$keys = sort keys %$data;
		my ($p1,$p2) = ($data->{$keys->[0]},$data->{$keys->[1]});
		$r = 263.5 -   6.8 * $p1 - 235.2 * $p2;
		$g = 268.6 - 229.6 * $p1 -   1.2 * $p2;
		$b =  45.2 +     0 * $p1 +   0.4 * $p2;
	} elsif ($gradient_type =~ /^byr$/i) {
		my $data; foreach my $state (@$in) {$data->{$state->{'state'}} = $state->{'probability'};}
		my $keys; @$keys = sort keys %$data;
		my ($p1,$p2) = ($data->{$keys->[0]},$data->{$keys->[1]});
		$r = 268.3 - 228.4 * $p1 -   0.0 * $p2;
		$g = 277.3 - 131.6 * $p1 - 230.8 * $p2;
		$b =  30.1 + 228.8 * $p1 +  14.4 * $p2;
	} else {fatalerror("error: gradient type is not recognized as a triangular gradient: $gradient_type\n");}
	return "rgb\($r\,$g\,$b\)";
}

sub calc_color_square_discrete {
	my $matrix = shift;
	my ($b,$c,$d,$a) = @$matrix;
	#my $y = $d - ( sin(3.141592/6) * ($a + $b + $c) );
	my ($side1,$side2,$side3);
	if ($a+$b) {$side1 = $a / ($a+$b);} else {$side1 = 0;}
	if ($b+$c) {$side2 = $b / ($b+$c);} else {$side2 = 0;}
	if ($c+$a) {$side3 = $c / ($c+$a);} else {$side3 = 0;}	
	my $hue1 = int($side1 * 120);# / 360;
	my $hue2 = int($side2 * 120 + 120);# / 360;
	my $hue3 = int($side3 * 120 + 240);# / 360;
	#message("    hues: $hue1 $hue2 $hue3\n");
	my ($r1,$g1,$b1) = hsv2rgb($hue1,1,1);
	my ($r2,$g2,$b2) = hsv2rgb($hue2,1,1);
	my ($r3,$g3,$b3) = hsv2rgb($hue3,1,1);
	#message("    rgb:  1: $r1 $g1 $b1 -- 2: $r2 $g2 $b2 -- 3: $r3 $g3 $b3\n");
	my ($r,$g,$k) = (($r1+$r2+$r3)/3,($g1+$g2+$g3)/3,($b1+$b2+$b3)/3);
	#message("    avg: $r $g $k\n");
	my ($h,$s,$v) = rgb2hsv($r,$g,$k);
	#message("    hsv: $h $s $v\n");
	$v = $v * $d;
	#message("    hsv2: $h $s $v\n");
	($r,$g,$k) = hsv2rgb($h,$s,$v);
	#message("    returning: $r $g $k\n");
	return "rgb\($r\,$g\,$k\)";
}

sub calc_color {
	my $val = shift;
	my $out;
	my ($min,$max) = ($values->{'gradmin'},$values->{'gradmax'});
	if ($gradient_type =~ /^bw$/i) {
		my $r = (255 / ($min - $max));
		my $a = (255 * $max) / ($max - $min);
		my $y = int(($r * $val) + $a);
		$out = "rgb\($y\,$y\,$y\)";
	} elsif ($gradient_type =~ /^royg$/i) {
		my $c = (120 / ($min - $max));
		my $a = (120 * $max) / ($max - $min);
		my $y = ($c * $val) + $a;
		my ($r,$g,$b) = hsv2rgb($y,1,1);
		$out = "rgb\($r\,$g\,$b\)";
	} elsif ($gradient_type =~ /^blg$/i) {
		my $a = (129 - 244) / ($max - $min);
		my $k = 244 - ($a * $min);
		my $y = ($a * $val) + $k;
		my ($r,$g,$b) = hsv2rgb($y,1,1);
		$out = "rgb\($r\,$g\,$b\)";
	} elsif ($gradient_type =~ /^bly$/i) {
		my $r = int((255 / ($max-$min)) * ($val - $min));
		my $g = int(((155 / ($max-$min)) * ($val - $min)) + 100);
		my $b = int((255 / ($min-$max)) * ($val - $max));
		$out = "rgb\($r\,$g\,$b\)";
	} elsif ($gradient_type =~ /^blw$/i) {
		my $r = int((((255-20)*$val)/($min-$max)) + 20 - (((255-20)*$max)/($min-$max)));
		my $g = int((((255-0)*$val)/($min-$max)) + 0 - (((255-0)*$max)/($min-$max)));
		my $b = int((((255-215)*$val)/($min-$max)) + 215 - (((255-215)*$max)/($min-$max)));
		$out = "rgb\($r\,$g\,$b\)";
	} elsif ($gradient_type =~ /^rainbow$/i) {
		my $c = (300 / ($min - $max));
		my $a = (300 * $max) / ($max - $min);
		my $y = ($c * $val) + $a;
		my ($r,$g,$b) = hsv2rgb($y,1,1);
		$out = "rgb\($r\,$g\,$b\)";
	} elsif ($gradient_type =~ /^g\[(.+?)\-(.+?)\]$/) {
		my ($low,$high) = ($1,$2);
		unless (($low <= 1) and ($low >= 0) and ($high <= 1) and ($high >= 0)) {fatalerror("Values given for gradient out of 0-1 range\n");}
		my $a = 255 * (1-$low);
		my $b = 255 * (1-$high);
		my $y = int($a + ((($b - $a)/($max - $min))*($val - $min)));
		$out = "rgb\($y\,$y\,$y\)";
	} elsif ($gradient_type =~ /^ocsst$/i) {  # ocean colors sea surface temperature approximation
		my ($r,$g,$b);
		my $t = $val;
		if ($t<1) {$r=-46*$t + 54} elsif ($t<20) {$r=8} elsif ($t<27) {$r=-0.0956*$t*$t*$t*$t + 8.2961*$t*$t*$t - 265.5*$t*$t + 3739*$t - 19639} elsif ($t<31) {$r=224} else {$r=0.8333*$t*$t*$t*$t - 108.76*$t*$t*$t + 5312.7*$t*$t -115145*$t + 934640}
		if ($t<3) {$g=8} elsif ($t<13) {$g=-0.7331*$t*$t + 35.011*$t - 93.933} elsif ($t<26) {$g=-0.076*$t*$t*$t*$t + 6.04567*$t*$t*$t - 175.04*$t*$t + 2178.2*$t - 9617.5} else {$g=-0.0758*$t*$t*$t*$t + 10.04*$t*$t*$t - 492.3*$t*$t + 10570*$t - 83619}
		if ($t<5) {$b=1.0648*$t*$t*$t - 10.04*$t*$t - 3.4299*$t + 238.34} elsif ($t<13) {$b = -0.1126*$t*$t + 18.593*$t + 12.978} elsif ($t<22) {$b=0.037*$t*$t*$t*$t - 3.1606*$t*$t*$t + 98.528*$t*$t -1353.1*$t + 7069.9} else {$b=8}
		$r = int($r); $g = int($g); $b = int($b);
		$out = "rgb\($r\,$g\,$b\)";
	} else {fatalerror("error: gradient type is not recognized as a linear gradient: $gradient_type\n");}
	$out =~ /^rgb\((\d+)\,(\d+)\,(\d+)\)$/;
	my ($h,$s,$v) =  rgb2hsv($1,$2,$3);
	$out = substr(hsv2rgb($h,$s,$v),1);
	return $out;
}


###################################################################################################################################
#   tree processing
###################################################################################################################################

sub add_node_number_tags {
	my $tree = shift;
	my $counter; $counter = 0;
	foreach my $node ($tree->get_nodes) {
		++$counter;
		$node->add_tag_value('node_nr',$counter);
	}
	return $tree;
}

sub maximum_root_to_tip_path_length {
	my $tree = shift;
	my $max; $max = 0;
	foreach my $node ($tree->get_leaf_nodes) {
		my $dist = dist_to_root($node);
		if ($dist > $max) {$max = $dist;}
	}
	return $max;
}

sub dist_to_root {
	my $node = shift;
	my $dist; $dist = 0;
	while (!(defined($node->id)) || !($node->id eq 'root_node')) {
		$dist += $node->branch_length;
		$node = $node->ancestor;
	}
	return $dist;
}

sub calc_nr_terminals {
	my $node = shift;
	return scalar(grep {$_->is_Leaf} $node->get_all_Descendents);
}

sub search_node {
	my $tree = shift;
	my $searchtype = shift;
	my $query = shift;
	if ($searchtype eq 'id') {
		foreach my $node ($tree->get_nodes) {
			if (($node->id) and ($node->id eq $query)) {return $node;}
		}
	} elsif ($searchtype eq 'desc') {
		my $node = search_node($tree,'id',$query->[0]);
		my $cont; $cont = 1;
		$node = $node->ancestor;
		while ($cont) {
			if (arraymatch($query,descendentlist($node))) {
				return $node;
			} elsif ($node->id eq 'root_node') {
				die "cannot find node with the given descendents\n";
			} else {
				$node = $node->ancestor;
			}
		}
	} else {
		die "currently the subroutine search_node does not support searching anything else than IDs\n";
	}
	use Data::Dumper; print Dumper $searchtype; print Dumper $query;
	die "no node found with $searchtype corresponding to $query\n";
}

sub descendentlist {
	my $in = shift;
	my $out;
	foreach my $node ($in->get_all_Descendents) {
		if ($node->is_Leaf) {push @$out,$node->id;}
	}
	return $out;
}


###################################################################################################################################
#   general calculations
###################################################################################################################################

sub double_discrete_states_to_matrix {
	my $states = shift;
	my $matrix;
	foreach my $s (@$states) {
		if ($s->{'state'} eq '00') {$matrix->[0] = $s->{'probability'};}
		if ($s->{'state'} eq '10') {$matrix->[1] = $s->{'probability'};}
		if ($s->{'state'} eq '01') {$matrix->[2] = $s->{'probability'};}
		if ($s->{'state'} eq '11') {$matrix->[3] = $s->{'probability'};}
	}
	return $matrix;
}

sub polar2cartesian {
	my $dist = shift;
	my $angle = shift;
	my $x = $dist * cos $angle;
	my $y = - $dist * sin $angle;
	return ($x,$y);
}

sub min2 {
	my $i = shift;
	my $j = shift;
	if ($i<$j) {return $i;} else {return $j;}
}

sub arraymatch {
	my $a1 = shift;
	my $a2 = shift;
	if (scalar(@$a1) != scalar(@$a2)) {return 0;}
	my $str;
	foreach my $e (@$a1) {$str->{$e} = 1;}
	foreach my $e (@$a2) {$str->{$e} = 1;}
	if (scalar(keys %$str) == scalar(@$a1)) {return 1;} else {return 0;}
}

sub minmax {
	my $matrix = shift;
	my $min; $min = $matrix->[0];
	my $max; $max = $matrix->[0];
	foreach my $element (@$matrix) {
		if ($element =~ /^[0-9\.]+$/) {
			if ($element > $max) {$max = $element;}
			if ($element < $min) {$min = $element;}
		}
	}
	return ($min,$max);
}

sub gradlimit {
	my $ul = shift;
	my $val = shift;
	my $magn = shift;
	if ($ul eq 'upper') {
		if ($val >= 0) {
			my $counter; $counter=1; while ($counter*(10**($magn-1)) < $val) {++$counter;}
			return ($counter*(10**($magn-1)));
		} else {
			my $counter; $counter=-1; while ($counter*(10**($magn-1)) >= $val) {--$counter;}
			return (($counter+1)*(10**($magn-1)));
		}
	} elsif ($ul eq 'lower') {
		if ($val >= 0) {
			my $counter; $counter=0; while ($counter*(10**($magn-1)) <= $val) {++$counter;}
			return (($counter-1)*(10**($magn-1)));
		} else {
			my $counter; $counter=0; while ($counter*(10**($magn-1)) > $val) {--$counter;}
			return ($counter*(10**($magn-1)));
		}
	}
}

sub order_of_magnitude {
	my $t1 = sprintf('%e',shift);
	$t1 =~ /.+e\+*0*(.+)/i;
	return $1;
}

sub all_keys_processed {
	my $keys = shift;
	my $processed = shift;
	foreach my $key (@$keys) {
		unless ($processed->{$key}) {return 0;}
	}
	return 1;
}

sub is_numeric {
	my $in = shift;
	unless ($in =~ /^[e\d\.\-]+$/) {return 0;}
	if ($in =~ /\-/) {my $t = $in; $t =~ s/[^\-]//g; if (length $t > 1) {return 0;}}
	if ($in =~ /\./) {my $t = $in; $t =~ s/[^\.]//g; if (length $t > 1) {return 0;}}
	if ($in =~ /e/) {my $t = $in; $t =~ s/[^e]//g; if (length $t > 1) {return 0;}}
	return 1;
}

sub is_valid_gradient_res {
	my $in = shift;
	unless ($in =~ /^\d+$/) {return 0;}
	if ($in < 1) {return 0;}
	return 1;
}

sub is_correct_state {
	my $in = shift;
	my $type = shift;
	if ($type =~ /^double_discrete$/i) {
		unless (length $in == 2) {return 0;}
		unless ($in =~ /^[01][01]$/) {return 0;}
	} elsif ($type =~ /^discrete$/i) {
		unless (length $in == 1) {return 0;}
		unless ($in =~ /^[01]$/) {return 0;}
	}
	return 1;
}

sub is_probability {
	my $in = shift;
	unless (is_numeric($in)) {return 0;}
	unless (($in <= 1) && ($in >= 0)) {return 0;}
	return 1;
}

sub get_support_average_of_descendents {
	my $node = shift;
	my $support_type = shift;
	my $a;
	foreach my $subnode ($node->each_Descendent) {
		foreach my $r ($subnode->get_tag_values('support')) {
			if (defined($r->{$support_type})) {push @$a,$r->{$support_type};}
		}
	}
	my $out;
	foreach my $b (@$a) {$out += $b/scalar(@$a);}
	return $out;
}

sub is_binary {
	my $in = shift;
	if ($in =~ /^[01]$/) {return 1;} else {return 0;}
}

###################################################################################################################################
#   xml parsing
###################################################################################################################################

sub phyloxml2trees {
	my $file = shift;
	my $trees;
	my $phyloxml = XML::DOM::Parser->new->parsefile($file)->getChildNodes->[0];
	unless ($phyloxml->getNodeName =~ /^phyloxml$/i) {fatalerror("not a valid phyloXML file\n  file:   $file\n");}
	my $phylogenycounter; $phylogenycounter = 0;
	foreach my $phylogeny (grep {$_->getNodeName =~ /^phylogeny$/i} $phyloxml->getChildNodes) {
		my ($tree,$name,$root);
		++$phylogenycounter;
		foreach my $sub ($phylogeny->getChildNodes) {
			if ($sub->getNodeName =~ /^name$/i) {
				$name = $sub->getChildNodes->[0]->getNodeValue;
			} elsif ($sub->getNodeName =~ /^clade$/i) {
				$root = dom2node($sub);
				$tree = Bio::Tree::Tree->new(-root => $root);
				push @$trees,$tree;
			}
		}
	}
	unless ($phylogenycounter) {fatalerror("no phylogenies found in phyloXML\n  file:   $file\n");}
	return $trees;
}

sub dom2node {
	my $in = shift;
	my $node = Bio::Tree::Node->new();
	foreach my $sub (grep {!($_->getNodeName eq '#text')} $in->getChildNodes) {
		if ($sub->getNodeName =~ /^name$/i) {
			$node->id($sub->getChildNodes->[0]->getNodeValue);
		} elsif ($sub->getNodeName =~ /^clade$/i) {
			$node->add_Descendent(dom2node($sub));
		} elsif ($sub->getNodeName =~ /^branch$/i) {
			my $sns = $sub->getChildNodes;
			foreach my $sn (grep {!($_->getNodeName eq '#text')} @$sns) {
				if ($sn->getNodeName =~ /^length$/i) {
					$node->branch_length($sn->getChildNodes->[0]->getNodeValue);
				} elsif ($sn->getNodeName =~ /^support$/i) {
					my $type = $sn->getAttributes->{'type'}->getChildNodes->[0]->getNodeValue;
					$node->add_tag_value('support',{$type => $sn->getChildNodes->[0]->getNodeValue});
				}
			}
		} elsif ($sub->getNodeName =~ /^custom$/i) {
			my $sns = $sub->getChildNodes;
			my ($tag);
			foreach my $sn (grep {!($_->getNodeName eq '#text')} @$sns) {
				if ($sn->getNodeName =~ /^name$/i) {
					$tag->{'name'} = $sn->getChildNodes->[0]->getNodeValue;
					unless (defined($tag->{'name'})) {fatalerror("Invalid input: One of the custom tags has no name\n");}
				} elsif ($sn->getNodeName =~ /^value$/i) {
					$tag->{'value'} = $sn->getChildNodes->[0]->getNodeValue;
					unless (defined($tag->{'value'})) {fatalerror("Invalid input: The value of tag ".$tag->{'name'}." is not defined for one of the nodes\n");}
					unless (is_numeric($tag->{'value'})) {fatalerror("Invalid input: The value of tag ".$tag->{'name'}." is not numeric for one of the nodes: ".$tag->{'value'}."\n");}
				} elsif ($sn->getNodeName =~ /^type$/i) {
					$tag->{'type'} = $sn->getChildNodes->[0]->getNodeValue;
					unless (defined($tag->{'type'})) {fatalerror("Invalid input: The character type of tag ".$tag->{'name'}." is not defined for one of the nodes\n");}
					unless ($recognized->{'vartype'}->{'all'}->{$tag->{'type'}}) {fatalerror("Invalid input: The character type of tag ".$tag->{'name'}." is not recognized for one of the nodes: ".$tag->{'type'}."\n");}
				} elsif ($sn->getNodeName =~ /^state$/i) {
					my $q;
					foreach my $in (grep {!($_->getNodeName eq '#text')} @{$sn->getChildNodes}) {
						if ($in->getNodeName =~ /^state$/) {
							$q->{'state'} = $in->getChildNodes->[0]->getNodeValue;
						} elsif ($in->getNodeName =~ /^probability$/) {
							$q->{'probability'} = $in->getChildNodes->[0]->getNodeValue;
						}
					}
					push @{$tag->{'state'}}, $q;
				}
			}
			if ($tag->{'state'}) {
				foreach my $q (@{$tag->{'state'}}) {
					unless (defined($q->{'state'})) {fatalerror("Invalid input: A state of tag ".$tag->{'name'}." is not defined for one of the nodes\n");}
					unless (is_correct_state($q->{'state'},$tag->{'type'})) {fatalerror("Invalid input: A state of tag ".$tag->{'name'}." is not recognized for one of the nodes: ".$q->{'state'}."\n");}
					unless (defined($q->{'probability'})) {fatalerror("Invalid input: A probability of tag ".$tag->{'name'}." is not defined for one of the nodes\n");}
					unless (is_probability($q->{'probability'})) {fatalerror("Invalid input: A probability of tag ".$tag->{'name'}." is of a wrong format for one of the nodes: ".$q->{'probability'}."\n");}
				}
			}
			if ($tag) {
				$node->add_tag_value($tag->{'name'},$tag);
			}
			unless ($variable_type) {
				if ((defined $colorinfo) && ($tag->{'name'} eq $colorinfo)) {
					$variable_type = $tag->{'type'};
				}
			}
		} else {
			message("  not parsing: ".$sub->getNodeName."\n");
		}
	}
	push @$taxorder,$node;
	return $node;
}

###################################################################################################################################
#   messages subroutines
###################################################################################################################################

sub message {   # prints a message to STDOUT
	my $msg = shift;
	print STDOUT $msg;
}

sub fatalerror {   # prints an error message and terminates the script
	error(shift);
	exit;
}

sub error {   # prints an error message to STDERR
	my $msg = shift;
	print STDERR "\n---begin error message---\n";
	print STDERR $msg;
	print STDERR "---end error message---\n";
}
