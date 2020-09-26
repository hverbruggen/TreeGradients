#!/usr/bin/perl

#
#	TreeExtender.pl     version 1.02
#
#   heroen.verbruggen@gmail.com
#
# 
#   Version history:
#   1.03  parser for continuous characters from a simple list added
#         order of taxa in newick file is now retained
#   1.02  Bayestraits multistate parsing added
#   1.01  initial release

use warnings;
use strict;
use IO::File;
use Bio::TreeIO;
use XML::DOM;
use XML::Writer;
use Sort::Array qw(Sort_Table Discard_Duplicates);
use Data::Dumper;

message("\nTreeExtender.pl\n\n");

###################################################################################################################################
#   global variables
###################################################################################################################################

my (
	# filenames
	$treein,		# input tree file
	$treeformat,	# input tree file format
	$treeout,		# output tree filename
	$f1,			# file 1 for ancestral state parsing
	$f2,			# file 2 for ancestral state parsing

	# command-line parameters
	$nodeid,		# defines what the node IDs represent
	$parser,		# type of parser (ancml, ape, bt, mesq, ...)
	$varname,		# name of variable
	$vartype,		# type of variable
	$joinprob,		# specifies whether and how double_discrete needs to be converted to discrete
	$multistates,	# sorted array containing the states that the multistate character being parsed can take
	
	# data container
	$trees,			# array of trees in Bio::Tree::Tree format
	$taxorder,		# whether taxon sorting is in effect
	
);

my $defaults = {
	'treeformat' => 'newick',
	'nodeid' => 'bootstrap',
};

my $recognized = {
	'suffixes' => {	'nwk' => 'newick',		# all of these should be in lower case
					'newick' => 'newick',
					'nhx' => 'nhx',
					'nexus' => 'nexus',
					'nex' => 'nexus',
					'tre' => 'nexus',
					},
	'treeformats' => {	'newick' => 1,		# all of these should be in lower case
						'phyloxml' => 1,
						},
	'node_id' => {	'posterior' => 1,		# all of these should be in lower case
					'bootstrap' => 1,
					'id' => 1,
					},
	'parsers' => {	'ancml' => 1,
					'ape' => 1,
					'bt' => 1,
					'mesq' => 1,
					'list' => 1,
					},
	'vartype' => {
		'mesq' => {	'discrete' => 1,
					'continuous' => 1,
					},
		'bt' => {	'double_discrete' => 1,
					'multistate' => 1,
					},
		'list' => {	'continuous' => 1,
					},
		'all' => {	'discrete' => 1,
					'continuous' => 1,
					'double_discrete' => 1,
					'multistate' => 1,
		},
	},
	'joinprob' => {
		1 => 1,
		2 => 1,
	},
};

###################################################################################################################################
#   command line parsing
###################################################################################################################################

	unless (	($ARGV[0]) and (substr($ARGV[0],0,1) eq "-") and
				defined($ARGV[1]) )
		{ usage(); }
	
	for (my $i=0; $i<scalar(@ARGV); $i+=2) {
		if ($ARGV[$i] =~ /^\-i$/i) {$treein=$ARGV[$i+1];}
		elsif ($ARGV[$i] =~ /^\-o$/i) {$treeout=$ARGV[$i+1];}
		elsif ($ARGV[$i] =~ /^\-tf$/i) {$treeformat=$ARGV[$i+1];}
		elsif ($ARGV[$i] =~ /^\-id$/i) {$nodeid=$ARGV[$i+1];}
		elsif ($ARGV[$i] =~ /^\-p$/i) {$parser=$ARGV[$i+1];}
		elsif ($ARGV[$i] =~ /^\-f1$/i) {$f1=$ARGV[$i+1];}
		elsif ($ARGV[$i] =~ /^\-f2$/i) {$f2=$ARGV[$i+1];}
		elsif ($ARGV[$i] =~ /^\-vn$/i) {$varname=$ARGV[$i+1];}
		elsif ($ARGV[$i] =~ /^\-vt$/i) {$vartype=$ARGV[$i+1];}
		elsif ($ARGV[$i] =~ /^\-jp$/i) {$joinprob=$ARGV[$i+1];}
		else { usage(); }
	}
	unless ($treein) {usage()}
	unless (-e $treein) {fatalerror("file not found\n  file:    $treein\n")}
	unless ($treeformat) {$treeformat = $defaults->{'treeformat'};}
	$treeformat = lc $treeformat;
	unless ($recognized->{'treeformats'}->{$treeformat}) {fatalerror("specified tree format not supported:\n  format:  $treeformat\n")}
	unless ($treeout) {$treeout = default_treeout_filename($treein)}
	unless ($nodeid) {$nodeid = $defaults->{'nodeid'};}
	$nodeid = lc $nodeid;
	unless ($recognized->{'node_id'}->{$nodeid}) {fatalerror("the specified value for the -id option is not recognized\n  specified value:  $nodeid\n  allowed values:   ".join(', ',keys(%{$recognized->{'node_id'}}))."\n");}
	if (defined($joinprob)) {unless($recognized->{'joinprob'}->{$joinprob}) {fatalerror("the specified value for the -jp option is not recognized\n  specified value:  $joinprob\n  allowed values:   ".join(', ',sort(keys(%{$recognized->{'joinprob'}})))."\n");}}
	if ($parser) {
		$parser = lc $parser;
		unless ($recognized->{'parsers'}->{$parser}) {fatalerror("the parser specified with the -p option is invalid\n  specified parser:  $parser\n  allowed values:   ".join(', ',keys(%{$recognized->{'parsers'}}))."\n");}
		if ($parser eq 'mesq') {
			unless ($f1) {fatalerror("when using the Mesquite parser, the file containing the Mesquite output must be specified with the -f1 option\n");}
			unless (-e $f1) {fatalerror("file not found\n  file:    $treein\n");}
			unless ($vartype) {fatalerror("when parsing Mesquite output, the variable type (discrete or continuous) must be specified with the -vt option\n");}
			unless ($recognized->{'vartype'}->{$parser}->{$vartype}) {fatalerror("the variable type specified with the -vt option was not recognized\n  specified value:   $vartype\n   allowed values for $parser parser:   ".join(', ',keys(%{$recognized->{'vartype'}->{$parser}}))."\n");}
		}
		if ($parser eq 'bt') {
			unless ($f1) {fatalerror("when using the BayesTraits parser, the file containing the BayesTraits output must be specified with the -f1 option\n");}
			unless (-e $f1) {fatalerror("file not found\n  file:    $f1\n");}
			unless ($f2) {fatalerror("when using the BayesTraits parser, the file containing character states for the terminal taxa must be specified with the -f2 option\n");}
			unless (-e $f2) {fatalerror("file not found\n  file:    $f2\n");}
			unless ($vartype) {fatalerror("when parsing mesquite output, the variable type must be specified with the -vt option\n  currently implemented variable types: ".join(', ',keys(%{$recognized->{'vartype'}->{$parser}}))."\n");}
			unless ($recognized->{'vartype'}->{$parser}->{$vartype}) {fatalerror("the variable type specified with the -vt option was not recognized\n  specified value:   $vartype\n   allowed values for $parser parser:   ".join(', ',keys(%{$recognized->{'vartype'}->{$parser}}))."\n");}
		}
		if ($parser eq 'list') {
			unless ($f1) {fatalerror("when using the simple list parser, the file containing the values must be specified with the -f1 option\n");}
			unless (-e $f1) {fatalerror("file not found\n  file:    $f1\n");}
			unless ($vartype) {fatalerror("when parsing a simple list, the variable type must be specified with the -vt option\n  currently implemented variable types: ".join(', ',keys(%{$recognized->{'vartype'}->{$parser}}))."\n");}
			unless ($recognized->{'vartype'}->{$parser}->{$vartype}) {fatalerror("the variable type specified with the -vt option was not recognized\n  specified value:   $vartype\n   allowed values for $parser parser:   ".join(', ',keys(%{$recognized->{'vartype'}->{$parser}}))."\n");}
		}
	}
	
sub usage {			# prints warning about incorrect command line usage
	print "\nusage:\n"; 
	print "\nmandatory parameters\n";
	print "   -i    input tree file\n";
	print "\noptional parameters\n";
	print "   -o    output tree file\n";
	print "            output tree is always saved in phyloxml format\n";
	print "            if no output filename is specified by user, it will be based on the input filename\n";
	print "   -tf   format of input tree file\n";
	print "            newick\n";
	print "            phyloxml\n";
	print "            [default: ",$defaults->{'treeformat'},"]\n";
	print "   -id   what the IDs of internal nodes in the newick/nexus trees represent\n";
	print "            ID         clade names\n";
	print "            bootstrap  bootstrap values\n";
	print "            posterior  posterior probabilities\n";
	print "            [default: ",$defaults->{'nodeid'},"]\n";
	print "   -p    file parser to use\n";
#	print "            ancml      Ancestor states for continuous traits using Maximum Likelihood\n";
#	print "            ape        Analysis of Phylogenetics and Evolution\n";
	print "            bt         BayesTraits\n";
	print "            mesq       Mesquite\n";
	print "            list       list file\n";
	print "\nMesquite parser-specific information\n";
	print "   -f1   file containing Mesquite ancestor states output\n";
	print "   -vt   variable type that is being parsed\n";
	print "            continuous\n";
	print "            discrete\n";
	print "   -vn   name of the variable that is being parsed (without spaces or funny characters)\n";
	print "            if unspecified, the program uses the variable name specified in the Mesquite output\n";
	print "\nBayesTraits parser-specific information\n";
	print "   -f1   file containing BayesTraits output\n";
	print "   -f2   file containing character states of terminal taxa (file used as input for BayesTraits)\n";
	print "   -vt   variable type that is being parsed\n";
	print "            double_discrete    simultaneous analysis of two discrete characters (dependent or independent)\n";
	print "            multistate         a single multistate character (maximum of three states)\n";
	print "   -jp   option that allows plotting the probability of a single discrete trait for double_discrete input\n";
	print "            1    plots the probabilities for the first character\n";
	print "            2    plots the probabilities for the second character\n";
	print "   -vn   name of the variable that is being parsed (without spaces or funny characters)\n";
	print "            if unspecified, the program uses 'BayesTraits' as the variable name\n";
	print "\nList file parser-specific information\n";
	print "   -f1   file containing list of states for each node\n";
	print "   -vt   variable type that is being parsed\n";
	print "            continuous\n";
	print "   -vn   name of the variable that is being parsed (without spaces or funny characters)\n";
	print "            if unspecified, the program uses 'list' as the variable name\n";
	print "\n";
	exit;
}


###################################################################################################################################
#   main program
###################################################################################################################################

$trees = open_trees($treein,$treeformat);
testtrees($trees);
unless ($treeformat eq 'phyloxml') {
	$trees = process_nodeID_support($trees);
}
if ($parser) {
	if (($parser =~ /^mesq$/i) and ($vartype =~ /^discrete$/i)) {
		fatalerror("Mesquite parser not working yet\n");
		$trees = mesquite_discrete_alltrees($f1,$trees);
	} elsif ($parser =~ /^bt$/i) {
		unless ($varname) {$varname = 'BayesTraits';}
		message("parsing BayesTraits $vartype\n");
		$trees = parse_bayestraits_alltrees($vartype,$f1,$f2,$trees);
		message("  done\n");
		if (($vartype =~ /^double_discrete$/i) && defined($joinprob)) {
			message("jp option activated: converting data to probabilities for the desired character alone\n");
			message("  calculating state probabilities for character $joinprob across states of character ".(3 - $joinprob)."\n");
			$trees = bayestraits_convert_double_discrete_to_discrete_alltrees($trees,$varname,$joinprob);
			message("  done\n");
		}
	} elsif ($parser =~ /^list$/i) {
		unless ($varname) {$varname = 'list';}
		message("parsing list $vartype\n");
		$trees = parse_list_alltrees($vartype,$f1,$trees);
		message("  done\n");
	}
}
trees2phyloxml($trees,$treeout);


###################################################################################################################################
#   IO subroutines
###################################################################################################################################

# this subroutine opens the input tree file
sub open_trees {
	my $file = shift;
	my $format = shift;
	my $trees;
	message("opening tree file\n");
	if ($format =~ /^phyloxml$/i) {
		$trees = phyloxml2trees($file);
	} else {
		$trees = open_with_bioperl($file,$format);
	}
	if ($trees) {
		print "  ",scalar(@$trees)," trees loaded\n";
	} else {
		fatalerror("TreeExtender.pl was unable to extract trees from input file\n  file: $file\n");
	}
	determine_taxorder($file,$format,$trees);
	return $trees;
}

sub determine_taxorder {
	message("determining taxon order\n");
	my $file = shift;
	my $format = shift;
	my $trees = shift;
	unless ($format =~ /^newick$/i) {
		message("  taxon order may not be preserved -- use newick file as input if you want to retain taxon order\n");
		return undef;
	}
	unless (open FH,$file) {fatalerror("cannot read from file $file\n")}
	my @a = <FH>; close FH;
	my $l = join "",@a;
	$l =~ s/[\r\n]//g;
	my @tr = split ';',$l;
	unless (scalar @tr == scalar @$trees) {fatalerror("number of trees in file does not match number of trees in memory\n");}
	for (my $i = 0; $i < scalar @tr; ++$i) {
		my $out;
		my $str = $tr[$i];
		$str =~ s/\)[\d\.\-e]+/\)/gi;
		$str =~ s/\:[\d\.\-e]+//gi;
		$str =~ s/[\(\)\,]/ /g;
		$str =~ s/\s+/ /g;
		$str =~ s/^\s+//g;
		$str =~ s/\s+$//g;
		my @p = split ' ',$str;
		my $counter; $counter = 0;
		foreach my $n (@p) {
			++$counter;
			$out->{$n} = $counter;
		}
		foreach my $n ($trees->[$i]->get_leaf_nodes) {
			unless ($out->{$n->id}) {
				message("  error encountered while determining taxon order -- taxon order may not be preserved\n");
				return undef;
			}
			$n->{'_taxorder'} = $out->{$n->id};
		}
	}
	message("  done for ".scalar(@tr)." trees\n");
	$taxorder = 1;
}

sub open_with_bioperl {
	my $file = shift;
	my $format = shift;
	my $trees;
	my $treeio = Bio::TreeIO->new(-file => $file, -format => $format);
	while (my $tree = $treeio->next_tree) {
		push @$trees, $tree;
	}
	unless ($trees) {fatalerror("BioPerl was unable to extract trees\n  file:    $file\n  format:  $format\n");}
	return $trees;
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
					#unless (is_correct_state($q->{'state'},$tag->{'type'})) {fatalerror("Invalid input: A state of tag ".$tag->{'name'}." is not recognized for one of the nodes: ".$q->{'state'}."\n");}
					unless (defined($q->{'probability'})) {fatalerror("Invalid input: A probability of tag ".$tag->{'name'}." is not defined for one of the nodes\n");}
					unless (is_probability($q->{'probability'})) {fatalerror("Invalid input: A probability of tag ".$tag->{'name'}." is of a wrong format for one of the nodes: ".$q->{'probability'}."\n");}
				}
			}
			if ($tag) {
				$node->add_tag_value($tag->{'name'},$tag);
			}
		} else {
			message("  not parsing: ".$sub->getNodeName."\n");
		}
	}
	return $node;
}

###################################################################################################################################
#   simple list file parsing
###################################################################################################################################

sub parse_list_alltrees {
	my $vartype = shift;
	my $file = shift;
	my $trees = shift;
	unless ($vartype eq 'continuous') {fatalerror("List file parsing is only implemented for continuous characters\n");}
	my $data = parse_list_file($file);
	foreach my $tree (@$trees) {
		$tree = annotate_tree_from_list($tree,$data)
	}
	return $trees;
}

sub annotate_tree_from_list {
	my $tree = shift;
	my $data = shift;
	foreach my $node ($tree->get_nodes) {
		my $taxalist;
		if ($node->is_Leaf) {
			$taxalist = $node->id;
		} else {
			my @nm;
			foreach my $t (grep {$_->is_Leaf} $node->get_all_Descendents) {
				push @nm,$t->id;
			}
			$taxalist = join ' ',sort @nm;
		}
		my $datum = get_list_value($data,$taxalist);
		if (defined $datum) {
			my $tag = {
				type => 'continuous',
				name => $varname,
				value => $datum,
			};
			$node->add_tag_value($tag->{name},$tag);
		} else {
			fatalerror("no data value not found for common ancestor of these taxa:\n$taxalist\n");
		}
	}
	return $tree;
}

sub get_list_value {
	my $data = shift;
	my $taxalist = shift;
	foreach my $datum (@$data) {
		if ($datum->{taxalist} eq $taxalist) {return $datum->{value}}
	}
	return undef;
}

sub parse_list_file {
	my $file = shift;
	unless (open FH,$file) {fatalerror("Cannot read from list file $file\n");}
	my @a = <FH>; close FH;
	my $out;
	foreach my $line (@a) {
		$line =~ s/[\r\n]//g;
		unless ($line =~ /^\s*$/) {
			if ($line =~ /^([\d\.\e\-]+)\t(.*)$/i) {
				my $val = $1; my $taxalist = $2;
				if ($taxalist =~ /\s/) {
					$taxalist = join ' ', sort split /\s+/,$taxalist;
				}
				push @$out,{taxalist => $taxalist,value => $val};
			} else {
				fatalerror("List file ill-formatted in this line\n$line\n");
			}
		}
	}
	return $out;
}

###################################################################################################################################
#   bayestraits parsing
###################################################################################################################################

sub parse_bayestraits_alltrees {
	my $vartype = shift;
	my $bt_output_file = shift;
	my $bt_input_data_file = shift;
	my $trees = shift;
	if (scalar @$trees > 1) {fatalerror("BayesTraits parsing is currently implemented for a single tree only\n");}
	foreach my $tree (@$trees) {
		if ($vartype =~ /^double_discrete$/i) {
			$tree = bayestraits_double_discrete($bt_output_file,$bt_input_data_file,$tree);
		} elsif ($vartype =~ /^multistate$/i) {
			$tree = bayestraits_multistate($bt_output_file,$bt_input_data_file,$tree);
		}
	}
	return $trees;
}

sub bayestraits_multistate {
	my $bt_output_file = shift;
	my $bt_input_data_file = shift;
	my $tree = shift;
	my $internal_states = parse_bt_output_multistate($bt_output_file);
	$tree = check_and_annotate_bt_multistate($internal_states,$tree);
	my $terminal_states = parse_bt_input_data_multistate($bt_input_data_file);
	$tree = check_and_annotate_bt_terminal_multistate($terminal_states,$tree);
	return $tree;
}

sub remove_duplicates_from_state_reconstructions {
	my $states = shift;
	my $strings;
	my $out;
	$out->{'Root'}=$out->{'Root'};
	foreach my $key (keys %$states) {
		unless ($key eq 'Root') {
			my $str = join ', ',sort(@{$states->{$key}->{'desc'}});
			if (defined $strings) {
				unless ($strings->{$str}) {
					$strings->{$str} = 1;
					$out->{$key} = $states->{$key};
				}
			} else {
				$strings->{$str} = 1;
				$out->{$key} = $states->{$key};
			}
		}
	}
	return $out
}

sub check_and_annotate_bt_multistate {
	my $states = shift;
	my $tree = shift;
	unless (scalar(keys %$states) - 1 == scalar(grep {!($_->is_Leaf)} $tree->get_nodes)) {
		$states = remove_duplicates_from_state_reconstructions($states);
		unless (scalar(keys %$states) - 1 == scalar(grep {!($_->is_Leaf)} $tree->get_nodes)) {
			fatalerror("number of nodes in tree does not correspond to number of nodes for which ancestral states have been reconstructed\n  in tree: ".scalar(grep {!($_->is_Leaf)} $tree->get_nodes)." -- reconstructed: ".(scalar(keys %$states)-1)." (".scalar(keys %$states)." - 1 for root node)\n");
		}
	}
	my $root_tag;
	$root_tag->{'name'} = $varname || 'BayesTraits';
	$root_tag->{'type'} = 'multistate';
	foreach my $state (keys %{$states->{'Root'}}) {
		unless ($state eq "desc") {
			$state =~ /P(\d+)/;
			push @{$root_tag->{'state'}}, {'state' => $1, 'probability' => $states->{'Root'}->{$state}};
		}
	}
	$tree->get_root_node->add_tag_value($root_tag->{'name'},$root_tag);
	foreach my $node (grep {!($_->is_Leaf)} $tree->get_nodes) {
		if ($node->ancestor) {
			my $tag;
			$tag->{'name'} = $varname || 'BayesTraits';
			$tag->{'type'} = 'multistate';
			my $continue; $continue = 1;
			foreach my $tagname ($node->get_all_tags) {
				if ($tagname eq $tag->{'name'}) {
					$continue = 0;
				}
			}
			if ($continue) {
				my $desc;
				foreach my $terminal (grep {$_->is_Leaf} $node->get_all_Descendents) {
					push @$desc,$terminal->id;
				}
				$desc = join ' ',sort(@$desc);
				my $found; $found = 0;
				foreach my $key (keys %$states) {
					unless ($key eq 'Root') {
						my $list = join ' ',sort(@{$states->{$key}->{'desc'}});
						if ($list eq $desc) {
							$found = 1;
							foreach my $state (keys %{$states->{$key}}) {
								unless ($state eq "desc") {
									$state =~ /P(\d+)/;
									push @{$tag->{'state'}}, {'state' => $1, 'probability' => $states->{$key}->{$state}};
								}
							}
							$node->add_tag_value($tag->{'name'},$tag);
						}
					}
				}
				unless ($found) {fatalerror("the node with these descendents was not found in the BayesTraits output\n  $desc\n");}
			}
		}
	}
	return $tree;
}

sub parse_bt_output_multistate {
	my $file = shift;
	print "  getting inferred ancestral character state probabilities from $file\n";
	unless (open FH,$file) {fatalerror("cannot read from file: $file\n");}
	my @a = <FH>; my $all = join '',@a; $all =~ s/[\r\n]/\n/g;
	close FH;
	my $out;
	unless ($all =~ /Model:\s+Multistate\s*/i) {fatalerror("you did not use multistate for your BayesTraits analysis\n");}
	unless ($all =~ /Analy*sis Type:\s+Maximum Likelihood\s/) {fatalerror("only maximum likelihood results can be parsed at present\n");}
	$all =~ s/Node:/;Node:/gs;
	$all =~ s/([NM]RCA):/;$1:/gs;
	$all =~ s/Fossil:/;Fossil:/gs;
	$all =~ s/Removed\s+taxa/;Removed taxa/igs;
	$all =~ s/Tree\s+information/;Tree Information/igs;
	my $array; @$array = split ';',$all;
	do {shift @$array} until (($array->[0] =~ /^\s*Node:/) or ($array->[0] =~ /^\s*Fossil:/) or ($array->[0] =~ /^\s*[NM]RCA:/));
	do {pop @$array} until (($array->[scalar(@$array)-1] =~ /^\s*Node:/) or ($array->[scalar(@$array)-1] =~ /^\s*Fossil:/) or ($array->[scalar(@$array)-1] =~ /^\s*[NM]RCA:/));
	foreach my $entry (@$array) {
		my $arr; @$arr = split /\n/,$entry;
		my $first = shift @$arr;
		$first =~ /:\s+([^\s]+)\s+/;
		my $name = $1;
		my $desc;
		foreach my $line (@$arr) {
			unless ($line =~ /^\s*$/) {
				$line =~ /^\s+\d+\s+([^\s]+)\s*$/;
				push @$desc,$1;
			}
		}
		$out->{$name}->{'desc'} = $desc;
	}
	$all =~ /\n(Tree No\t.+?)\n(.+)$/s;
	my $headers; @$headers = split /\t/,$1;
	unless (defined $multistates) {
		my $hash;
		foreach my $header (@$headers) {
			if ($header =~ /(.+)\s+P\((\d+)\)/i) {
				$hash->{$2} = 1;
			}
		}
		@$multistates = sort keys %$hash;
	}
	my $data_text = $2;
	my $data_lines; @$data_lines = split ' ___',$data_text;
	foreach my $line (@$data_lines) {
		unless ($line =~ /^\s*$/) {
			$line =~ s/[\n\r]//g;
			my $data; @$data = split /\t/,$line;
			unless (scalar @$headers == scalar @$data) {fatalerror("file ill-formatted: $file\n  the number of values in the following line didn't match the number of labeled data columns\n$line\n");}
			for (my $i = 0; $i < scalar(@$headers); ++$i) {
				my $header = $headers->[$i];
				my $datum = $data->[$i];
				if ($datum =~ /([\.\d]+)e(\-*\d+)/i) {$datum = convert_from_scientific_format($datum);}
				if ($header =~ /(.+)\s+P\((\d+)\)/i) {
					$out->{$1}->{'P'.$2} = $datum;
				}
			}
		}
	}
	return $out;
}


sub check_and_annotate_bt_terminal_multistate {
	my $states = shift;
	my $tree = shift;
	unless (scalar(keys %$states) == scalar($tree->get_leaf_nodes)) {fatalerror("number of taxa in tree and in file with character states of terminal taxa does not match\n");}
	my $present_states; foreach my $taxon (keys %$states) {$present_states->{$states->{$taxon}} = 1;}
	unless ((scalar keys %$present_states) < 4) {fatalerror("the present version of this script allows a maximum of three character states for BayesTraits multistate parsing\n");}
	foreach my $taxon ($tree->get_leaf_nodes) {
		if (defined($states->{$taxon->id})) {
			my $state = $states->{$taxon->id};
			my $tag;
			$tag->{'name'} = $varname || 'BayesTraits';
			$tag->{'type'} = 'multistate';
			foreach my $state (@$multistates) {
				my $p; if ($state eq $states->{$taxon->id}) {$p = 1} else {$p = 0;}
				push @{$tag->{'state'}}, {'state' => $state, 'probability' => $p};
			}
			$taxon->add_tag_value($tag->{'name'},$tag);
		} else {
			fatalerror("taxon ".$taxon->id.", which is present in the tree, was not found in the file with character states of terminal taxa\n");
		}
	}
	return $tree;
}

sub parse_bt_input_data_multistate {
	my $file = shift;
	print "  getting character states of terminal taxa from $file\n";
	unless (open FH,$file) {fatalerror("cannot read from file: $file\n");}
	my @a = <FH>;
	close FH;
	my $out;
	foreach my $line (@a) {
		$line =~ s/[\r\n]//g;
		unless ($line =~ /^\s*$/) {
			if ($line =~ /^\s*([^\s]+)\s+(\d+)\s*$/) {
				$out->{$1} = $2;
			} else {
				fatalerror("file ill-formatted: $file\n  cannot parse the following line\n$line\n");
			}
		}
	}
	return $out;
}

sub bayestraits_double_discrete {
	my $bt_output_file = shift;
	my $bt_input_data_file = shift;
	my $tree = shift;
	my $terminal_states = parse_bt_input_data_double_discrete($bt_input_data_file);
	$tree = check_and_annotate_bt_terminal_states($terminal_states,$tree);
	my $internal_states = parse_bt_output_double_discrete($bt_output_file);
	$tree = check_and_annotate_bt_internal_states($internal_states,$tree);
	return $tree;
}

sub check_and_annotate_bt_internal_states {
	my $states = shift;
	my $tree = shift;
	unless (scalar(keys %$states) - 1 == scalar(grep {!($_->is_Leaf)} $tree->get_nodes)) {
		fatalerror("number of nodes in tree does not correspond to number of nodes for which ancestral states have been reconstructed\n  in tree: ".scalar(grep {!($_->is_Leaf)} $tree->get_nodes)." -- reconstructed: ".(scalar(keys %$states)-1)." (".scalar(keys %$states)." - 1 for root node)\n");
	}
	my $root_tag;
	$root_tag->{'name'} = $varname || 'BayesTraits';
	$root_tag->{'type'} = 'double_discrete';
	push @{$root_tag->{'state'}}, {'state' => '00', 'probability' => $states->{'Root'}->{'P00'}};
	push @{$root_tag->{'state'}}, {'state' => '01', 'probability' => $states->{'Root'}->{'P01'}};
	push @{$root_tag->{'state'}}, {'state' => '10', 'probability' => $states->{'Root'}->{'P10'}};
	push @{$root_tag->{'state'}}, {'state' => '11', 'probability' => $states->{'Root'}->{'P11'}};
	$tree->get_root_node->add_tag_value($root_tag->{'name'},$root_tag);
	foreach my $node (grep {!($_->is_Leaf)} $tree->get_nodes) {
		if ($node->ancestor) {
			my $tag;
			$tag->{'name'} = $varname || 'BayesTraits';
			$tag->{'type'} = 'double_discrete';
			my $continue; $continue = 1;
			foreach my $tagname ($node->get_all_tags) {
				if ($tagname eq $tag->{'name'}) {
					$continue = 0;
				}
			}
			if ($continue) {
				my $desc;
				foreach my $terminal (grep {$_->is_Leaf} $node->get_all_Descendents) {
					push @$desc,$terminal->id;
				}
				$desc = join ' ',sort(@$desc);
				my $found; $found = 0;
				foreach my $key (keys %$states) {
					unless ($key eq 'Root') {
						my $list = join ' ',sort(@{$states->{$key}->{'desc'}});
						if ($list eq $desc) {
							$found = 1;
							push @{$tag->{'state'}}, {'state' => '00', 'probability' => $states->{$key}->{'P00'}};
							push @{$tag->{'state'}}, {'state' => '01', 'probability' => $states->{$key}->{'P01'}};
							push @{$tag->{'state'}}, {'state' => '10', 'probability' => $states->{$key}->{'P10'}};
							push @{$tag->{'state'}}, {'state' => '11', 'probability' => $states->{$key}->{'P11'}};
							$node->add_tag_value($tag->{'name'},$tag);
						}
					}
				}
				unless ($found) {fatalerror("the node with these descendents was not found in the BayesTraits output\n  $desc\n");}
			}
		}
	}
	return $tree;
}

sub parse_bt_output_double_discrete {
	my $file = shift;
	print "  getting inferred ancestral character state probabilities from $file\n";
	unless (open FH,$file) {fatalerror("cannot read from file: $file\n");}
	my @a = <FH>; my $all = join '',@a; $all =~ s/[\r\n]/ ___ /g;
	close FH;
	my $out;
	unless (($all =~ /Model:\s+Discrete\s/) || ($all =~ /Model:\s+Discete\s/)) {fatalerror("you did not use a double discrete (dependent or independent) for your BayesTraits analysis\n");}
	unless (($all =~ /Analysis Type:\s+Maximum Likelihood\s/) || ($all =~ /Analsis Type:\s+Maximum Likelihood\s/)) {fatalerror("only maximum likelihood results can be parsed at present\n");}
	$all =~ / ___ Node:\s+(.*)Tree Information/;
	my $array; @$array = split / ___ Node:\s+/,$1;
	foreach my $entry (@$array){
		$entry =~ /([^\s]+)\s+[\.\d]+\s* ___ (.*)/;
		my $name = $1;
		my $desc_text = $2;
		$desc_text =~ s/\s+___ //g;
		$desc_text =~ s/\s+\d+\s*/ /g;
		$desc_text =~ s/^\s*(.*)/$1/g;
		$desc_text =~ s/\s+/ /g;
		my $desc; @$desc = split /\s/,$desc_text;
		$out->{$name}->{'desc'} = $desc;
	}
	$all =~ / ___ (Tree No\t.+?) ___ (.+)$/;
	my $headers; @$headers = split /\t/,$1;
	my $data_text = $2;
	my $data_lines; @$data_lines = split ' ___',$data_text;
	foreach my $line (@$data_lines) {
		unless ($line =~ /^\s*$/) {
			my $data; @$data = split /\t/,$line;
			unless (scalar @$headers == scalar @$data) {fatalerror("file ill-formatted: $file\n  the number of values in the following line didn't match the number of labeled data columns\n$line\n");}
			for (my $i = 0; $i < scalar(@$headers); ++$i) {
				my $header = $headers->[$i];
				my $datum = $data->[$i];
				if ($datum =~ /([\.\d]+)e(\-*\d+)/i) {$datum = convert_from_scientific_format($datum);}
				if ($header =~ /(.+) - P\(([0-1])\,([0-1])\)/) {
					$out->{$1}->{'P'.$2.$3} = $datum;
				} elsif ($header =~ /(.+) - T([1-2]) - P\(([0-1])\)/) {
					my ($a,$b,$c) = ($1,$2,$3);
					$out->{$a}->{'T'.$b.'P'.$c} = $datum;
					
					if ((scalar(keys %{$out->{$a}}) > 4) or (($a =~ /\s*Root\s*/i) and (scalar(keys %{$out->{$a}}) == 4))) {
						$out->{$a}->{'P00'} = $out->{$a}->{'T1P0'} * $out->{$a}->{'T2P0'};
						$out->{$a}->{'P11'} = $out->{$a}->{'T1P1'} * $out->{$a}->{'T2P1'};
						$out->{$a}->{'P01'} = $out->{$a}->{'T1P0'} * $out->{$a}->{'T2P1'};
						$out->{$a}->{'P10'} = $out->{$a}->{'T1P1'} * $out->{$a}->{'T2P0'};
						delete $out->{$a}->{'T1P0'};
						delete $out->{$a}->{'T1P1'};
						delete $out->{$a}->{'T2P0'};
						delete $out->{$a}->{'T2P1'};
					}
				}
			}
		}
	}
	return $out;
}

sub check_and_annotate_bt_terminal_states {
	my $states = shift;
	my $tree = shift;
	unless (scalar(keys %$states) == scalar($tree->get_leaf_nodes)) {fatalerror("number of taxa in tree and in file with character states of terminal taxa does not match\n");}
	foreach my $taxon ($tree->get_leaf_nodes) {
		if ($states->{$taxon->id}) {
			my $tag;
			$tag->{'name'} = $varname || 'BayesTraits';
			$tag->{'type'} = 'double_discrete';
			push @{$tag->{'state'}}, {'state' => '00', 'probability' => $states->{$taxon->id}->{'P00'}};
			push @{$tag->{'state'}}, {'state' => '01', 'probability' => $states->{$taxon->id}->{'P01'}};
			push @{$tag->{'state'}}, {'state' => '10', 'probability' => $states->{$taxon->id}->{'P10'}};
			push @{$tag->{'state'}}, {'state' => '11', 'probability' => $states->{$taxon->id}->{'P11'}};
			$taxon->add_tag_value($tag->{'name'},$tag);
		} else {
			fatalerror("taxon ".$taxon->id.", which is present in the tree, was not found in the file with character states of terminal taxa\n");
		}
	}
	return $tree;
}

sub parse_bt_input_data_double_discrete {
	my $file = shift;
	print "  getting character states of terminal taxa from $file\n";
	unless (open FH,$file) {fatalerror("cannot read from file: $file\n");}
	my @a = <FH>;
	close FH;
	my $out;
	foreach my $line (@a) {
		$line =~ s/[\r\n]//g;
		unless ($line =~ /^\s*$/) {
			if ($line =~ /^\s*([^\s]+)\s+([0-1]+)\s+([0-1]+)\s*$/) {
				my $taxon = $1;
				if ($2 && $3) {
					$out->{$taxon}->{'P00'} = 0;
					$out->{$taxon}->{'P01'} = 0;
					$out->{$taxon}->{'P10'} = 0;
					$out->{$taxon}->{'P11'} = 1;
				} elsif ($2 && !$3) {
					$out->{$taxon}->{'P00'} = 0;
					$out->{$taxon}->{'P01'} = 0;
					$out->{$taxon}->{'P10'} = 1;
					$out->{$taxon}->{'P11'} = 0;
				} elsif (!$2 && $3) { 
					$out->{$taxon}->{'P00'} = 0;
					$out->{$taxon}->{'P01'} = 1;
					$out->{$taxon}->{'P10'} = 0;
					$out->{$taxon}->{'P11'} = 0;
				} elsif (!$2 && !$3) {
					$out->{$taxon}->{'P00'} = 1;
					$out->{$taxon}->{'P01'} = 0;
					$out->{$taxon}->{'P10'} = 0;
					$out->{$taxon}->{'P11'} = 0;
				}                   
			} else {
				fatalerror("file ill-formatted: $file\n  cannot parse the following line\n$line\n");
			}
		}
	}
	return $out;
}

sub bayestraits_convert_double_discrete_to_discrete_alltrees {
	my $trees = shift;
	my $varname = shift;
	my $char = shift;
	if (scalar @$trees > 1) {fatalerror("conversion of double_discrete to discrete is currently implemented for a single tree only\n");}
	foreach my $tree (@$trees) {
		$tree = bayestraits_convert_double_discrete_to_discrete($tree,$varname,$char);
	}
	$vartype = 'discrete';
	return $trees;
}

sub bayestraits_convert_double_discrete_to_discrete {
	my $tree = shift;
	my $varname = shift;
	my $char = shift;
	foreach my $node ($tree->get_nodes) {
		my $prob;
		my $orig; @$orig = $node->get_tag_values($varname);
		$orig = $orig->[0];
		unless ($node->remove_tag($varname)) {message("  unable to remove tag $varname from node\n");}
		my $tag;
		$tag->{'name'} = $orig->{'name'};
		$tag->{'type'} = 'discrete';
		my ($state0,$state1); ($state0,$state1) = (0,0);
		if ($char == 1) {
			foreach my $state (@{$orig->{'state'}}) {
				if ($state->{'state'} =~ /^0.$/) {
					$state0 += $state->{'probability'};
				} elsif ($state->{'state'} =~ /^1.$/) {
					$state1 += $state->{'probability'};
				}
			}
		} elsif ($char == 2) {
			foreach my $state (@{$orig->{'state'}}) {
				if ($state->{'state'} =~ /^.0$/) {
					$state0 += $state->{'probability'};
				} elsif ($state->{'state'} =~ /^.1$/) {
					$state1 += $state->{'probability'};
				}
			}
		}
		push @{$tag->{'state'}},{'state' => 0, 'probability' => $state0};
		push @{$tag->{'state'}},{'state' => 1, 'probability' => $state1};
		$node->add_tag_value($varname,$tag);
	}
	return $tree;
}

###################################################################################################################################
#   mesquite parsing
###################################################################################################################################

sub mesquite_discrete_alltrees {
	my $file = shift;
	my $trees = shift;
	foreach my $tree (@$trees) {
		$tree = mesquite_discrete($file,$tree);
	}
	return $trees;
}

sub mesquite_discrete {
	my $file = shift;
	my $tree = shift;
	my $nodes = mesquite_parse_tree_drawing($file);
	my $character;
	open(FH,$file) or fatalerror("unable to open file\n  file:   $file\n");
	my @filedump = <FH>; close FH;
	while (!($filedump[0] =~ /^Character:\s/)) {
		shift @filedump;
	}
	$filedump[0] =~ /^Character:\s(.*)$/;
	if ($varname) {
		$character = $varname;
	} else {
		$character = $1; $character =~ s/[\n\r]//g; $character =~ s/\s/_/g;
	}
	while (!($filedump[0] =~ /^node\s\d+\:/)) {
		shift @filedump;
	}
	while ($filedump[0] =~ /^node\s(\d+)\:(.*)$/) {
		my $nodenr = $1;
		$nodes->{'node'}->{$nodenr}->{'values'}->{$character}->{'type'} = 'discrete';
		my $inf; $inf = $2; $inf =~ s/\*//g;
		my @pairs = split(', ',$inf);
		foreach my $pair (@pairs) {
			$pair =~ /\s*(\d+)\:([\-\.\dE]+)/;
			my $state = $1;
			my $value = $2;
			if ($value =~ /([\.\d]+)E([\-\d]+)/i) {
				$value = $1 * (10**$2);
			}
			$nodes->{'node'}->{$nodenr}->{'values'}->{$character}->{'states'}->{'s_'.$state} = $value;
		}
		shift @filedump;
	}
	return mesquite_nodes2tree_discrete($nodes,$tree,$character);
}

sub mesquite_nodes2tree_discrete {
	my $nodes = shift;
	my $tree = shift;
	my $character = shift;
	foreach my $tr_node ($tree->get_nodes) {
		my $st_node = identify_node($tr_node,$nodes);
		if ($st_node->{'values'}->{$character}->{'type'} eq 'discrete') {
			my $tag;
			$tag->{'name'} = $character;
			$tag->{'type'} = 'discrete';
			foreach my $state (keys %{$st_node->{'values'}->{$character}->{'states'}}) {
				$tag->{'value'}->{$state} = $st_node->{'values'}->{$character}->{'states'}->{$state};
			}
			$tr_node->add_tag_value($character,$tag);
		}
	}
	return $tree;
}

sub identify_node {
	my $tr_node = shift;
	my $nodes = shift;
	my ($terminals,$internals);
	foreach my $st_node (@{$nodes->{'node'}}) {
		if ($st_node->{'descendents'}) {
			push @$internals,$st_node;
		} else {
			push @$terminals,$st_node;
		}
	}
	if ($tr_node->is_Leaf) {
		foreach my $st_node (@$terminals) {
			if ($st_node->{'taxon_label'} eq $tr_node->id) {
				return $st_node;
			}
		}
		fatalerror("the taxon ".$tr_node->id.", which is present in the tree, was not found in the Mesquite data\n");
	} else {
		my $desc_hash;
		foreach my $terminal (grep {$_->is_Leaf} $tr_node->get_all_Descendents) {
			$desc_hash->{$terminal->id} = 1;
		}
		foreach my $st_node (@$internals) {
			if ((scalar keys %$desc_hash) == (scalar @{$st_node->{'descendents'}->{'descendent'}})) {
				my $counter; $counter = 0;
				foreach my $desc (@{$st_node->{'descendents'}->{'descendent'}}) {
					if ($desc_hash->{$desc}) {
						++$counter;
					}
				}
				if ($counter == scalar keys %$desc_hash) {
					return $st_node;
				}
			}
		}
		fatalerror("the clade containing the following taxa was found in the tree but not in the Mesquite data\n".join("\n  ",keys(%$desc_hash))."\n");
	}
}

sub mesquite_parse_tree_drawing {
	my $file = shift;
	my ($nodes,$tree_drawing);
	# open file and skip to "Tree with node numbers"
		open(FH,$file) or fatalerror("unable to open file\n  file:   $file\n");
		my @filedump = <FH>; close FH;
		while (!($filedump[0] =~ /Tree with node numbers:/)) {
			shift @filedump;
		}
		shift @filedump;
	# read "Tree with node numbers" into $tree_drawing
		while (!($filedump[0] =~ /^\s*$/)) {
			push @$tree_drawing,$filedump[0];
			shift @filedump;
		}
	# parse labels of terminal taxa
		my ($tiplabels,$tiplabels_inv) = get_tip_labels($tree_drawing);
		foreach my $node_nr (keys %$tiplabels) {
			$nodes->{'node'}->{$node_nr}->{'type'} = 'terminal';
			$nodes->{'node'}->{$node_nr}->{'taxon_label'} = $tiplabels->{$node_nr};
		}
	# parse topology
		my $edges = get_edges($tree_drawing);
		foreach my $parent (keys %$edges) {
			$nodes->{'node'}->{$parent}->{'type'} = 'internal';
			my $list;
			$list->{$parent} = 1;
			$list = mesquite_descendents_list($list,$edges);
			foreach my $desc (keys %$list) {
				push @{$nodes->{'node'}->{$parent}->{'descendents'}->{'descendent'}},$tiplabels->{$desc};
			}
		}
	return $nodes;
}

sub mesquite_descendents_list {
	my $in = shift;
	my $edges = shift;
	my $out;
	foreach my $nodenr (keys %$in) {
		if ($edges->{$nodenr}) {
			delete $in->{$nodenr};
			foreach my $nr1 (@{$edges->{$nodenr}}) {
				my $exp; $exp->{$nr1} = 1;
				my $get = mesquite_descendents_list($exp,$edges);
				foreach my $nr2 (keys %$get) {
					$out->{$nr2} = 1;
				}
			}
		} else {
			$out->{$nodenr} = 1;
		}
	}
	return $out;
}


###################################################################################################################################
#   xml generating
###################################################################################################################################

# the following subroutine
#   opens the IO::File
#   takes care of the big structural XML elements
#   calls the node2xml subroutine on the root node of every tree
sub trees2phyloxml {
	message("generating XML file\n");
	my $trees = shift;
	my $file = shift;
	my $output = IO::File->new(">$file");
	my $xml = XML::Writer->new(DATA_MODE => 'true', DATA_INDENT => 2, OUTPUT => $output);
	$xml->xmlDecl('UTF-8');
	$xml->startTag('phyloxml');
	my $treecounter; $treecounter = 0;
	foreach my $tree (@$trees) {
		++$treecounter;
		$xml->startTag('phylogeny');
		$xml->startTag('name');
		if ($tree->id) {
			$xml->characters($tree->id);
		} else {
			$xml->characters('no_name_'.$treecounter);
		}
		$xml->endTag('name');
		node2xml($tree->get_root_node,$xml);
		$xml->endTag('phylogeny');
	}
	$xml->endTag('phyloxml');
	message("  done\n");
}

# this subroutine generates an xml structure for a single node and calls itself to nest the xml structure of any daughter nodes within the xml structure of the present node
sub node2xml {
	my $node = shift;
	my $xml = shift;
	# start present node
		$xml->startTag('clade');
	# annotate name of present node
		if ($node->id) {
			if ($node->is_Leaf || !($recognized->{'support_types'}->{$nodeid})) {
				$xml->startTag('name');
				$xml->characters($node->id);
				$xml->endTag('name');
			}
		}
	# annotate branch properties for present node
		if ($node->branch_length or $node->has_tag('support')) {
			$xml->startTag('branch');
			# annotate branch length
				if ($node->branch_length) {
					$xml->startTag('length');
					$xml->characters($node->branch_length);
					$xml->endTag('length');
				}
			# annotate node support
				if ($node->has_tag('support')) {
					my @values = $node->get_tag_values('support');
					my $support = $values[0];
					foreach my $type (keys %$support) {
						$xml->startTag('support', 'type' => $type);
						$xml->characters($support->{$type});
						$xml->endTag('support');
					}
				}
			$xml->endTag('branch');
		}
	# annotate all tags other than support for present node
		my $tagnames; @$tagnames = $node->get_all_tags;
		if ($tagnames) {
			my @tags = $node->get_all_tags;
			foreach my $tag (@tags) {
				my @values = $node->get_tag_values($tag);
				my $str = $values[0];
				unless ($tag eq 'support') {
					$xml->startTag('custom');
					$xml->dataElement('name',$tag);
					$xml->dataElement('type',$str->{'type'});
					if (($str->{'type'} =~ /^discrete$/i) || ($str->{'type'} =~ /^double_discrete$/i) || ($str->{'type'} =~ /^multistate$/i)) {
						foreach my $state (@{$str->{'state'}}) {
							$xml->startTag('state');
							$xml->dataElement('state',$state->{'state'});
							$xml->dataElement('probability',$state->{'probability'});
							$xml->endTag('state');
						}
					} elsif ($str->{'type'} =~ /^continuous$/i) {
						$xml->dataElement('value',$str->{'value'});
					} else {
						fatalerror("phyloXML saving of variables of state ",$str->{'type'}," is not supported\n")
					}
					$xml->endTag('custom');
				}
			}
		}
	# process subnodes
		if ($node->each_Descendent) {
			foreach my $subnode (order_subnodes($node->each_Descendent)) {
				node2xml($subnode,$xml);
			}
		}
	# close present node
		$xml->endTag('clade');
}

sub order_subnodes {
	my @in = @_;
	unless (defined $taxorder) {return @in}
	my @arr;
	my $counter; $counter = -1;
	foreach my $node (@in) {
		++$counter;
		my $o;
		if ($node->is_Leaf) {$o = $node->{'_taxorder'}}
		else {my @td = grep {$_->is_Leaf} $node->get_all_Descendents; $o = $td[0]->{'_taxorder'}}
		push @arr,$o.' -- '.$counter;
	}
	@arr = Sort_Table(
		cols      => '2',
		field     => '1',
		sorting   => 'ascending',
		structure => 'csv',
		separator => ' -- ',
		data      => \@arr,
	);
	my @out;
	foreach my $line (@arr) {
		$line =~ /^\d+\s\-\-\s(\d+)$/;
		push @out,$in[$1];
	}
	return @out;
}

###################################################################################################################################
#   tree processing
###################################################################################################################################

# this subroutine handles node support:
#   if there is node support in the $node->bootstrap, than this is used
#   otherwise, if the $node->id is numeric, it is taken to be node support and annotated as such
#   i'm also duplicating the support information as a custom tag for extra flexibility
sub process_nodeID_support {
	my $trees = shift;
	message("processing node ID values if present\n");
	my ($nodecounter,$treecounter); ($nodecounter,$treecounter) =(0,0);
	foreach my $tree (@$trees) {
		++$treecounter;
		foreach my $node (grep {!($_->is_Leaf)} $tree->get_nodes) {
			++$nodecounter;
			if (defined($node->id) && (length($node->id)>0)) {
				unless ($nodeid eq 'id') {
					unless ($node->id =~ /^[\.\d]+$/) {error("one of the support values was not numeric: ".$node->id."\n");}
					$node->add_tag_value('support',{$nodeid => $node->id});
					$node->{'_id'} = "";
				}
			}
		}
	}
	message("  $nodecounter node IDs belonging to $treecounter trees processed\n");
	return $trees;
}


###################################################################################################################################
#   various subroutines
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

# this subroutine parses the input tree filename to come up with a suitable output filename if the user didn't specify one
#   if there are no dots in the filename, the extension .xml is simply added to the input tree filename
#   otherwise, if the file extension of the input tree filename is a standard tree file extension, this is replaced by .xml
#   if it isn't a standard tree file extension, the extension .xml is added to the input tree filename
sub default_treeout_filename {
	my $infile = shift;
	if ($infile =~ /\./) {
		if ($infile =~ /(.+)\.([^\.]+)/) {
			my ($prefix,$suffix) = ($1,$2);
			$suffix = lc $suffix;
			if ($recognized->{'suffixes'}->{$suffix}) {
				return $prefix.'.xml';
			} else {
				return $infile.'.xml';
			}
		}
	} else {
		return $infile.'.xml';
	}
}

# this subroutine does some basic tests to see if trees meet some minimum requirements
sub testtrees {
	my $trees = shift;
	my $treecounter; $treecounter = 0;
	foreach my $tree (@$trees) {
		++$treecounter;
		my $hash;
		foreach my $node ($tree->get_leaf_nodes) {
			if ($hash->{$node->id}) {fatalerror("Tree number $treecounter has multiple occurences of the same taxon name. This is not allowed.\n  Taxon:   ".$node->id."\n")}
			$hash->{$node->id} = 1;
		}
	}
}

sub remove_leading_dashes_and_zeros {
	my $in = shift;
	my $out; $out = $in;
	while ($out =~ /^[0\-]/) {$out = substr($out,1);}
	return $out;
}

sub convert_from_scientific_format {
	my $datum = shift;
	$datum =~ /([\.\d]+)e(\-*\d+)/i;
	my $base = $1; my $exp = $2; my $dir;
	if ($exp =~ /\-/) {$dir = 'd';} else {$dir = 'm';}
	$exp = remove_leading_dashes_and_zeros($exp);
	for (my $i = 0; $i < $exp; ++$i) {
		if ($dir eq 'd') {$base /= 10;} else {$base *= 10;}
	}
	return $base;
}

sub is_probability {
	my $datum = shift;
	unless ($datum =~ /[\d\.\e\-]+/i) {return 0;}
	if ($datum =~ /e/i) {
		unless ($datum =~ /[\d\.]+e\-*\d+/i) {return 0;}
	}
	if ($datum > 1 || $datum < 0) {return 0;}
	return 1;
}