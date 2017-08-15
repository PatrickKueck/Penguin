#!/usr/bin/perl
use strict;
use Encode;
use warnings;

# v1.34 -> v1.35 bug in generating quartets from clade-definition taxa infile order (line 1628,1629)
# v1.35 -> v1.36 bug in -p definition. If -p undefined, wrong parameter value adapted from next line in %default_parameter_setting (line 202), singleton sort command in phyquart analysis
# v1.36 -> v1.37 revised file print options
# v1.37 -> v1.38 revised 4basal infile print out
# v1.38 -> v1.39 source code improvement in subroutine &perform_quartet_analyses
# v1.39 -> v1.391 adapted to new print function of python 3.0 due to the generation of the p4 inputfile: print commands are now in brackets
# v1.391-> v1.392 bug fixed due to p4 inputfile using python 2.0: removement of incorrect bracket prints in p4 output file
# v1.40 -> v1.41 Implementation of k optional normalization scheme of identified split scores
############################################## START GLOBAL SETTING ###########################################

#######################
# # Scriptname and Version Number must be separated by 1 Underscore.
# No other underscore allowed due version number extraction in subroutine &preface
# UPDATED Versions:
# from 1.2beta -> 1.3beta : Script extension from only nucleotide data to nucleotide and amino acid data (Mainly P4 input file setup changes)
# from 1.3beta -> 1.4beta : Implementation of conmpletely spin divided multiple clan analysis (-r option), Revised Graphic Output, New AA substitution models, Spin dependent split signal of each single quartet analysis normalised by the lowest value (set to zero) given all 6 spin scores (previously: lowest of the three best spin topologies normalised to zero)
my $scriptname	= 'Penguin_v1.41beta' ;
#######################

#######################
# Author Info
my %info_of_author			=(
								'name'			=>	'Patrick Kueck'			,
								'institute'		=>	'ZFMK Bonn, GER'		,
								'email'			=>	'pkueck[at]web.de'		,
								'homepage'		=>	'www.zfmk.de'			,
								'lastUpdate'	=>	'March, 2017'			,
								'p4'			=>	'http://p4.nhm.ac.uk/'	,
) ;

#######################
# colors used for graphical split und barplot svg print out
# only three colors are allowed to be not printed out
my %grafik_colours			= (
								'0'							=>	'DodgerBlue'	,	# split and brachart color of     best supported quartet topology -> blue 30,144,255
								'1'							=>	'Tomato'		,	# split and brachart color of 2nd best supported quartet topology -> red 255,99,71
								'2'							=>	'LimeGreen'		,	# split and brachart color of 3rd best supported quartet topology -> green 191,255,0
								#'triangle'					=>	'Cornsilk'		,	# triangle background color 255,248,220
								'triangle'					=>	'WhiteSmoke'	,	# triangle background color
								#'triangle'					=>	'Seashell'		,	# triangle background color
								#'back'						=>	'AliceBlue'		,	# triangle background color
								'back'						=>	'Cornsilk'		,	# triangle background color
								'header'					=>	'RoyalBlue'		,	# header color
								'barleft'					=>	'WhiteSmoke'	,	# barplot color of lower supported spin dependent trees
								'qdots'						=>	'DodgerBlue'	,	# triangle dot color for single quartet presentation
								'clan1'						=>	'DodgerBlue'	,	# color trinagle taxon dots of first clan taxa
								'clan2'						=>	'BlueViolet'	,	# color trinagle taxon dots of second clan taxa
								'clan3'						=>	'LimeGreen'		,	# color trinagle taxon dots of third clan taxa
								'clan4'						=>	'brown'			,	# color trinagle taxon dots of fourth clan taxa
) ;

#######################



#######################
# titel text for each graphic svg output file, coded by hash assigned string
my %graphic_title			=(
								'tritax'					=>	"Penguin Mean Split Support per Taxa"		,	# Title for SVG Triangle output for splitsupport of single taxa analyses
								'trimed'					=>	"Penguin Median Split Support per Taxa"		,	# Title for SVG Triangle output for splitsupport of single taxa analyses
								'triall'					=>	"Penguin Split Support per Quartet"			,	# Title for SVG Triangle output for splitsupport of single quartet analyses
								'splimean'					=>	"Penguin Mean Split Support (Overall)"		,	# Title for SVG Split output for mean splitsupport of single taxa analyses
								'splimedi'					=>	"Penguin Median Split Support (Overall)"	,	# Title for SVG Split output for median splitsupport of single taxa analyses
								'barDs'						=>	"Penguin Spin Score Ranking Distribution"	,	# Title for SVG Ranking Distribution of Each Quartet Tree given single spin quartet analysis
								'barDb'						=>	"Penguin Tree Score Ranking Distribution"	,	# Title for SVG Ranking Distribution of Each Quartet Tree given best spin quartet analysis
) ;
#######################



#######################
# output filenames , coded by hash assigned string
# Additional print out possible for:
# txtssupport, barN, txtbsupport, txtaddprint
# comment out hash code assigned 'none' value line to activate additional print out
my %text_title				=(
								'txtrejected'				=>	$scriptname.'_qsingle_rej_pos.txt'				,	# Text list result file about rejected quartets
								'txtqsupport'				=>	$scriptname.'_qsingle_spl_sup.txt'				,	# Text list presenting single quartet split scores (used for 4basal)
								'txttsupport'				=>	$scriptname.'_tsingle_mean_sup.txt'				,	# Text list present mean split scores for each sequence and for each quartet topology
								'txttsupmedi'				=>	$scriptname.'_tsingle_median_sup.txt'			,	# Text list present median split scores for each sequence and for each quartet topology
								'txtssupport'				=>	$scriptname.'_ssingle_spl_sup.txt'				,	# Text list present mean split scores for each subclade relationship of given quartets
								'txtinfo'					=>	$scriptname.'_result_info.txt'					,	# Main info result file
								'tmpinfo'					=>	$scriptname.'_result_info.tmp'					,	# Temporary main info result file
								'triall'					=>	$scriptname.'_triangle_qsingle_split_sup.svg'	,	# Triangle output for splitsupport of single quartet analyses
								'tritax'					=>	$scriptname.'_triangle_tsingle_mean_sup.svg'	,	# Triangle output for splitsupport of single taxa analyses
								'trimed'					=>	$scriptname.'_triangle_tsingle_median_sup.svg'	,	# Triangle output for splitsupport of single taxa analyses
								'splimean'					=>	$scriptname.'_split_qmean_sup.svg'				,	# Split output for mean splitsupport of single taxa analyses
								'splimedi'					=>	$scriptname.'_split_qmedian_sup.svg'			,	# Split output for median splitsupport of single taxa analyses
								'barDs'						=>	$scriptname.'_spin_tree_order_distribution.svg'	,	# SVG Ranking Distribution of Each spin Quartet Tree given single quartet analysis
								'barDb'						=>	$scriptname.'_best_tree_order_distribution.svg'	,	# SVG Ranking Distribution of Each best spin or overall Quartet Tree given single quartet analysis
								'txtDb'						=>	$scriptname.'_best_tree_order_distribution.txt'	,	# Text print of quartets, sorted by best, second best and worst support
								'txtaddprint'				=>	$scriptname.'_detailed_split_calc.txt'			,	# Additional print out for each single quartet analyses, presenting detailed calculation steps

								# Comment out for corresponding result print out (as defined above)
								#'txtrejected'				=>	'none'										,	# Text list result file about rejected quartets
								#'txtqsupport'				=>	'none'										,	# Text list presenting single quartet split scores (used for 4basal)
								#'txttsupport'				=>	'none'										,	# Text list present mean split scores for each sequence and for each quartet topology
								#'txttsupmedi'				=>	'none'										,	# Text list present median split scores for each sequence and for each quartet topology
								'txtssupport'				=>	'none'										,	# If none, no print of 4Basal input file. For print comment out this line
								#'txtinfo'					=>	'none'										,	# Main info result file
								#'triall'					=>	'none'										,	# Triangle output for splitsupport of single quartet analyses
								#'tritax'					=>	'none'										,	# Triangle output for splitsupport of single taxa analyses
								#'trimed'					=>	'none'										,	# Triangle output for splitsupport of single taxa analyses
								#'splimean'					=>	'none'										,	# Split output for mean splitsupport of single taxa analyses
								#'splimedi'					=>	'none'										,	# Split output for median splitsupport of single taxa analyses
								'barDs'					=>	'none'										,	# SVG Ranking Distribution of Each spin Quartet Tree given single quartet analysis
								'barDb'					=>	'none'										,	# SVG Ranking Distribution of Each best spin or overall Quartet Tree given single quartet analysis
								#'txtDb'					=>	'none'										,	# If none, no print of taxa combinations found for best and second best split support. For print comment out this line
								'txtaddprint'				=>	'none'										,	# If none, no print of detailed split calculations. For print comment out this line
);
#######################



#######################
# Hash lists of possible P4 ML models, assigned to corresponding sequence type
my %hoh_ml_models_of_seq_type	=(
								'nuc'						=>	{
																	'JC'		=>	'nuc'	,
																	'F81'		=>	'nuc'	,
																	'K2P'		=>	'nuc'	,
																	'HKY'		=>	'nuc'	,
																	'GTR'		=>	'nuc'	,
								}	,

								'aa'						=>	{
																	'd78'		=>	'aa'	,
																	'jtt'		=>	'aa'	,
																	'wag'		=>	'aa'	,
																	'mtrev24'	=>	'aa'	,
																	'lg'		=>	'aa'	,
																	'blosum62'	=>	'aa'	,
								}
) ;
#######################



#######################
# Standard values of grafik output considering subtree circles, terminal branches, and internal node circle radius
my %graphic_setup			=(
								'stroke_width_branches'			=>	0.25		,	#
								'stroke_color_branches'			=>	'grey'		,	#
								'axis_stroke_opacity'			=>	'1.000000'	,	#
								'axis_fill_opacity'				=>	'0.000000'	,	#
								'axis_fill_color'				=>	"#000000"	,	#
								'axis_stroke_color'				=>	"#000000"	,	#
								'text_font_size'				=>	30			,	#
								'text_stroke'					=>	'black'		,	#
								'text_stroke_2'					=>	'orange'	,	#
								'text_font_weight'				=>	'normal'	,	#
								'text_stroke_width'				=>	0			,	#
								'text_fill'						=>	'black'		,	#
								'text_fill_2'					=>	'orange'	,	#
								'text_font_family'				=>	'Arial'		,	#
								'text_curser_fill'				=>	'RoyalBlue'	,	#
								'scale_color'					=>	'black'		,	#
								'scale_stroke_width'			=>	1			,	#
								'scale_std_value'				=>	0.1			,	#
								'scale_font_size'				=>	10			,	#
								'legend_font_size'				=>	10			,	#
								'legend_font_weight'			=>	'normal'	,	#
								'legend_stroke_width'			=>	1			,	#
								'stroke_width_branches'			=>	2			,	#
								'split_factor'					=>	200			,	#
								'font_size'						=>	10			,	#
								'font_family'					=>	'Arial'		,	#
								'stroke-width'					=>	1			,	#
								'circle_radius'					=>	1.25		,	#
								'circle_radius_mean'			=>	1.25		,	#
								'circle_stroke_color'			=>	'blue'		,	#
								'circle_stroke_color_mean'		=>	'red'		,	#
								'circle_stroke_width'			=>	0.1			,	#
								'circle_fill_color'				=>	'blue'		,	#
								'circle_fill_color_mean'		=>	'red'		,	#
								'circle_curser_fill'			=>	'black'		,	#
								'bar_stroke_width'				=>	1			,	#
								'bar_stroke_color'				=>	"#000000"	,	#
								'bar_fill_color'				=>	"#BEBEBE"	,	#
								'bar_stroke_opacity'			=>	'1.000000'	,	#
								'bar_stroke_fill'				=>	'1.000000'	,	#
								'raster_stroke_width'			=>	0.5			,	#
								'raster_stroke_color'			=>	"#000000"	,	#
								'raster_stroke_opacity'			=>	'0.700000'	,	#
								'raster_fill_color'				=>	"#000000"	,	#
								'raster_fill_opacity'			=>	'0.000000'	,	#
);
#######################



#######################
# Resultfolders
my %name_of_subfolder		=(
								'svg'						=>	$scriptname.'_SVG'	,	# Subfolder SVG output
								'txt'						=>	$scriptname.'_TXT'	,	# Subfolder TXT output

								# If all result files have to be printed to the script home folder comment out following two lines
								#'svg'						=>	'.'	,	# Subfolder SVG output
								#'txt'						=>	'.'	,	# Subfolder TXT output
);
#######################



#######################
# print header
print	"\n\t___________________________________________________",
		"\n\t---------------------------------------------------",
		"\n\t",
		"\n\t", $scriptname ,
		"\n\tSynapomorph Split Signal Identification in Quartets",
		"\n\t",
		"\n\t---------------------------------------------------\n"
;
#######################

############################################## END GLOBAL SETTING ###########################################






############################################## START ARGV READ IN ###########################################

#######################
# Define default parameter options
# store parameter options & infile name in %$href_par_set_of_opt
# key: parameter; value: set up value
my %default_parameter_setting = (
						'N' => 0		,	# def: 0		indel,amb handling_msa					skip from 0 (reject in single quartets, default) to 1 (overall alignment exclusion)
						'c' => 0		,	# def: 0		character definition					skip from 0 (states, default) to 1 (classes)
						'a'	=> 0.5		,	# def: 0.5		<float>	(alpha value)					skip from '0.5' (default) to defined value
						'I'	=> 0.3		,	# def: 0.3		<float>	(pinv value)					skip from '0.3' (default) to defined value
						'm'	=> 'GTR'	,	# def: GTR		<string> (subst model)					skip from 'GTR' (default) to defined value (JC & GTR allowed)
						'x'	=> 'wag'	,	# default model for amino acid data. If sequence data is identified as aa data, waq will be new default value for option 'm'
						'l' => 10000	,	# def: 10000	<integer> (Max Quartet limit)			skip from 10000 (default) to defined value
						'M'	=> 500		,	# def: 500		<integer> (Min. N. of L seq)			skip from 500 (default) to defined value
						'p'	=> ''		,	# undef			<string> (clade infile)					skip from undef to input clade definition file (txt)
						'i' => ''		,	# undef			<string> (MSA infile)					skip from undef to msa input file (.fas or .phy)
						'u' => 0		,	# def: 0		reducing N quartets randomly			skip from 0 -> ask to draw quartets randomly if N possible quartets is higher as N allowed quartets to 1 -> just do it without questioning (option has been implemented to avoid terminal questioning in process pipelines)
						'r' => 0		,	# def: 0		consideration spin directions			skip from 0 -> no separate spin directive analyses of multiple taxon analyses to 1 -> spin directive analyses
						'k' => 0		,	# def: 0		step of split score malization			skip from 0 -> mormalization by setting the lowest of the three best spin dependent scores to zero (in subroutine &start_phyquart) =>1 : normalization by setting the lowest of all 6 spin dependent split scores to zero
					) ;

my %parameter_setting_of_option	=	%default_parameter_setting ;
#######################



#######################
# ARGV -> READ IN script options, clade definition infile, and msa infile name
# All option values are stored in %parameter_setting_of_option
# key: option command, value: option value
&argv_handling(
								\$text_title{tmpinfo}			,	#	Temporary main info result file					In -> Defined, OUT -> Unchanged
								\$scriptname					,	#	Scriptname										In -> Defined, OUT -> Unchanged
								\%parameter_setting_of_option	,	#	key: parameter code; value: parameter value		In -> Defined, OUT -> Changed
								\%info_of_author				,	#	key: author code; value: info					In -> Defined, OUT -> Unchanged
								\%default_parameter_setting		,	#	key: parameter code; value: default value		In -> Defined, OUT -> Unchanged
) ;
#######################

############################################## END ARGV READ IN ##############################################






############################################## START INFILE READIN ###########################################
# MSA Data readin, data check, exclusion of defined sites and recoding of quartet pattern

my %data_of_infile_property ;			# key: property; value : property value
# $data_of_infile_property{type}		= sequence tpye ('nuc' || 'aa')
# $data_of_infile_property{length}		= original sequence length (e.g.: 10000)
# $data_of_infile_property{remain_pos}	= Number of remaining positions after exclusion of defined sites

my %unreduced_sequence_of_taxa	= () ;
# key: taxonname;
# value sequence

my %hol_sequence_states_of_taxon	= () ;
# key1: taxonname;
# key2: seqposition number
# value site state


#######################
# MSA Infile Check (incl. Data Read IN)
&msa_infile_check (

								\%parameter_setting_of_option	,	#	key: option command, value: option value												In -> Defined, OUT -> Unchanged
								\%data_of_infile_property		,	#	key: property; value : property value													In -> Undefined, OUT -> defined
								\%unreduced_sequence_of_taxa	,	#	key: taxonname; value sequence															In -> Undefined, OUT -> defined
								\%hol_sequence_states_of_taxon	,	#	key1: taxon name; key2: state position number; value: sequence state at that position	In -> Undefined, OUT -> defined
								\$text_title{tmpinfo}			,	#	temporary info outfile name																In -> Defined, OUT -> Unchanged
								\$scriptname					,	#	scriptname																				In -> Defined, OUT -> Unchanged
								\%info_of_author				,	#	key: author code; value: info															In -> Defined, OUT -> Unchanged
								\%default_parameter_setting		,	#	key: parameter code; value: default value												In -> Defined, OUT -> Unchanged
) ;
#######################



#######################
# P4 ML model definition check, defined model must be congruent to identified sequence type
#unless ( $hoh_ml_models_of_seq_type{$data_of_infile_property{type}}{$parameter_setting_of_option{m}}	eq $data_of_infile_property{type} ){
unless ( $hoh_ml_models_of_seq_type{$data_of_infile_property{type}}{$parameter_setting_of_option{m}} ){

	my	 $outinfo	=	$text_title{tmpinfo} ;
	open OUTinfo,	">>$outinfo" or die "OUTFILE-ERROR: Cannot open info outfile ", $outinfo, "in subroutine in main routine!\n" ;

	print 			"\n\tPARAMETER-WARNING: Defined P4 ML substituion model ", $parameter_setting_of_option{m}, " inappropriate for ", $data_of_infile_property{type}, " data!\n";
	print OUTinfo	"\n\tPARAMETER-WARNING: Defined P4 ML substituion model ", $parameter_setting_of_option{m}, " inappropriate for ", $data_of_infile_property{type}, " data!\n";

	if		( $data_of_infile_property{type} eq 'nuc'	){	$parameter_setting_of_option{m}	=	$default_parameter_setting{m} }
	elsif	( $data_of_infile_property{type} eq 'aa' 	){	$parameter_setting_of_option{m}	=	$default_parameter_setting{x} }

	print 			"\n\tChanged to default model ", $parameter_setting_of_option{m}, "\n" ;
	print OUTinfo	"\n\tChanged to default model ", $parameter_setting_of_option{m}, "\n" ;

	close OUTinfo
}
#######################



#######################
# Clade Infile Check (incl. Data Read IN)
my %hol_taxa_of_subgroup = &clade_check(

								\%unreduced_sequence_of_taxa	,	#	key: taxonname; value sequence					In -> Defined, OUT -> Unchanged
								\%parameter_setting_of_option	,	#	key: option command, value: option value		In -> Defined, OUT -> Unchanged
								\$text_title{tmpinfo}			,	#	temporary info outfile name						In -> Defined, OUT -> Unchanged
								\$scriptname					,	#	scriptname										In -> Defined, OUT -> Unchanged
								\%info_of_author				,	#	key: author code; value: info					In -> Defined, OUT -> Unchanged
								\%default_parameter_setting		,	#	key: parameter code; value: default value		In -> Defined, OUT -> Unchanged
) ;

## test print
#for my $subclade ( sort keys %hol_taxa_of_subgroup ){
#
#	print "\n\tSubclade: ", $subclade, "\n\tTaxa: " ;
#	my @taxa = @{$hol_taxa_of_subgroup{$subclade}};
#	for my $taxon ( @taxa ){ print " $taxon" }
#	print "\n"
#}
############################################## END INFILE READIN ###########################################






############################################## START OVERALL SITE EXCLUSION ################################
my %rejected_site_positions ; # key: site position numbers which are rejected in further split analysis
if ( $parameter_setting_of_option{N} == 1 ){
			
			print "\n\tExclude MSA Site Positions..." ;
			
			
			
			#######################
			# Define allowed sequence states
			# all other sequence states are rejected from further split analyses
			my %allowed_states ; # include all allowed site characters
			
			if			( $data_of_infile_property{type} =~ 'nuc'	){ for ( qw/A C G T/ 										)	{ $allowed_states{$_}++ } } # all nuc states allowed
			elsif		( $data_of_infile_property{type} =~ 'aa'	){ for ( qw/A C G T N Y R W S K M D V H I E L Q F P/		)	{ $allowed_states{$_}++ } } # all aa states allowed
			else	{ die "\n\tBUG-ERROR: Cannot assign data exclusion parameter N and A in main routine!\n\tPlease, report BUG to system developer!\n\t" 	}
			#######################
			
			
			
			#######################
			# Exlusion of indel site positions over all alignment sequences if...
			# ...-N = 1 option or -A = 0
			&site_exclusion (
			
								\%unreduced_sequence_of_taxa,		# key: taxonname; value sequence															In -> Defined, OUT -> Unchanged
								\%data_of_infile_property,			# key: property; value : property value														In -> Defined, OUT -> Changed
								\%parameter_setting_of_option,		# key: option command, value: option value													In -> Defined, OUT -> Unchanged
								\%rejected_site_positions,			# key: site position numbers which are rejected in further split analysis					In -> Undefined, OUT -> defined
								\%allowed_states,					# key: allowed character states; value: 1													In -> Defined, OUT -> Unchanged
								\%hol_sequence_states_of_taxon,		# key1: taxon name; key2: state position number; value: sequence state at that position		In -> Defined, OUT -> Unchanged
			) ;
			#######################
			
			
			
			#######################
			# If remaining site positions smaller as allowed sequence length abort analysis
			my $N_rejected_site_positions	=	0 + keys %rejected_site_positions ;
			my $remaining_pos = $data_of_infile_property{length} - $N_rejected_site_positions ;
			if ( $remaining_pos < $parameter_setting_of_option{M} ){
				
				my		$outinfo	=	$text_title{tmpinfo} ;
				open	OUTinfo,	">>$outinfo" or die "OUTFILE-ERROR: Cannot open info outfile ", $outinfo, "in main routine!\n" ;
				print	OUTinfo		"\n\tCannot Analyse Data Set\n\tRemaining Overall Sequence Length ", $remaining_pos, "<", $parameter_setting_of_option{M}, "bp\n\n" ;
				close	OUTinfo		;
				
				print "\n\tCannot analyse data set\n\tRemaining alignment length ", $remaining_pos, "<", $parameter_setting_of_option{M}, "bp\n\n"; &help( \'M', \$scriptname, \%parameter_setting_of_option, \%info_of_author, \%default_parameter_setting )
			}
			#######################
			
			
			
			#######################
			# open info file
			# terminal & file print out of choosen parameter settings
			my	 $outinfo	=	$text_title{tmpinfo} ;
			open OUTinfo, ">>$outinfo" or die "OUTFILE-ERROR: Cannot open info outfile ", $outinfo, "in main routine!\n" ;
			
			print OUTinfo 	"\n\tRejected Overall Alignment Positions:\t", $N_rejected_site_positions ;
			print 			"\n\tRejected Overall Alignment Positions:\t", $N_rejected_site_positions ;
			
			print OUTinfo	"\n\tRemaining Alignment Positions:\t", $remaining_pos, "\n" ;
			print 			"\n\tRemaining Alignment Positions:\t", $remaining_pos, "\n" ;
			
			close OUTinfo ;
			#######################
}
############################################## END OVERALL SITE EXCLUSION #################################






############################################## START QUARTET GENERATION ################################

#######################
# Generate all possible quartets between defined subgroup taxa
my @quartets = &generate_quartets(

	\%hol_taxa_of_subgroup,	# key: subclade name; value: list of subclade assigned taxa		-> IN (not changed)
) ;

## test print
#my $counter = 1 ;
#print "\n\tGenerated Quartets:\n\t" ;
#for ( @quartets ){ print "test Q", $counter, ": ", $_, "\n\t"; $counter++ }
#print "\n" ;
####
#######################



#######################
# Check maximum Number of Quartets
if ( @quartets > $parameter_setting_of_option{l} ){
	
	my $answer	= 0 ;
	
	unless ( $parameter_setting_of_option{u} == 1 ){
		
		print	"\n\t---------------------------------------------------\n" ,
				"\n\tWARNING: N Quartets Succeed Max. Limit of Allowed Quartets\n\t",
				"Press <q>	<enter>\tto Quit\n\t",
				"Press <r>	<enter>\tto Choose ", $parameter_setting_of_option{l} ," Quartet(s) Randomly\n\t",
				"Press		<enter>\tto Proceed\n\n\tCommand: ";
		
		chomp ( $answer = <STDIN> ) ;
		
		if		( $answer =~ /q/i )	{ print "\n\tQuit ", $scriptname ,"!\n\n"; exit }
	}
	
	if	(	(	$answer							=~	/r/i	)	||
			(	$parameter_setting_of_option{u}	==	1		)	){
		
		print "\n\t", $parameter_setting_of_option{l} ," Quartet(s) Drawn by Random!\n\n";
		
		my $range 				=	@quartets-1 ;
		my @new_quartets						;
		my $quartet_counter		=	1 			;
		my %seen_new_quartet					;
		
		while ( $quartet_counter <= $parameter_setting_of_option{l} ){
			
			my $random_number = int(rand($range));
			unless ( $seen_new_quartet{$quartets[$random_number-1]} ){
				
				push @new_quartets, $quartets[$random_number-1] ;
				$seen_new_quartet{$quartets[$random_number-1]}++;
				$quartet_counter++
			}
		}
		
		@quartets = @new_quartets ;
		( %seen_new_quartet, @new_quartets, $quartet_counter, $range ) = () ;
	}
	else	{ print "\n\tProceed ", $scriptname ,"!\n\n" }
}
#######################



#######################
# Print terminal info about generated quartets
my $counter = 1 ;
print "\n\tGenerated Quartet(s):\n\t" ;

for ( @quartets ){  print "q", $counter, ":\t", $_, "\n\t"; $counter++ } $counter = () ;

print "---------------------------------------------------\n" ;
#######################

############################################## END QUARTET GENERATION ################################






############################################## START QUARTET ANALYSES ################################

#######################
# Start single quartet analysis
&perform_quartet_analyses(
							\@quartets							,	# list of generated quartets (each equal $taxon1.":".$taxon2.":".$taxon3.":".$taxon4)		In -> Defined, OUT -> Unchanged
							\%unreduced_sequence_of_taxa		,	# key: taxonname; value sequence															In -> Defined, OUT -> Unchanged
							\%data_of_infile_property			,	# key: property; value : property value														In -> Defined, OUT -> Unchanged
							\%parameter_setting_of_option		,	# key: option command, value: option value													In -> Defined, OUT -> Unchanged
							\%rejected_site_positions			,	# key: site position numbers which are rejected in further split analysis					In -> Defined, OUT -> Changed
							\%hol_taxa_of_subgroup				,	# key: subclade name; value: list of subclade assigned taxa									In -> Defined, OUT -> Unchanged
							\%hol_sequence_states_of_taxon		,	# key1: taxon name; key2: state position number; value: sequence state at that position		In -> Defined, OUT -> Unchanged
							\%grafik_colours					,	# key: topology support (0,1, or 2); value: assigned svg color								In -> Defined, OUT -> Unchanged
							\%graphic_title						,	# key: graphic outfilename; value: graphic title											In -> Defined, OUT -> Unchanged
							\%text_title						,	# key: output code of txt files; value: output filename										In -> Defined, OUT -> Unchanged
							\%graphic_setup						,	# key: svgoption; value: option setup														In -> Defined, OUT -> Unchanged
							\%name_of_subfolder					,	# key: subfolder code; value: subfolder name												In -> Defined, OUT -> Unchanged
						) ;
#######################

############################################## END QUARTET ANALYSES ################################






############################################## START SCRIPT CLOSING ################################

#######################
# Rename temporary inputfile
unless  ( $text_title{txtinfo} eq 'none' ){ rename ( $text_title{tmpinfo}, $name_of_subfolder{txt}."/".$text_title{txtinfo} ); print	"\n\tPrint TXT ", $text_title{txtinfo} }
else	{ unlink $text_title{tmpinfo} }
#######################



#######################
# print end
print	"\n\t---------------------------------------------------",
		"\n\n\t" , $scriptname, " Analysis Completed !",
		"\n\tAll Results Printed to ", $scriptname, " Homefolder",
		"\n\t",
		"\n\tThank You & Good Bye !" ,
		"\n\t",
		"\n" ,
		"               __                       \n" ,
		"            -=(o '.                     \n" ,
		"               '.-.\\                   \n" ,
		"               /|  \\\\                 \n" ,
		"               '|  ||                   \n" ,
		"                _\\_):,_                  " ,
		"\n\t___________________________________________________",
		"\n\t---------------------------------------------------\n\n",
;
#######################

exit ;
############################################## END SCRIPT CLOSING ################################






################################################################################################ SUBROUTINES ################################################################################################################################

sub help{

	my	$sref_command			= $_[0]	;	#	help print of specific command. if command = 'h' print everything		In -> Defined, OUT -> Unchanged
	my	$sref_script_name		= $_[1]	;	#	scriptname																In -> Defined, OUT -> Unchanged
	my	$href_setup_of_option	= $_[2]	;	#	key: parameter; value: set up value										In -> Defined, OUT -> Unchanged
	my	$href_authorinfo		= $_[3]	;	#	key: info code; value: info												In -> Defined, OUT -> Unchanged
	my	$href_default_set		= $_[4]	;	#	key: parameter code; value: default value								In -> Defined, OUT -> Unchanged

	my $website_P4						= $href_authorinfo->{p4} ;

	if		( $$sref_command eq 'h' ){

		print	"\n\tUsage (MAC || Linux):" ,
				"\n\tperl ", $$sref_script_name, " -i [msa infile] -p [clan infile] -[command]",
				"\n",
				"\n\tUsage (Windows):" ,
				"\n\t", $$sref_script_name, " -i [msa infile] -p [clan infile] -[command]",
				"\n\t",
				"\n\t[msa infile]\t:\tMultiple Sequence Alignment Infile (.phy || .fas)",
				"\n\t[clan infile]\t:\tClan Definition Infile (.txt)",
				"\n\t[command]\t:\tAdditional Parameter Options...",
				"\n\t",
				"\n\t [M] <integer>\t:\tMinimum Allowed Sequence Length For Quartets",
				"\n\t [l] <integer>\t:\tMaximum Number Generated Quartets",
				"\n\t [m] <string>\t:\tML Substitution Model (NUC || AA)",
				"\n\t [a] <float>\t:\tAlpha Start Value Of ML Estimation",
				"\n\t [I] <float>\t:\tpINV Start Value Of ML Estimation",
				"\n\t [N]\t\t:\tExclude Indels Of The Entire Alignment",
				"\n\t [c]\t\t:\tTranslate Characters To Class States",
				"\n\t [r]\t\t:\tAdditive Root Info In Multi-Taxon Analyses",
				"\n\t [k]\t\t:\tFinal Quartet Score Weighting Scheme",
				"\n\t [u]\t\t:\tElimination Of Script Queries ",
				"\n\t",
				"\n\tFor A Detailed Help Menu About Parameters Type:",
				"\n\tperl ", $$sref_script_name, " -h [command], e.g.:",
				"\n\tperl ", $$sref_script_name, " -h i",
				"\n\t",
				"\n\tFor License Information Type:",
				"\n\tperl ", $$sref_script_name, " -P",
				"\n\t",
				"\n\t", $$sref_script_name, " Needs The P4 Python Packages",
				"\n\tInstalled. For Download And A Detailed Explanation",
				"\n\tOf P4 Visit: ", $website_P4,
				"\n\t---------------------------------------------------\n\n",
		;
	}

	elsif	( $$sref_command eq 'i' ){

		print	"\n\tMultiple Sequence Alignment Infile:" ,
				"\n\t-----------------------------------",
				"\n\tAllowed Formats:",
				"\n\t- FASTA  (.fas)",
				"\n\t- PHYLIP (.phy) (Strict || Relaxed)",
				"\n\t",
				"\n\tAllowed Sequence Conditions:",
				"\n\t\t- Sequences Of Equal Length",
				"\n\t\t- Nucleotide States",
				"\n\t\t- Amino Acid States",
				"\n\t\t- Ambiguity States",
				"\n\t\t- Indel/GAP States ('-')",
				"\n\t",
				"\n\tAllowed Sequence Names:",
				"\n\t\t- Alphanumeric Signs",
				"\n\t\t- Underscores ('_')",
				"\n\t\t- Only Unique Sequence Names",
				"\n\t",
				"\n\tSpecified Via '-i' Option, e.g.:",
				"\n\t\t\"-i MSA_infile_name.fas\"",
				"\n\t---------------------------------------------------\n\n",
		;
	}

	elsif	( $$sref_command eq 'p' ){

		print	"\n\tClan Definition Infile:" ,
				"\n\t-----------------------------------",
				"\n\tAllowed Format:",
				"\n\t- TEXT (.txt), e.g.:",
				"\n\t",
				"\n\t\t--",
				"\n\t\tClan1,SeqName1,SeqName2,SeqName3...",
				"\n\t\tClan2,SeqNameA,SeqNameB,SeqNameC...",
				"\n\t\tClan3,SeqName4,SeqName5,SeqName6...",
				"\n\t\tClan4,SeqNameD,SeqNameE,SeqNameF...",
				"\n\t\t--",
				"\n\t",
				"\n\tClade Definition:",
				"\n\t\t- Four Clans Must be Defined in Separate Lines",
				"\n\t\t- Clan Defined By First Code In Each Line",
				"\n\t\t- Only Uniquely Defined Clan Codes Allowed",
				"\n\t",
				"\n\tSequence Assignment to Clans:",
				"\n\t\t- Sequence Names of Given MSA (Comma separated)",
				"\n\t\t- Must Be In The Same Line As Assigned Clan Code",
				"\n\t\t- No Whitespace Allowed",
				"\n\t\t- Only Unique Assigned Sequence Names Allowed",
				"\n\t\t- Sequence Names Must Be Identic With MSA (Case Sensitive)",
				"\n\t",
				"\n\tDefinition of Possible Quartet Combinations",
				"\n\tBetween Defined Clan Sequences. Must Not",
				"\n\tBe Specified If MSA Consists Of Only Four",
				"\n\tSequences. Only Alphanumeric Signs and",
				"\n\tUnderscores Allowed.",
				"\n\t",
				"\n\tSpecified Via '-p' Option (If Alignment > 4 Seq.), e.g.:",
				"\n\t\t\"-p clan_infile_name.txt\"",
				"\n\t---------------------------------------------------\n\n",
		;
	}

	elsif	( $$sref_command eq 'a' ){

		print	"\n\tML Alpha Shape Start Parameter:" ,
				"\n\t-----------------------------------",
				"\n\tAllowed Format:",
				"\n\t- Float Number (Default: ", $href_default_set->{a} ,")",
				"\n\t",
				"\n\tML Alpha Shape Parameter Of Rate Heterogeinity.",
				"\n\tStart Parameter For P4 Implemented ML Estimation",
				"\n\tOf Potentially Convergent Evolved Split Pattern",
				"\n\tFrequencies.",
				"\n\t",
				"\n\tNOTE: For ML Pattern Estimation Without ASRV Type:",
				"\n\t\t\"-a 100\"",
				"\n\t",
				"\n\tSpecified Via '-a' Option, e.g.:",
				"\n\t\t\"-a 1.0\"",
				"\n\t---------------------------------------------------\n\n",
		;
	}

	elsif	( $$sref_command eq 'I' ){

		print	"\n\tML pINV Start Parameter:" ,
				"\n\t-----------------------------------",
				"\n\tAllowed Format:",
				"\n\t- Float Number Between 0 and 1 (Default: ", $href_default_set->{I} ,")",
				"\n\t",
				"\n\tML Proportion Of Invariabel Site Estimation.",
				"\n\tStart Parameter For P4 Implemented ML Estimation",
				"\n\tOf Potentially Convergent Evolved Split Pattern",
				"\n\tFrequencies.",
				"\n\t",
				"\n\tSpecified Via '-I' Option, e.g.:",
				"\n\t\t\"-I 0.15\"",
				"\n\t---------------------------------------------------\n\n",
		;
	}

	elsif	( $$sref_command eq 'l' ){

		print	"\n\tMaximum Number Generated Quartets:" ,
				"\n\t-----------------------------------",
				"\n\tAllowed Format:",
				"\n\t- Positive Integer Number (Default: ", $href_default_set->{l} ,")",
				"\n\t",
				"\n\tMaximum Number Of Generated Quartets Between Sequences",
				"\n\tOf Predefined Clans. ", $$sref_script_name, " Generates",
				"\n\tAll Possible Quartets Between Clans. However, The Number",
				"\n\tOf Possible Quartets Increases Strongly With The",
				"\n\tNumber Of Clan Assigned Sequences. ", $$sref_script_name,
				"\n\tNeeds Around 2 Seconds To Analyse A Single Quartet.",
				"\n\tTherefore, The Total Computation Time Can Get Very Long",
				"\n\tIf Number Of Quartets Is High.",
				"\n\t",
				"\n\tSpecified Via '-l' Option, e.g.:",
				"\n\t\t\"-l 5000\"",
				"\n\t",
				"\n\tNOTE: If The Number Of Possible Quartets Exceeds",
				"\n\tThe Maximum Number Of Allowed Quartets,",
				"\n\t", $$sref_script_name, " Stops And Asks About Further ",
				"\n\tQuartet Processing:",
				"\n\t",
				"\n\t\t1: Quit ", $$sref_script_name, ": <q enter>",
				"\n\t\t2: Select N Allowed Quartets Randomly: <r enter>",
				"\n\t\t3: Analyse All Possible Quartets: <enter>",
				"\n\t---------------------------------------------------\n\n",
		;
	}

		elsif	( $$sref_command eq 'u' ){

		print	"\n\tElimination Of Script Queries:" ,
				"\n\t-----------------------------------",
				"\n\tThe '-u' Option Opposes Possible Script Queries",
				"\n\tWhich Would Stop Automatic Process Pipelines.",
				"\n\tThe Maximum Number Of Allowed Quartets Is Drawn",
				"\n\tBy Random If Number Of Possible Quartets Succeeds",
				"\n\tThe Number Of Allowed Quartets.",
				"\n\t",
				"\n\tSpecified Via '-u' Option:",
				"\n\t\t\"-u\"",
				"\n\t---------------------------------------------------\n\n",
		;
	}

	elsif	( $$sref_command eq 'M' ){

		print	"\n\tMinimum Sequence Length for Quartets:" ,
				"\n\t-----------------------------------",
				"\n\tAllowed Format:",
				"\n\t- Positive Integer Number (Default (bp): ", $href_default_set->{M} ,")",
				"\n\t",
				"\n\t", $$sref_script_name, " Analyses Only Quartet",
				"\n\tSequences Given A Sequence Length Greater/Equal",
				"\n\tTo The Defined Minimum Number Of Base Pairs (bp).",
				"\n\t", $$sref_script_name, " Performs Best Using Long",
				"\n\tSequences. Therefore, The Higher The Sequence Length",
				"\n\tThe Better The Split Pattern Estimation For",
				"\n\tGiven Quartet Analyses.",
				"\n\tQuartet Sequence Lengths Below Defined Threshold",
				"\n\tAfter Excluding Undesired Sequence Sites",
				"\n\tAre Rejected From The Overall Quartet Analysis.",
				"\n\tUndesired Sequence Sites Include Ambiguity Site States",
				"\n\tAnd Indel Events (Indel Events If Specified).",
				"\n\t",
				"\n\tSpecified Via '-M' Option, e.g.:",
				"\n\t\t\"-M 2000\"",
				"\n\t---------------------------------------------------\n\n",
		;
	}

	elsif	( $$sref_command eq 'N' ){

		print	"\n\tExclude Indel/Ambiguity Sites Over Complete Alignment:" ,
				"\n\t-----------------------------------",
				"\n\tTwo Possible Options To Deal With Unallowed",
				"\n\tSite Positions ('-'):",
				"\n\t",
				"\n\t\t1: Reject INDEL Sites From The Complete Alignment",
				"\n\t\t2: Reject INDEL Sites From Single Quartets (Default)",
				"\n\t",
				"\n\tExclusion From The Complete Sequence Alignment Can",
				"\n\tBe Specified Via The '-N' Option. With The '-N' Option",
				"\n\t",
				"\n\tNOTE: Rejected Site Position Are Irrevocably Rejected ",
				"\n\tIn All Single Quartet Analyses. Under Default Site Position",
				"\n\tAre Excluded Separately From Scratch In Each Single Quartet",
				"\n\tAnalysis.",
				"\n\t",
				"\n\tSpecified Via '-N' Option:",
				"\n\t\t\"-N\"",
				"\n\t---------------------------------------------------\n\n",
		;
	}

	elsif	( $$sref_command eq 'm' ){

		print	"\n\tML Substitution Model Parameter:" ,
				"\n\t-----------------------------------",
				"\n\tAllowed Format:",
				"\n\t- String\n",
				"\n\t(Default NUC): ", $href_default_set->{m} ,
				"\n\t(Default AA): ", $href_default_set->{x} ,
				"\n\t",
				"\n\tParameter For P4 Implemented ML Estimation",
				"\n\tOf Potentially Convergent Evolved Split Pattern",
				"\n\tFrequencies. Following Models Are Implemented...",
				"\n\t",
				"\n\t- ML Substitution Models For Nucleotide Data:",
				"\n\t",
				"\n\t\t1: GTR",
				"\n\t\t2: HKY",
				"\n\t\t3: K2P",
				"\n\t\t4: F81",
				"\n\t\t5: JC" ,
				"\n\t",
				"\n\t- ML Substitution Models For Amino Acid Data:",
				"\n\t",
				"\n\t\t1: wag",
				"\n\t\t2: jtt",
				"\n\t\t3: d78",
				"\n\t\t4: mtrev24",
				"\n\t\t5: lg",
				"\n\t\t6: blosum62",
				"\n\t",
				"\n\tBE AWARE of CASE SENSITIVE MODEL DEFINITION!",
				"\n\t",
				"\n\tSpecified Via '-m' Option, e.g.:",
				"\n\t\t\"-m JC\"",
				"\n\t---------------------------------------------------\n\n",
		;
	}

	elsif	( $$sref_command eq 'c' ){

		print	"\n\tTranslate NUC States to Class States:" ,
				"\n\t-----------------------------------",
				"\n\tWith The '-c' Option All Nucleotide Characters",
				"\n\tAre Translated And Analysed As Purines And Pyrimidines,",
				"\n\tWhich Is A More Conservative Split Analyses",
				"\n\tDue To The Reduction Of Characters From Four To Two.",
				"\n\t",
				"\n\tSpecified Via '-c' Option:",
				"\n\t\t\"-c\"",
				"\n\t---------------------------------------------------\n\n",
		;
	}

	elsif	( $$sref_command eq 'r' ){

		print	"\n\tAdditive Root Info in a Multi-Taxon Analysis:" ,
				"\n\t-------------------------------------------",
				"\n\tWith the '-r' option split support values",
				"\n\tof single quartets are evaluated separately",
				"\n\tin multi-taxon analyses due to different root",
				"\n\tdirective assumptions. Final split support",
				"\n\trefers to the in total best supporting spin",
				"\n\tassumption of each of the three possible",
				"\n\tquartet topologies and not on a spin mixture",
				"\n\tof best supporting scores of single quartets.",
				"\n\t",
				"\n\tWithout the '-r' option, best split support values",
				"\n\tof each quartets of a multi-taxon analysis",
				"\n\tare summarised for each possible quartet tree",
				"\n\tindependent of given spin assumptions.",
				"\n\t",
				"\n\tSpecified via '-r' option:",
				"\n\t\t\"-r\"",
				"\n\t---------------------------------------------------\n\n",
		;
	}

	elsif	( $$sref_command eq 'k' ){

		print	"\n\tFinal Quartet Weighting Scheme:" ,
				"\n\t-------------------------------------------",
				"\n\tFor each quartet tree best split scores of each",
				"\n\tof the three quartet relationships are weighted",
				"\n\tby setting the lowest of the best spin support",
				"\n\tscores of the three possible quartet trees to zero.",
				"\n\tWith the -k option split scores of each",
				"\n\tof the three quartet relationships are weighted",
				"\n\tby setting the lowest of the six spin dependent",
				"\n\tsupport scores to zero.",
				"\n\t",
				"\n\tspecified via '-k' option:",
				"\n\t\t\"-k\"",
				"\n\t---------------------------------------------------\n\n",
		;
	}

	elsif	( $$sref_command eq 'P' ){

		&preface( \$$sref_script_name, \%$href_authorinfo )
	}

	else{ die "BUG-ERROR: Cannot assign command ", $$sref_command, " in subroutine &help!\n\tPlease, report BUG to system developer!\n\n" }

	exit;
}

sub preface{

	my	$sref_script_name		= $_[0]	;	#	scriptname						In -> Defined, OUT -> Unchanged
	my	$href_aut_info			= $_[1]	;	#	key: info code; value: info		In -> Defined, OUT -> Unchanged

	my	@script_name_parts		= split "_", $$sref_script_name ;

	print	"\n\tVersion     : ", $script_name_parts[1]											,
			"\n\tLanguage    : PERL"															,
			"\n\tLast Update : ",	$href_aut_info->{lastUpdate}								,
			"\n\tAuthor      : ",	$href_aut_info->{name}, ", ", $href_aut_info->{institute}	,
			"\n\te-mail      : ",	$href_aut_info->{email}										,
			"\n\tHomepage    : ",	$href_aut_info->{homepage}									,
			"\n\t",
			"\n\tThis program is free software; you can distribute it ",
			"\n\tand/or modify it under the terms of the GNU General Public ",
			"\n\tLicense as published by the Free Software Foundation ; ",
			"\n\teither version 2 of the License, or (at your option) any ",
			"\n\tlater version.",
			"\n\t",
			"\n\tThis program is distributed in the hope that it will be",
			"\n\tuseful, but WITHOUT ANY WARRANTY; without even the",
			"\n\timplied warranty of MERCHANTABILITY or FITNESS FOR A",
			"\n\tPARTICULAR PURPOSE. See the GNU General Public License for",
			"\n\tmore details.",
			"\n\t",
			"\n\tYou should have received a copy of the GNU General Public",
			"\n\tLicense along with this program; if not, write to the Free",
			"\n\tSoftware Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139,",
			"\n\tUSA.",
			"\n\t---------------------------------------------------\n\n"
	;

	exit;
}

sub argv_handling{

	# &argv_handling
	# my %parameter_setting_of_option = &argv_handling( \$info_outfile_name ) ;

	#######################
	# store name of info outfile as string referenz for
	# print out of choosen parameter setting
	my	$sref_inf_out			= $_[0]	;	#	temporary info print outfile name				In -> Defined, OUT -> Unchanged
	my	$sref_script			= $_[1]	;	#	scriptname										In -> Defined, OUT -> Unchanged
	my	$href_par_set_of_opt	= $_[2]	;	#	key: parameter code; value: parameter value		In -> Defined, OUT -> Changed
	my	$href_info_author		= $_[3]	;	#	key: author code; value: info					In -> Defined, OUT -> Unchanged
	my	$href_def_par_setting	= $_[4]	;	#	key: parameter code; value: default value		In -> Defined, OUT -> Unchanged
	#######################

	##############################################
	# Script start under defined options:
	# perl 'scriptname.pl' -i 'quartet_msa_infile.fas' -'parameter_option1' 'value_if_required'  -'parameter_option2'...<enter>
	##############################################

	##############################################
	# parameter options:
	# N						indel handling		skip from 0 (remain, default) to 1 (reject) (overall alignment exclusion)
	# n						indel handling		skip from 1 (reject, default) to 0 (remain) (single exclusion for each 4 taxon quartet)
	# p	<string>			clade infile		skip from '' (undefined, default) to infile name 		# .txt allowed !
	# c						state definition	skip from 0 (states, default) to 1 (classes)
	# i	<string>			msa infile			skip from '' (undefined, default) to infile name 		# .phy .fas .fasta .FASTA allowed !
	# a	<float>				alpha value			skip from '1.0' (default) to defined value
	# I	<float>				pinv value			skip from '0.3' (default) to defined value
	# m	<string>			subst model			skip from 'GTR' (default) to defined value (JC & GTR allowed)
	# M	<number>			Min. N. of L seq	skip from 500 (default) to defined value
	# l <number>			Max Quartet limit	skip from 10000 (default) to defined value
	# r 					spin directions		skip from 0 (default) -> don't analyse spin directions separately to 1 -> analyse spin directions separately in multiple taxon analyses
	# h						help menu			skip to help menu
	# P						Preface				skip to preface text
	##############################################
	
	##############################################
	# help menu options
	# -h					->	general help menu
	# -h -parameter option	->	detailed help menu of parameter option
	##############################################
	
	#####################################################################
	## START ARGV handling
	## e.g.: 'perl m-ice.pl -i 'infile_name' -option value <enter>
	# multiple options possible, some options allow defined values others not
	my ( $commandline ) = join "", @ARGV ;
	
	if ( $commandline ){
		
		#######################
		# Check opion commands
		# Assign user defined settings to %par_set_of_opt
		$commandline					=~	s/ |\s+//g ;
		my		@commands				=	split "-", $commandline ;
		shift	@commands;							# Remove first empty element due to split "-", which is the first sign in the line
		my		%taxon_counter	= () ;	# key: taxonname; value counter of taxonname
		
		for my $command ( @commands ){
			
			my	@command_signs	= split	"", $command ;
					$command				= shift	    @command_signs ;
			my	$parameter			= join	"", @command_signs ;
			
			if 		( (	$command =~ /^p/					) && ( $parameter =~	/\.txt$/																	) )	{ $href_par_set_of_opt->{$command}	= $parameter		}
			elsif 	( 	$command =~ /^N$|^c$|^u$|^r$|^k$/	)																										{ $href_par_set_of_opt->{$command}	= 1					}
			elsif 	( (	$command =~ /^i$/ 					) && ( $parameter =~	/\.fas$|\.phy$|\.fasta$/i													) )	{ $href_par_set_of_opt->{$command}	= $parameter		}
			elsif 	( (	$command =~ /^a$/					) && ( $parameter =~	/^\d+.\d+$/																	) )	{ $href_par_set_of_opt->{$command}	= $parameter		}
			elsif 	( (	$command =~ /^I$/					) && ( $parameter =~	/^\d+.\d+$/																	) )	{ $href_par_set_of_opt->{$command}	= $parameter		}
			elsif 	( (	$command =~ /^m$/					) && ( $parameter =~	/^JC$|^GTR$|^F81$|^K2P$|^HKY$|^d78$|^jtt$|^wag$|^mtrev24$|^lg$|^blosum62$/	) )	{ $href_par_set_of_opt->{$command}	= $parameter		}
			elsif 	( (	$command =~ /^M$|^l$/				) && ( $parameter =~	/^\d+$/																		) )	{ $href_par_set_of_opt->{$command}	= $parameter		}
			
			elsif 	( 	$command =~ /^h$/i					){
			
				if ( $parameter =~ /^p$|^N$|^c$|^i$|^a$|^I$|^m$|^M$|^l$|^P$|^u$|^r$|^k$/ )	{ &help( \$parameter, \$$sref_script, \%$href_par_set_of_opt, \%$href_info_author, \%$href_def_par_setting ) }
				else 																		{ &help( \'h', \$$sref_script, \%$href_par_set_of_opt, \%$href_info_author, \%$href_def_par_setting ) }
			}
			
			elsif 	( 	$command =~ /^P$/	){ &preface( \$$sref_script, \%$href_info_author ) }
			
			else 		{ print "\n\tPARAMETER-ERROR: Unallowed command or value: -", $command, " ",$parameter, "!\n"; &help( \'h', \$$sref_script, \%$href_par_set_of_opt, \%$href_info_author, \%$href_def_par_setting ) }
		}
		#######################
		
		
		
		#######################
		# Check sequence definition
		# if sequences defined by -i store infilename without path extension in $par_set_of_opt{z}
		if		( $href_par_set_of_opt->{i} ){
		
			my	@infile_path				= split( /\/|\\/, $href_par_set_of_opt->{i} ) ;
				$href_par_set_of_opt->{z}	= $infile_path[-1] ;
				@infile_path				= ()
		}
		
		# if no sequences defined go to help
		else	{ print "\n\tPARAMETER-ERROR: Undefined sequence alignment infile!\n" ; &help( \'i', \$$sref_script, \%$href_par_set_of_opt, \%$href_info_author, \%$href_def_par_setting ) }
		#######################
		
		
		
		#######################
		# Check if clade definition file has been defined
		# Exclude path extension of clade deifinition file and store name in $par_set_of_opt{y}
		if		( $href_par_set_of_opt->{p} ){
		
			my	@infile_path				= split( /\/|\\/, $href_par_set_of_opt->{p} ) ;
				$href_par_set_of_opt->{y}	= $infile_path[-1] ;
				@infile_path				= ()
		}
		else{ $href_par_set_of_opt->{y} = 'undefined' }
		#######################
		
		
		
		#######################
		# Check sequence definition
		# if defined inv site proportion not between 0 and 1:
		if		(	(	$href_par_set_of_opt->{I} < 0 )	||
					(	$href_par_set_of_opt->{I} > 1 )	){
			
			print "\n\tPARAMETER-ERROR: Invariabel site proportion must be between 0 and 1!\n" ; &help( \'I', \$$sref_script, \%$href_par_set_of_opt, \%$href_info_author, \%$href_def_par_setting )
		}
		#######################
		
		
		
		#######################
		# if alpha greater/equal 100 set P4 ML model to non ASRV
		if		( $href_par_set_of_opt->{a} >= 100 ){
		
			$href_par_set_of_opt->{a}	=	100	;
			$href_par_set_of_opt->{I}	=	0	;
		}
		#######################
		
		( @commands, $commandline ) = () ;
	}
	else{ &help( \'h', \$$sref_script, \%$href_par_set_of_opt, \%$href_info_author, \%$href_def_par_setting ) }
	
	## END ARGV handling
	#####################################################################
	
	#####################################################################
	# START PRINT of parameter settings to info file & terminal
	
	#######################
	# define head and botton print out of parameter settings
	my $print_setting_head	= "\n\tPenguin Setup:\n\t---------------------------------------------------\n\t" ;
	my $print_setting_end	= "---------------------------------------------------\n" ;
	#######################
	
	#######################
	# define print out format of single parameter settings
	my %hoh_print_of_parameter = (
									'N'	=>
										{
											0 => "Indel/Amb Sites:\t\trejected (single)\n\t" ,
											1 => "Indel/Amb Sites:\t\trejected (overall)\n\t" ,
										} ,
									'c'	=>
										{
											1 => "Pattern Handling:\t\tclasses\n\t" ,
											0 => "Pattern Handling:\t\tstates\n\t" ,
										} ,
									'r'	=>
										{
											0 => "Additive Root Info:\t\tmixture\n\t" ,
											1 => "Additive Root Info:\t\tspin directive\n\t" ,
										} ,
									'k'	=>
											{
												0 => "Score Weighting:\t\tlowest tree score = 0\n\t" ,
												1 => "Score Weighting:\t\tlowest spin score = 0\n\t" ,
											} ,
									'a'	=>	"\n\tStart Alpha (ML):\t\t".$href_par_set_of_opt->{a} ,
									
									'I'	=>	"\n\tStart pINV (ML):\t\t".$href_par_set_of_opt->{I} ,
									
									'm'	=>	"\n\tSubstitution Model:\t\t".$href_par_set_of_opt->{m} ,
									
									'l'	=>	"\n\n\tMaximum Limit Quartets:\t\t".$href_par_set_of_opt->{l} ,
									
									'M'	=>	"\n\tMinimum Sequence Length:\t".$href_par_set_of_opt->{M} ,
									
									'p'	=>	"\n\n\tClan Definition Infile:\t\t".$href_par_set_of_opt->{y} ,
									
									'i'	=>	"\n\tMSA Quartet Infile:\t\t".$href_par_set_of_opt->{z}."\n\t",
								) ;
	#######################
	
	#######################
	# open info file
	# terminal & file print out of choosen parameter settings
	open OUTinfo, ">>$$sref_inf_out" or die "OUTFILE-ERROR: Cannot open info outfile ", $$sref_inf_out, "in subroutine &argv_handling!\n" ;
	
	print OUTinfo 	$print_setting_head ;
	print 			$print_setting_head ;
	
	for my $option ( qw/N c r k/ ){
		
		print OUTinfo	$hoh_print_of_parameter{$option}{$href_par_set_of_opt->{$option}} ;
		print 			$hoh_print_of_parameter{$option}{$href_par_set_of_opt->{$option}}
	}
	
	for my $option ( qw/m a I l M p i/ ){
		
		print OUTinfo	$hoh_print_of_parameter{$option} ;
		print 			$hoh_print_of_parameter{$option}
	}
	
	
	print OUTinfo	$print_setting_end ;
	print 			$print_setting_end ;
	
	close OUTinfo ;
	#######################
	
	#######################
	# empty print variabels which are unnecessary for further steps
	( %hoh_print_of_parameter, $print_setting_head, $print_setting_end ) = () ;
	#######################

	# END PRINT of parameter settings to info file & terminal
	#####################################################################
}

sub msa_infile_check{
	
	# &msa_infile_check
	# &msa_infile_check ( \%parameter_setting_of_option, \%data_of_infile_property, \%unreduced_seq_of_tax, \%hol_sequence_states_of_taxon) ;
	
	my	$href_setting_of_option				= $_[0]	;	#	key: option command, value: option value												In -> Defined, OUT -> Unchanged
	my	$href_info_of_property				= $_[1]	;	#	key: property; value : property value													In -> Undefined, OUT -> defined
	my	$href_unred_seq_of_taxon			= $_[2]	;	#	key: taxonname; value sequence -> defined via ARGV (-q option)							In -> Undefined, OUT -> defined
	my	$href_hol_sequence_states_of_taxon	= $_[3]	;	#	key1: taxon name; key2: state position number; value: sequence state at that position	In -> Undefined, OUT -> defined
	my	$sref_info_outfile_tmp				= $_[4]	;	#	temporary info outfile name																In -> Defined, OUT -> Unchanged
	my	$sref_scriptname					= $_[5]	;	#	scriptname																				In -> Defined, OUT -> Unchanged
	my	$href_info_author					= $_[6]	;	#	key: author code; value: info															In -> Defined, OUT -> Unchanged
	my	$href_def_par_setting				= $_[7]	;	#	key: parameter code; value: default value												In -> Defined, OUT -> Unchanged
	
	print "\n\tREAD IN MSA Infile ", $href_setting_of_option->{z} ,"..." ;
	
	#######################
	# read in msa infile and store info in %unred_seq_of_tax
	# key: taxon name; value: taxon associated unreduced sequence
	# check of multiple taxon names and
	# non-alphanumeric signs in taxonnames
	if 		( $href_setting_of_option->{z} =~ /\.fas$|\.fasta$/i	){ %$href_unred_seq_of_taxon = &readin_fasta	( \%$href_setting_of_option, \$$sref_scriptname, \%$href_info_author, \%$href_def_par_setting ) }
	elsif 	( $href_setting_of_option->{z} =~ /\.phy$|\.phylip$/i	){ %$href_unred_seq_of_taxon = &readin_phylip	( \%$href_setting_of_option, \$$sref_scriptname, \%$href_info_author, \%$href_def_par_setting ) }
	
	# for my $taxon (sort keys %unred_seq_of_tax){ print "\nunreduced sequence taxon ", $taxon, ":\n", $unred_seq_of_tax{$taxon}, "\n" }
	#######################
	
	##############################################
	# START data check
	my %hol_sequence_states_of_taxon ;
	$href_info_of_property->{type} = 'nuc' ;
	
	#######################
	# check seq lengths -> store msa length to $href_info_of_property->{length}
	# identify data type -> srore data type to $href_info_of_property->{type}
	# check of forbidden sequence signs
	for my $taxon ( sort keys %$href_unred_seq_of_taxon ){
		
		#######################
		# check equal sequence lengths
		unless ( $href_info_of_property->{length} ){ $href_info_of_property->{length} = length $href_unred_seq_of_taxon->{$taxon} }
		else{ unless ( $href_info_of_property->{length} == length $href_unred_seq_of_taxon->{$taxon} ){ print "\n\tINFILE-ERROR: Unequal sequence lengths in msa infile!\n"; &help( \'i', \$$sref_scriptname, \%$href_setting_of_option, \%$href_info_author, \%$href_def_par_setting ) } }
		# print "\nsequence length taxon ", $taxon, ": ", $href_info_of_property->{length}, "\n";
		#######################
		
		#######################
		# substitute lower cased site states to upper cased site states
		# substitute U (Uracil) to T (Thymine)
		$href_unred_seq_of_taxon->{$taxon} =~ s/(\w+)/\U$1/gi ;
		$href_unred_seq_of_taxon->{$taxon} =~ s/U/T/g ;
		# print "\nunreduced sequence after changing lower cases to upper cases and U to T taxon: $taxon\n", $unred_seq_of_tax{$taxon}, "\n";
		#######################
		
		#######################
		# store single site states in identic order as list in hol @hol_sequence_states_of_taxon
		# key: taxon names; values: list elements site states of corresponding taxon
		my @seq_states = split "", $href_unred_seq_of_taxon->{$taxon} ;
		for my $state_pos ( 0 .. @seq_states-1 ){
			
			#######################
			# Check of forbidden sequence signs
			unless ( $seq_states[$state_pos] =~ /A|C|G|T|N|Y|R|W|S|K|M|D|V|H|I|E|L|Q|F|P|B|X|-|\?/ ){
				
				print "\n\tINFILE-ERROR: Forbidden sign (", $seq_states[$state_pos], ") in sequence of ", $taxon, "!\n"; &help( \'i', \$$sref_scriptname, \%$href_setting_of_option, \%$href_info_author, \%$href_def_par_setting )
			}
			#######################
			
			#######################
			# Identify data  type
			# if sequence includes characters only defined for amino acid
			# change data type from nuc to aa
			# afterwards stop checking data type
			unless ( $href_info_of_property->{type} eq 'aa' ){
				
				if ( $seq_states[$state_pos] =~ /E|L|Q|F|P|X/i ){ $href_info_of_property->{type} = 'aa' }
			}
			#######################
			
			$href_hol_sequence_states_of_taxon->{$taxon}[$state_pos] = $seq_states[$state_pos]
		}
		#######################
	}
	#######################
	
	
	
	#######################
	# open info file
	# terminal & file print out of choosen parameter settings
	open	OUTinfo,	">>$$sref_info_outfile_tmp" or die "Outfile-ERROR: Cannot open info outfile ", $$sref_info_outfile_tmp, "in subroutine &msa_infile_check!\n" ;
	
	print 	OUTinfo		"\n\tMSA Sequence Type:\t", $href_info_of_property->{type} ;
	print 				"\n\tMSA Sequence Type:\t", $href_info_of_property->{type} ;
	
	print 	OUTinfo		"\n\tMSA Sequence Length:\t", $href_info_of_property->{length}, "bp\n" ;
	print 				"\n\tMSA Sequence Length:\t", $href_info_of_property->{length}, "bp\n" ;
	
	close	OUTinfo ;
	#######################
	
	# END data check
	##############################################
}

sub readin_fasta{
	
	# &readin_fasta
	# my %unreduced_msa_seq_of_taxon = &readin_fasta ;
	
	my	$href_set_opt		= $_[0]	;	#	key: option command, value: option value	In -> Defined, OUT -> Unchanged
	my	$sref_script		= $_[1]	;	#	scriptname									In -> Defined, OUT -> Unchanged
	my	$href_inf_author	= $_[2]	;	#	key: author code; value: info				In -> Defined, OUT -> Unchanged
	my	$href_def_par_set	= $_[3]	;	#	key: parameter code; value: default value	In -> Defined, OUT -> Unchanged
	
	####################################
	# READ IN fasta interleaved and non-interleaved formated files
	# store single taxa and their associated sequence in has %seq_of_tax
	# key: taxon name ; value: taxon associated sequence
	my 	%seq_of_tax ; my $taxon ; my %seen_fasta_taxon ;
	my $fas_file	= $href_set_opt->{i} ;
	
	open INfas, $fas_file or die "\n\t!INFILE-ERROR!: Cannot Read IN ", $fas_file, "!\n";
	while ( my $line = <INfas> ){
		
		chomp	$line ;
		
		#######################
		# identification of taxon names
		# delete spaces within taxon names
		# check if taxon names appear multiple times
		# check if taxon names consist of only alphanumeric signs
		if ( $line =~ /^\>/ )	{
			
			( $taxon = $line ) =~ s/^\>|\s+// ;
			
			if ( ( $taxon =~ /\w+/ ) && ( $seen_fasta_taxon{$taxon} ) ){
				
				print "\n\t!MSA INFILE-ERROR!: Taxon ", $taxon, " appears multiple times!\n"; &help( \'i', \$$sref_script, \%$href_set_opt, \%$href_inf_author, \%$href_def_par_set )
			}
			
			elsif( ( $taxon =~ /\w+/ ) ){ $seen_fasta_taxon{$taxon}++ }
			
			else { print "\n\t!INFILE-ERROR!: Taxon ", $taxon, " includes non-alphanumeric signs!\n"; &help( \'i', \$$sref_script, \%$href_set_opt, \%$href_inf_author, \%$href_def_par_set ) }
		}
		#######################
		
		#######################
		# store single sequence lines of identified taxon as hash value of taxon
		elsif ( $seq_of_tax{$taxon} ){ $seq_of_tax{$taxon} .= $line }
		else { $seq_of_tax{$taxon} = $line }
		#######################
	}
	
	close INfas ;
	####################################
	
	#######################
	# return taxa and their associated sequence in has %seq_of_tax
	return %seq_of_tax ;
	#######################
}

sub readin_phylip{
	
	# &readin_fasta
	# my %unreduced_msa_seq_of_taxon = &readin_phylip ( \$infilename_with_path ) ;
	
	my	$href_set_opt		= $_[0]	;	#	key: option command, value: option value	In -> Defined, OUT -> Unchanged
	my	$sref_script		= $_[1]	;	#	scriptname									In -> Defined, OUT -> Unchanged
	my	$href_inf_author	= $_[2]	;	#	key: author code; value: info				In -> Defined, OUT -> Unchanged
	my	$href_def_par_set	= $_[3]	;	#	key: parameter code; value: default value	In -> Defined, OUT -> Unchanged
	
	my (
		%seq_of_tax,
		%seen_phylip_taxon
	) ;
	
	
	########################################################################
	## READ IN phylip interleaved and non-interleaved formated files
	## store single taxa and their associated sequence in has $href_seq_of_tax
	my $phy_file	= $href_set_opt->{i} ;
	open INphy , $phy_file or die "\n\t!INFILE-ERROR!: Cannot find ", $phy_file, "!\n" ;
	chomp (my @all_lines_phy = <INphy>) and close INphy ;
	
	####################################
	## exclude empty lines (indelible simulated files often include empty lines)
	my @cleaned_lines_phy ;my $linenumber = 1 ;
	for my $line (@all_lines_phy){ if ( $line =~ /\w+|\?|\-/ ){ push @cleaned_lines_phy, $line; $linenumber++} }
	@all_lines_phy = @cleaned_lines_phy ;
	@cleaned_lines_phy = () ;
	####################################
	
	####################################
	## extract the first line (infoline) and determine the number of taxa ($info_line[0])
	( my $infoline	= shift @all_lines_phy ) =~ s/\s+/ / ;
	my	@infos_line	=  split " ", $infoline ;
	unless ( ( $infos_line[0] =~ /\d+/ ) && ( $infos_line[1] =~ /\d+/ ) ){
	
		print "\n\tINFILE-ERROR: Infile not in phylip format!\n\tMissing first line with taxon and state number!\n";
		&help ( \'i', \$$sref_script, \%$href_set_opt, \%$href_inf_author, \%$href_def_par_set )
	}
	####################################
	
	
	
	####################################
	## phylip files can be in interleaved and non-interleaved format
	## interleaved:	tax1 ACGT...	# first part of all taxa
	## 				tax2 ACGT...
	##								#''space line''
	##				tax1 CCCC...	# Second part of all taxa
	##				tax2 GGGG...
	## to concatenate sequence parts correctly store single lines in separate hashkeys
	## until number of taxa is reached ($c equal $infos_line[0]). afterwards remove the following spaceline and concatenate
	## next lines to their corresponding taxon sequences inferred from the first round and so on...
	## If phylip file is in non-interleaved format, the while lopp stops automatically after the first foreach loop
	my 	%seq_phy = () ;
	while ( @all_lines_phy ){
		
		my $N_lines = @all_lines_phy ; #print "\nremain: ", $N_lines, "\n";
		for ( my $c=1; $c<=$infos_line[0]; $c++ ){ my $seq_line_phy = shift @all_lines_phy ; push ( @{$seq_phy{$c}} , $seq_line_phy ) }
		
		#shift @all_lines_phy ;
	}
	####################################
	#exit;
	
	
	####################################
	## join single sequence parts of each taxon (if interleaved there are multiple key values)
	## taxonnames are in the same line as sequenceinformation (separated by one or multiple whitespaces), therefore
	## substitute multiple whitespaces into one whitespace, join all sequence parts to one string (only important if interleaved format)
	## split taxonnames from sequenceinformation and store both in the hashreference %href_seq_of_tax (key: taxon; value: complete sequence)
	for my $line_c ( sort {$a<=>$b} keys %seq_phy ){ #print "\n", $line_c;
		
		my	@seq_single_parts				=	exists($seq_phy{$line_c}) ? @{$seq_phy{$line_c}} :( ) ;
			$seq_single_parts[0]			=~	s/\s+/:::/ ;
		my	$seq_complete					=	join	""		, @seq_single_parts ;
			$seq_complete					=~	s/\s+//g ;
			@seq_single_parts				=	split	":::"		, $seq_complete ;
		my	$taxon							=	shift @seq_single_parts ; #print "\ntaxon: ", $taxon, "\t", length $seq_single_parts[0], "\n" ;
		
		if ( ( $taxon =~ /\w+/ ) && ( $seen_phylip_taxon{$taxon} ) ){
		
			print "\n\t!INFILE-ERROR!: Taxon ", $taxon, " appears multiple times in msa infile!\n"; &help ( \'i', \$$sref_script, \%$href_set_opt, \%$href_inf_author, \%$href_def_par_set )
		}
		
		elsif ( $taxon =~ /\w+/ ){ $seen_phylip_taxon{$taxon}++ }
		
		else { print "\n\t!INFILE-ERROR!: Taxon ", $taxon, " includes non-alphanumeric signs!\n"; &help ( \'i', \$$sref_script, \%$href_set_opt, \%$href_inf_author, \%$href_def_par_set ) }
			
			$seq_of_tax{$taxon}			=	$seq_single_parts[0] ;
			@seq_single_parts			= () ;
	}
	%seq_phy = () ; #exit;
	########################################################################
	
	#######################
	# return taxa and their associated sequence in has %seq_of_tax
	return %seq_of_tax ;
	#######################
}

sub clade_check{
	
	################################
	# READIN of hashvariabel with keys: sequence names; values: sequence name assigned sequence;;; Name of the clade definition infile
	# RETURN of hashlist with keys: subclade name; values: list of subclade assigned sequence names which have been found in given alignment info of $href_tax_of_msa
	################################
	
	my	$href_tax_of_msa			= $_[0]	;	#	key: taxonname; value sequence				In -> Defined, OUT -> Unchanged
	my	$href_setting_of_option		= $_[1]	;	#	key: option command, value: option value	In -> Defined, OUT -> Unchanged
	my	$sref_info_outfile_tmp		= $_[2]	;	#	temporary info outfile name					In -> Defined, OUT -> Unchanged
	my	$sref_scriptname			= $_[3]	;	#	scriptname									In -> Defined, OUT -> Unchanged
	my	$href_inf_author			= $_[4]	;	#	key: author code; value: info				In -> Defined, OUT -> Unchanged
	my	$href_def_par_setting		= $_[5]	;	#	key: parameter code; value: default value	In -> Defined, OUT -> Unchanged
	
	my (
			%hol_taxa_of_subclade,		# key: cladename; value: list of clade assigned taxa
			%seen_subclade_name,		# key: cladename; value: N occurence
			%seen_subclade_sequence		# key: taxon name; value: N occurence
	) ;
	
	
	
	#######################
	# open info file
	open	OUTinfo,	">>$$sref_info_outfile_tmp" or die "OUTFILE-ERROR: Cannot open info outfile ", $$sref_info_outfile_tmp, "in subroutine &clade_check!\n" ;
	#######################
	
	
	
	####################################
	# if clade definition file defined, assign each subclade
	if ( $href_setting_of_option->{p} ){
		
		####################################
		# READ IN txt file with defined clade subgroups
		# Allowed Format:
		#
		# Subcladecode1,seqname1,seqname2,seqname3,seqname4...\n
		# Subcladecode2,seqname5,seqname6,seqname7,seqname8...\n
		# Subcladecode3,seqname9,seqname10,seqname11,seqname12...\n
		# Subcladecode4,seqnameX,seqnameY,seqnameW,seqnameZ...\n
		#
		# Defined seqnames must be unique in and between each subclade!!
		print			"\n\tREAD IN Clan File: ", $href_setting_of_option->{y} ,"..." ;
		print 	OUTinfo	"\n\tREAD IN Clan File: ", $href_setting_of_option->{y} ,"..." ;
		
		
		my $clade_infile	=	$href_setting_of_option->{p} ;
		open IN,	"<$clade_infile" or die "\n\tINFILE-ERROR: Cannot Read IN ", $clade_infile, "!\n" ;
		
		while ( my $line = <IN> ){
			
			chomp $line ;
			my @lineparts	= split ",", $line ;
			my $cladename	= shift @lineparts ;
			
			
			for ( 0 .. @lineparts-1 ){
				
				#############
				# Check if defined subclade sequence names are given in msa infile, otherwise skip these seqs from further analyses
				if ( $href_tax_of_msa->{$lineparts[$_]} ){ push @{$hol_taxa_of_subclade{$cladename}}, $lineparts[$_] } else { print "\n\tINFILE-WARNING: Cannot find defined sequence name ", $lineparts[$_], " in msa infile!\n" }
				#############
				
				#############
				# Check that each defined subclade sequence name appears only once in given clade infile
				if ( $seen_subclade_sequence{$lineparts[$_]} ){ print "\n\tINFILE-ERROR: Defined sequence name ", $lineparts[$_], " appears multiple times in clan infile!\n"; &help( \'p', \$$sref_scriptname, \%$href_setting_of_option, \%$href_inf_author ,\%$href_def_par_setting ) } else { $seen_subclade_sequence{$lineparts[$_]}++ }
				#############
			}
			
			#############
			# Check that each defined subclade name appears only once in given clade infile
			if ( $seen_subclade_name{$cladename} ){ print "\n\tINFILE-ERROR: Defined clan ", $cladename, " appears multiple times in clan infile!\n"; &help( \'p', \$$sref_scriptname, \%$href_setting_of_option, \%$href_inf_author,\%$href_def_par_setting ) } else { $seen_subclade_name{$cladename}++ }
			#############
		}
		
		close IN ;
	}
	####################################
	
	
	
	####################################
	# if no clade definition file defined, take each taxon as single subclade
	else{
		
		print			"\n\tGenerate Clans from MSA Infile...\n" ;
		print	OUTinfo	"\n\tGenerate Clans from MSA Infile...\n" ;
		
		for my $taxon ( sort keys %$href_tax_of_msa ){
			
			
			unless ( $seen_subclade_name{$taxon} ){
				
				push @{$hol_taxa_of_subclade{$taxon}}, $taxon;
				$seen_subclade_name{$taxon}++
			}
			
			else { print "\n\tINFILE-ERROR: Defined taxon ", $taxon, " appears multiple times in msa infile!\n"; &help( \'p', \$$sref_scriptname, \%$href_setting_of_option, \%$href_inf_author,\%$href_def_par_setting ) }
		}
	}
	####################################
	
	
	
	#######################
	# terminal & file print out of choosen parameter settings
	my $clade_counter = 0 ;
	print 	OUTinfo		"\n\n\tDefined Clans:" ; for ( sort keys %seen_subclade_name ){ $clade_counter++; print OUTinfo	"\n\tClan ", $clade_counter, ":\t", $_ } $clade_counter = 0 ;
	print 				"\n\n\tDefined Clans:" ; for ( sort keys %seen_subclade_name ){ $clade_counter++; print 			"\n\tClan ", $clade_counter, ":\t", $_ }
	
	print 	OUTinfo		"\n" ;
	print 				"\n" ;
	
	close	OUTinfo ;
	#######################
	
	
	
	####################################
	# Check if number of defined subclades equal 4
	my $N_subclades = keys %seen_subclade_name ; print "\n\tN Clans: ", $N_subclades ;
	if ( $N_subclades != 4 ){ print "\n\tINFILE-ERROR: Defined number of taxa clans must be equal 4!\n"; &help( \'p', \$$sref_scriptname, \%$href_setting_of_option, \%$href_inf_author,\%$href_def_par_setting ) }
	####################################
	
	
	
	#######################
	# return subclades and their associated sequence as hash %hol_taxa_of_subclade
	return %hol_taxa_of_subclade ;
	#######################
}

sub site_exclusion{
	
	# &site_exclusion
	# &site_exclusion( \%hol_sequence_states_of_taxon, \%unreduced_sequence_of_taxa, \%data_of_infile_property ) ;
	
	my $href_unred_seq_of_taxon				= $_[0] ;	# key: taxonname; value sequence														-> IN (not changed)
	my $href_info_of_property				= $_[1] ;	# key: property; value : property value													-> IN/OUT (changed)
	my $href_setting_of_option				= $_[2] ;	# key: option command, value: option value												-> IN (not changed)
	my $href_rejected_site_positions		= $_[3] ;	# key: site position numbers which are rejected in further split analysis				-> OUT (defined)
	my $href_allowed_states					= $_[4] ;	# key: allowed character states; value: 1												-> IN (not changed)
	my $href_hol_sequence_states_of_taxon	= $_[5] ;	# key1: taxon name; key2: state position number; value: sequence state at that position	-> IN (not changed)
	
	
	#######################
	# store sorted taxon names as list
	my @taxa		= sort keys %$href_unred_seq_of_taxon ;
	#######################
	
	
	
	#####################################################################
	# START exclusion of uninformative sites, pattern recoding, site character frequency recoding
	# identify uninformative sites:
	# -ambiguities: nuc -> /N|Y|R|W|S|K|M|D|V|H|B/; aa -> /X/
	# -gaps (optionally)                               -> reject: $href_setting_of_option->{N} or $href_setting_of_option->{n} == 1
	# store remaining site positions in @aref_remaining_site_positions
	
	
	#######################
	# print terminal info
	my		$string_taxa	= join ":", @taxa ;
	
	if		( $href_setting_of_option->{N} == 0	){ print "\n\n\tQuartet Indel & Ambig Exclusion" }
	elsif	( $href_setting_of_option->{N} == 1	){ print "\n\n\tMSA Indel & Ambig Exclusion" }
	else	{	die "\n\tBUG-ERROR: Cannot assign exclusion parameter n & A in subroutine &site_exclusion\n\tPlease, report BUG to system developer!\n"		}
	
	$string_taxa	= (  ) ;
	#######################
	
	
	#######################
	# analyse quartet states of each sequence positions
	# if site positions identified as rejected skip directly to loop begin START:
	START:
	for my $seq_pos ( 0 .. $href_info_of_property->{length}-1 ){
		
		##############################################
		# START exclude site positions of undiserable ambiguity and indel states
		for my $taxon ( @taxa ){
			
			unless ( $href_allowed_states->{$href_hol_sequence_states_of_taxon->{$taxon}[$seq_pos]} ){
				
				#######################
				# -exclude sites of ambiguity states
				# -exclude indel site positions if optionally set to reject -> $href_setting_of_option->{N} == 0
				$href_rejected_site_positions->{$seq_pos}++; next START
			}
			#######################
		}
		# END exclude site positions of undiserable ambiguity and indel states
		##############################################
	}
	
	@taxa = () ;
	#####################################################################
}

sub generate_quartets{
	
	################################
	# READIN hashlist; key values: subclade names; values: list of subclade defined taxon names
	# RETURN of listvariabel containing all possible quartet combinations between taxa of predefined subclades
	################################
	
	my $href_hol_taxa_of_subclades	= $_[0] ; # key: subclade name; value: list of subclade assigned taxa		-> IN (not changed)
	
	my @subclades = sort keys %$href_hol_taxa_of_subclades ;
	
	## test print
	#for ( @subclades ){ print $_, "\n" }
	####
	
	my @generated_quartets ;
	
	until ( @subclades < 4){
		
		my $clade1		= shift @subclades ;
		my @cladetaxa1	= @{$href_hol_taxa_of_subclades->{$clade1}} ;
		
		## test print
		#print "\n\tsubclade ", $clade1, ":" ;
		#for ( @cladetaxa1 ){ print "\n\t$_" } print "\n" ;
		####
		
		my $k = 0 ;
		unless ( $k > @subclades - 3 ){
			
			my @cladetaxa2	= @{$href_hol_taxa_of_subclades->{$subclades[$k]}} ;
			my @cladetaxa3	= @{$href_hol_taxa_of_subclades->{$subclades[$k+1]}} ;
			my @cladetaxa4	= @{$href_hol_taxa_of_subclades->{$subclades[$k+2]}} ;
			
			for my $taxon1 ( @cladetaxa1 ){
				
				for my $taxon2 ( @cladetaxa2 ){
					
					for my $taxon3 ( @cladetaxa3 ){
						
						for my $taxon4 ( @cladetaxa4 ){
							
							my @sampled = sort ($taxon1,$taxon2,$taxon3,$taxon4);
							push @generated_quartets, $sampled[0].":".$sampled[1].":".$sampled[2].":".$sampled[3] ;
						}
					}
				}
			}
			
			$k++
		}
	}
	
	return @generated_quartets
}

sub perform_quartet_analyses{

	my	$aref_quartets						= $_[0]		;	# list of generated quartets (each equal $taxon1.":".$taxon2.":".$taxon3.":".$taxon4)		In -> Defined, OUT -> Unchanged
	my	$href_seq_of_tax					= $_[1]		;	# key: taxonname; value sequence															In -> Defined, OUT -> Unchanged
	my	$href_data_of_infile_property		= $_[2]		;	# key: property; value : property value														In -> Defined, OUT -> Unchanged
	my	$href_parameter_setting_of_option	= $_[3]		;	# key: option command, value: option value													In -> Defined, OUT -> Unchanged
	my	$href_rejected_site_positions		= $_[4]		;	# key: site position numbers which are rejected in further split analysis					In -> Defined, OUT -> Changed
	my	$href_hol_taxa_of_subgroup			= $_[5]		;	# key: subclade name; value: list of subclade assigned taxa									In -> Defined, OUT -> Unchanged
	my	$href_hol_sequence_states_of_taxon	= $_[6]		;	# key1: taxon name; key2: state position number; value: sequence state at that position		In -> Defined, OUT -> Unchanged
	my	$href_grafik_colours				= $_[7]		;	# key: topology support (0,1, or 2); value: assigned svg color								In -> Defined, OUT -> Unchanged
	my	$href_graphic_title					= $_[8]		;	# key: graphic outfilename; value: graphic title											In -> Defined, OUT -> Unchanged
	my	$href_txt_filenames					= $_[9]		;	# key: output code of txt files; value: output filename										In -> Defined, OUT -> Unchanged
	my	$href_graphic_setup					= $_[10]	;	# key: svgoption; value: option setup														In -> Defined, OUT -> Unchanged
	my	$href_name_subfolder				= $_[11]	;	# key: subfolder code; value: subfolder name												In -> Defined, OUT -> Unchanged
	
	
	my (
			%hoh_dist_1st_cla1_cla2				,	# Weighting scores substituted to distance scores; key1: sister-taxon1; key2: sister-taxon2; value: distance value (1/split score)
			%hoh_dist_1st_2nd_cla1_cla2 		,	# Weighting scores substituted to distance scores; key1: sister-taxon1; key2: sister-taxon2; value: distance value (1/split score)
			%hoh_seen_quartet					,	# key1: support number (0,1,2 -> 0 element: best supported quartet ), key2: topology; value: count number
			%hol_values_of_quartet				,	# hashlist key: subclade topology, value: list of split scores observed from corresponding quartets
			%hol_taxvalues_of_subclade			,	# key: best subclade topology ; value: taxon quartet => best split value ;
			%rejected_quartet_of_number			,	# key: quartet analysis number value: taxon quartet string
			%seen_quartets_all					,	# key: quartet topology; value: number of occurence
			%hol_val_of_top_of_tax				,	# key1: clade topology, key2: seq name; value: list of split score obtained for given taxa - topology combination
			$remaining_pos_all					,	# sum of overall remaining site position of quartets whose remaining sequence length is above threshold M
			$rejected_pos_all					,	# sum of overall rejected site position of quartets whose remaining sequence length is above threshold M
			%seen_taxa							,	# key: sequence name; value: N occurence in quartet analyses
			%sum_value_of_spin_dir				,	# key: spin directive clan topology; value: total score after single quartet analyses
			%spin_value_of_quartet_topology 	,	# key: quartet_topology; value: list of analysed spin_values of each quartet
			$N_startquartets					,	# Number of single quartet analyses
			%hoh_counter_order_top				,	# key1: quartet tree; key2: order position; value: N occurence
	) ;
	
	
	#######################
	# Generate Subfolder for print OUT
	for ( keys %$href_name_subfolder ){ mkdir $name_of_subfolder{$_} }
	#######################
	
	
	#######################
	# Open info file for info print out
	if ( ( $href_parameter_setting_of_option->{N} == 0 ) && ( $href_txt_filenames->{txtrejected} ne 'none' ) ){
	
		my		$outrej		= $href_txt_filenames->{txtrejected} ;
		open	OUTrej,		">$href_name_subfolder->{txt}/$outrej"	or die "OUTFILE-ERROR: Cannot open info outfile ", $href_name_subfolder->{txt}, "/", $outrej, "in subroutine &perform_quartet_analyses!\n" ;
	}
	
	my		$outtmp		= $href_txt_filenames->{tmpinfo} ;
	open	OUTinfo,	">>$outtmp"	or die "OUTFILE-ERROR: Cannot open info outfile ", $outtmp, "in subroutine &perform_quartet_analyses!\n" ;
	#######################
	
	
	################################
	# Store subclades as hashvalue for each sequence name (key), e.g. key: Libellula value: Hexapoda
	my %subclade_of_taxon ; # key: taxon name; value: assigned taxon clade
	for my $subclade ( sort keys %$href_hol_taxa_of_subgroup ){
	
		my		@cladetaxa	= @{$href_hol_taxa_of_subgroup->{$subclade}} ;
		for (	@cladetaxa ){ $subclade_of_taxon{$_} = $subclade }
	}
	
	$N_startquartets = @$aref_quartets ;
	
	### test print
	#for ( keys %subclade_of_taxon ){ print $_, "\t", $subclade_of_taxon{$_}, "\n" }
	################################
	
	
	#################################################################################################################################################### START Single Quartet Analyses
	LOOP:
	for my $quartet_number ( 1 .. @$aref_quartets ){
		
		#######################
		# store quartet taxa as list variabel (sorted by taxon name)
		# terminal print
		my $quartet_taxa	= @$aref_quartets[$quartet_number-1] ;
		my @q_taxa				= split ":", $quartet_taxa ;
		
		if ( ( $href_parameter_setting_of_option->{N} == 0 ) && ( $href_txt_filenames->{txtrejected} ne 'none' ) ){
		
			print	OUTrej		"Split Analysis q", $quartet_number, ":\t", $quartet_taxa ;
		}
		print					"\n\tStart Split Analysis q", $quartet_number, " (of ", $N_startquartets, "): ", $quartet_taxa ;
		#######################
		
		
		
		#######################
		# store sorted taxa of original msa as list in @sorted_taxa
		# and define split associated topolgies
		# as well as singleton split associated taxa
		my $split_topology1	= '(('.$q_taxa[0].','.$q_taxa[1].'),('.$q_taxa[2].','.$q_taxa[3].'))' ;	# newick topology 1
		my $split_topology2	= '(('.$q_taxa[0].','.$q_taxa[2].'),('.$q_taxa[1].','.$q_taxa[3].'))' ;	# newick topology 2
		my $split_topology3	= '(('.$q_taxa[0].','.$q_taxa[3].'),('.$q_taxa[1].','.$q_taxa[2].'))' ;	# newick topology 3
		my @split_topologies	= ( $split_topology1, $split_topology2, $split_topology3  ) ;
		#######################
		
		
		
		#######################
		# Assign quartet sequences as new hashvalue
		my %quartet_seq_of_tax ;	# key: quartet taxon; value: sequence of taxon
		for my $taxon ( @q_taxa ){ $quartet_seq_of_tax{$taxon} = $href_seq_of_tax->{$taxon} }
		#######################
		
		
		
		############################################################################################ START Site Exclusion of single quartet sequence data
		my %rejected_site_quartet_pos	= %$href_rejected_site_positions ;
		if ( $href_parameter_setting_of_option->{N} == 0 ){
					
					
					#######################
					# Define allowed sequence states
					# all other sequence states are rejected from further split analyses
					my %allowed_states ; # include all allowed site characters
					
					if			( $href_data_of_infile_property->{type} =~ 'nuc'	){ for ( qw/A C G T/																	){ $allowed_states{$_}++ } } # all nuc states allowed
					elsif		( $href_data_of_infile_property->{type} =~ 'aa'		){ for ( qw/A C G T N Y R W S K M D V H I E L Q F P/	){ $allowed_states{$_}++ } } # all aa states allowed
					else { die "\n\tBUG-ERROR: Cannot assign data type in subroutine &perform_quartet_analyses!\n\tPlease, report BUG to system developer!\n" }
					#######################
					
					
					
					#######################
					# Exlusion of indel site positions over all alignment sequences if...
					# ...-n = 1 option or -A = 1
					
					&site_exclusion (
										
										\%quartet_seq_of_tax 					,	# key: quartet taxon; value: sequence of taxon												In -> Defined, OUT -> Unchanged
										\%$href_data_of_infile_property			,	# key: property; value : property value														In -> Defined, OUT -> Changed
										\%$href_parameter_setting_of_option		,	# key: option command, value: option value													In -> Defined, OUT -> Unchanged
										\%rejected_site_quartet_pos				,	# key: site position numbers which are rejected in further split analysis					In -> Defined, OUT -> Changed
										\%allowed_states						,	# key: allowed character states; value: 1													In -> Defined, OUT -> Unchanged
										\%$href_hol_sequence_states_of_taxon	,	# key1: taxon name; key2: state position number; value: sequence state at that position		In -> Defined, OUT -> Unchanged
					) ;
					#######################
					
					
					
					#######################
					# open info file
					# terminal & file print out of choosen parameter settings
					my $N_rej_site_pos_overall	=	keys %$href_rejected_site_positions ;
					my $N_rej_site_pos_total		=	keys %rejected_site_quartet_pos ;
					my $N_rejected_pos_quartet	=	0 + $N_rej_site_pos_total - $N_rej_site_pos_overall ;
					my $remaining_pos_total			=	$href_data_of_infile_property->{length} - $N_rej_site_pos_total ;
					
					
					
					if ( $href_txt_filenames->{txtrejected} ne 'none' ){ print OUTrej 	"\tRejected Quartet Positions:\t",		$N_rejected_pos_quartet , "\tRemaining Quartet Positions:\t",	$remaining_pos_total }
					print 																"\n\tRejected Quartet Positions:\t",	$N_rejected_pos_quartet , "\n\tRemaining Quartet Positions:\t",	$remaining_pos_total ;
					#######################
		}
		#######################
		
		
		
		##############################################
		# Identification of remaining site positions
		my @remaining_site_positions ;
		
		for ( 0 .. $href_data_of_infile_property->{length}-1 ){ unless ( $rejected_site_quartet_pos{$_} ){ push @remaining_site_positions, $_ } }
		
		my $remaining_pos = @remaining_site_positions ;
		if ( $remaining_pos < $href_parameter_setting_of_option->{M} ){
			
			print		"\n\n\tCannot Analyse Quartet: ", $quartet_taxa ,"\n\tRemaining Sequence Length:\t", $remaining_pos, "<", $href_parameter_setting_of_option->{M},
							"bp\n\t----\n\n" ;
			
			if ( ( $href_parameter_setting_of_option->{N} == 0 ) && ( $href_txt_filenames->{txtrejected} ne 'none' ) ){
				
				print OUTrej	"\trejected!\n"
			}
			
			$rejected_pos_all	+= $remaining_pos ;
			$rejected_quartet_of_number{$quartet_number}	=	$quartet_taxa ; next LOOP
		}
		#else { $remaining_pos_all	+= $remaining_pos; if ( ( $href_parameter_setting_of_option->{N} == 0 ) && ( $href_txt_filenames->{txtrejected} ne 'none' ) ){ print OUTrej	"\n"} }
		############################################################################################ END Site Exclusion of single quartet sequence data
		
		
		
		##############################################
		# generate reduced sequence
		my	%reduced_seq_of_tax ;	# key: quartet taxon; value: reduced sequence
		
		for my $seq_pos ( @remaining_site_positions ){
			
			##############################################
			# Assign site characters to @pattern for recoding
			my @pattern ;
			
			$reduced_seq_of_tax{$q_taxa[0]}	.=	$href_hol_sequence_states_of_taxon->{$q_taxa[0]}[$seq_pos] ;
			$reduced_seq_of_tax{$q_taxa[1]}	.=	$href_hol_sequence_states_of_taxon->{$q_taxa[1]}[$seq_pos] ;
			$reduced_seq_of_tax{$q_taxa[2]}	.=	$href_hol_sequence_states_of_taxon->{$q_taxa[2]}[$seq_pos] ;
			$reduced_seq_of_tax{$q_taxa[3]}	.=	$href_hol_sequence_states_of_taxon->{$q_taxa[3]}[$seq_pos] ;
			##############################################
		}
		##############################################
		
		
		
		#################################################################################################################### START P4 handling
		
		########################################################## START PRINT OUT QUARTET ALIGNMENT
		# Print OUT quartet sequence alignment in phylip format for ML pattern estimation
		my $ml_infile_name ;
		&create_msa(
							\%reduced_seq_of_tax,					# key: quartet taxon; value: reduced sequence		-> IN (not changed)
							\$ml_infile_name						# Reduced quartet alignment filename				-> IN (not changed)
		) ;
		########################################################## END PRINT OUT QUARTET ALIGNMENT
		
		
		
		########################################################## START PRINT OF P4 SCRIPT
		# Print OUT p4 pattern script for ML pattern estimation
		my $p4_infile_name ;
		&write_p4_script(
							\$ml_infile_name,							# p4 alignment infile				In -> Defined, OUT -> Unchanged
							\$split_topology1,							# quartet topology (tax1,tax2)		In -> Defined, OUT -> Unchanged
							\$split_topology2,							# quartet topology (tax1,tax3)		In -> Defined, OUT -> Unchanged
							\$split_topology3,							# quartet topology (tax1,tax4)		In -> Defined, OUT -> Unchanged
							\$href_parameter_setting_of_option->{a},	# gamma start value					In -> Defined, OUT -> Unchanged
							\$href_parameter_setting_of_option->{I},	# pinv start value					In -> Defined, OUT -> Unchanged
							\$href_parameter_setting_of_option->{m},	# subst. model						In -> Defined, OUT -> Unchanged
							\$p4_infile_name,							# name of the p4 infile				In -> Undefined, OUT -> defined
							\$href_data_of_infile_property->{type}		# identified seq. type (nuc or aa)	In -> Defined, OUT -> Unchanged
		) ;
		########################################################## END PRINT OF P4 SCRIPT
		
		#exit;
		
		########################################################## START P4 PROGRAM
		# Execute ML pattern estimation using subscript p4
		print "\n\n\tStart P4..." ;
		my $program			= "p4";
		my $p4_outfile	= "p4_result.txt" ;
		system ( "$program $p4_infile_name >$p4_outfile" );
		unlink ( $p4_infile_name, $ml_infile_name ) ;
		########################################################## END P4 PROGRAM
		
		#exit;
		
		########################################################## START EXTRACTION OF P4 RESULTS
		# START result extractions of ML pattern evaluation program -> p4
		my (
				%hoh_found_N_of_topo_of_pattern						, # key1: topology ( e.g. ((S1,L5),(S2,L6)) ); key2: recoded site pattern ( e.g. A ); value: N site pattern in original data ( e.g. 200 )
				%hoh_expected_N_of_topo_of_pattern					, # key1: topology ( e.g. ((S1,L5),(S2,L6)) ); key2: recoded site pattern ( e.g. A ); value: N expected site pattern via ML branch length estimation of given topology ( e.g. 200 )
				$p4_error_message									, # defined by 1 if p4 stops with an error prompt
		) ;
		&read_in_p4_results(
		
							\%hoh_found_N_of_topo_of_pattern		,	# key1: topology; key2: recoded site pattern; value: N site pattern in original data	In -> Undefined, OUT -> defined
							\%hoh_expected_N_of_topo_of_pattern		,	# key1: topology; key2: recoded site pattern; value: expected ML N site pattern			In -> Undefined, OUT -> defined
							\$p4_outfile							,	# name of the p4 outfile																In -> Defined, OUT -> Unchanged
							\%$href_parameter_setting_of_option		,	# key: option command, value: option value												In -> Defined, OUT -> Unchanged
							\%$href_data_of_infile_property			,	# key: property; value : property value													In -> Defined, OUT -> Unchanged
							\$p4_error_message						,	# defined by 1 if p4 stops with an error prompt											In -> Undefined, OUT -> changed
		) ;
		
		# test print
		#print "\ntest:\n",			$split_topology_1,
		#					"\n",	$hoh_found_N_of_topo_of_pattern{$split_topology_1}{K},
		#					"\n",	$hoh_found_N_of_topo_of_pattern{$split_topology_1}{L},
		#					"\n",	$hoh_found_N_of_topo_of_pattern{$split_topology_1}{M},
		#					"\n",	$hoh_found_N_of_topo_of_pattern{$split_topology_1}{N}
		#;exit;
		###
		
		if ( $p4_error_message ){
			
			print "\n\tP4-INPUT-ERROR: ", $p4_error_message ,"\n\tQuartet combination q", $quartet_number ," rejected due to inappropriate P4 model assumptions!\n\n" ;
			if ( ( $href_parameter_setting_of_option->{N} == 0 ) && ( $href_txt_filenames->{txtrejected} ne 'none' ) ){ print OUTrej	"\trejected due to inappropriate P4 model assumptions!\n"}
			
			$rejected_pos_all	+= $remaining_pos ;
			$rejected_quartet_of_number{$quartet_number}	=	$quartet_taxa ; next LOOP
		}
		else { $remaining_pos_all	+= $remaining_pos; if ( ( $href_parameter_setting_of_option->{N} == 0 ) && ( $href_txt_filenames->{txtrejected} ne 'none' ) ){ print OUTrej	"\n"} }
		########################################################## END EXTRACTION OF P4 RESULTS
		
		#################################################################################################################### END P4 handling
		
		#exit;
		
		#################################################################################################################### START SPLIT ANALYSIS
		
		########################################################## T22 -> &determine_quartet_topology_11
		# determine quartet topology by estimating overall pattern distances for each quartet between expected pattern distribution (ML using subscript best_pattern_match_tree) and observed pattern distribution
		# T22: Nap = N observed split signal - N observed signal plesiomorph (Np) - N expected convergent signal (Nk) -> each value reduced by proportion (shortest number of singelton * 4) / total number of singeltons (separately for observed and expected pattern)
		# Convergent split signal reduced of given topology reduced for each of the two others AFTER calculating the mean convergent signal value
		print "\n\tSplit Pattern Evaluation..." ;
		
		#############################
		# Determine Split Topology by comparing ML estimated pattern distributions and idetified pattern distributions
		
		&start_phyquart(
									
									\$quartet_number							,	# quartet number																		In -> Defined, OUT -> Unchanged
									\@split_topologies							,	# list of split topologies (1,2,3)														In -> Defined, OUT -> Unchanged
									\%hoh_found_N_of_topo_of_pattern			,	# key1: topology; key2: recoded site pattern; value: N site pattern in original data	In -> Defined, OUT -> Unchanged
									\%hoh_expected_N_of_topo_of_pattern			,	# key1: topology; key2: recoded site pattern; value: expected ML N site pattern			In -> Defined, OUT -> Unchanged
									\%spin_value_of_quartet_topology			,	# key: quartet_topology; value: list of spin_values										In -> Undefined, OUT -> Defined
									\$href_parameter_setting_of_option->{f}		,	# print code (1 -> print all) (0 - print only important output)							In -> Defined, OUT -> Unchanged
									\$quartet_taxa								,	# string of quartet taxa for calculation print out ('t1:t2:t3:t4')						In -> Defined, OUT -> Unchanged
									\$href_txt_filenames->{txtaddprint}			,	# name single quartet calculation														In -> Defined, OUT -> Unchanged
									\$href_name_subfolder->{txt}				,	# key: subfolder code; value: subfolder name											In -> Defined, OUT -> Unchanged
									\$href_parameter_setting_of_option->{z}		,	# msa infile name																		In -> Defined, OUT -> Unchanged
									\%subclade_of_taxon							,	# key: taxon name; value: defined subclade name of taxon								In -> Defined, OUT -> Unchanged
									\%sum_value_of_spin_dir						,	# key: spin directive clan topology; value: total score after single quartet analyses	In -> Undefined, OUT -> Defined
									\%hoh_counter_order_top						,	# key1: quartet tree; key2: order position; value: N occurence							In -> Undefined, OUT -> Defined
		) ;
		
		print " done.\n\t----\n"
		#############################
	}
	#################################################################################################################### END SINGLE SPLIT ANALYSIS OF GENERATED QUARTETS
	
	
	
	#################################################################################################################### START Determination of Rejected Quartets
	
	#############################
	# Determine info about overall rejected quartets
	my $N_rejected_quartets	= 0 + keys %rejected_quartet_of_number ;
	my $N_final_quartets		= @$aref_quartets - $N_rejected_quartets ;
	
	##############
	# Data Check. DIE if no quartet can be analysed
	if ($N_final_quartets	== 0 ){
		
		print	OUTinfo		"\n\tDATA-ERROR: No quartet left with sequence length above ", $href_parameter_setting_of_option->{M} ,"bp!\n\n" ;
		print 				"\n\tDATA-ERROR: No quartet left with sequence length above ", $href_parameter_setting_of_option->{M} ,"bp!\n\n" ; exit
	}
	##############
	
	##############
	# Print rejected quartets on the terminal and in the rejected info file
	if ( %rejected_quartet_of_number ){
		
		print "\n\t", $N_rejected_quartets, " Rejected Quartets:"	; for ( sort {$a<=>$b} keys %rejected_quartet_of_number ){ print	"\n\tq"	, $_, ": ", $rejected_quartet_of_number{$_} } print "\n\t"	;
		
		if ( ( $href_parameter_setting_of_option->{N} == 0 ) && ( $href_txt_filenames->{txtrejected} ne 'none' ) ){
			
			print OUTrej	"\n\nRejected Quartets:"; for ( sort {$a<=>$b} keys %rejected_quartet_of_number ){ print OUTrej	"\nq"	, $_, ": ", $rejected_quartet_of_number{$_} }
		}
	}
	##############
	
	##############
	# Calculate mean sequence length of rejected quartets
	unless ( $rejected_pos_all ){ $rejected_pos_all = 'none' } else{ $rejected_pos_all /= $N_rejected_quartets }
	##############
	
	##############
	# Calculate mean sequence length of remaining quartets
	my $mean_remain_pos	=	$remaining_pos_all/$N_final_quartets ;
	##############
	
	##############
	# Info txt print of N remaining quartets, mean seq lengths of rej and remaining quartets
	print	OUTinfo	"\n\tN Start Quartets:\t\t"								, $N_startquartets		,
					"\n\tN Rejected Quartets:\t\t"							, $N_rejected_quartets 	,
					"\n\tN Conducted Quartet-Analyses:\t"					, $N_final_quartets		,
					"\n\tMean Number Site Positions Remaining Quartets:\t"	, $mean_remain_pos		,
					"\n\tMean Number Site Positions Rejected Quartets:\t"	, $rejected_pos_all		;
	#############################
	
	#################################################################################################################### END Determination of Rejected Quartets
	
	
	
	#################################################################################################################### START OVERALL SPLIT SCORE ASSIGNMENT
	
	#############################
	# Analysis of overall single quartet scores under seperate consideration of different spins (-r option)
	if ( $href_parameter_setting_of_option->{r} == 1 ){
		
		&spin_directive_clan_evaluation(
									
									\@$aref_quartets							,	# list of generated quartets (each equal $taxon1.":".$taxon2.":".$taxon3.":".$taxon4)													In -> Defined, OUT -> Unchanged
									\%sum_value_of_spin_dir						,	# key: spin directive clan topology; value: total score after single quartet analyses													In -> Defined, OUT -> Unchanged
									\%spin_value_of_quartet_topology			,	# key: quartet_topology; value: list of spin_values																						In -> Defined, OUT -> Changed
									\%hoh_seen_quartet							,	# key1: support number (0,1,2 -> 0 element: best supported quartet ), key2: topology; value: count number								In -> Undefined, OUT -> defined
									\%hol_taxvalues_of_subclade					,	# key1: subclade topology; key2: support code (0->best, 1-> 2nd best); value1: sorted taxonnames (comma separated); value2: split score	In -> Undefined, OUT -> defined
									\%hol_val_of_top_of_tax						,	# key1: subclade topology; key2: taxonname; value: list of taxon observed split scores													In -> Undefined, OUT -> defined
									\%rejected_quartet_of_number				,	# key1: quartet number; value: counter of rejection (0 or 1)																			In -> Defined, OUT -> Unchanged
									\$href_txt_filenames->{txtssupport}			,	# name summarized split support, former M-ICE symetric outfile																			In -> Defined, OUT -> Unchanged
									\$href_parameter_setting_of_option->{z}		,	# msa infile name																														In -> Defined, OUT -> Unchanged
									\$href_name_subfolder->{txt}				,	# key: subfolder code; value: subfolder name																							In -> Defined, OUT -> Unchanged
									\%seen_taxa									,	# key: sequence name; value: N occurence in quartet analyses																			In -> Undefined, OUT -> defined
									\$href_parameter_setting_of_option->{k}		,	# normalization setup (0-> lowest of all 6 spin scores, 1-> lowest of best spins)														In -> Ddefined, OUT -> Unchanged
		) ;
		
		print	OUTinfo	"\n\n\tSpin Dependend Multiple-Clan Analysis.",
						"\n\tIdentified Character Alteration along Internal Branch: Left to Right" ;
	}
	
	else {
		
		&non_directive_clan_evaluation(
									
									\@$aref_quartets							,	# list of generated quartets (each equal $taxon1.":".$taxon2.":".$taxon3.":".$taxon4)													In -> Defined, OUT -> Unchanged
									\%spin_value_of_quartet_topology			,	# key: quartet_topology; value: list of spin_values																						In -> Defined, OUT -> Changed
									\%hoh_seen_quartet							,	# key1: support number (0,1,2 -> 0 element: best supported quartet ), key2: topology; value: count number								In -> Undefined, OUT -> defined
									\%hol_taxvalues_of_subclade					,	# key1: subclade topology; key2: support code (0->best, 1-> 2nd best); value1: sorted taxonnames (comma separated); value2: split score	In -> Undefined, OUT -> defined
									\%hol_val_of_top_of_tax						,	# key1: subclade topology; key2: taxonname; value: list of taxon observed split scores													In -> Undefined, OUT -> defined
									\%rejected_quartet_of_number				,	# key1: quartet number; value: counter of rejection (0 or 1)																			In -> Defined, OUT -> Unchanged
									\$href_txt_filenames->{txtssupport}			,	# name summarized split support, former M-ICE symetric outfile																			In -> Defined, OUT -> Unchanged
									\$href_parameter_setting_of_option->{z}		,	# msa infile name																														In -> Defined, OUT -> Unchanged
									\$href_name_subfolder->{txt}				,	# key: subfolder code; value: subfolder name																							In -> Defined, OUT -> Unchanged
									\%seen_taxa									,	# key: sequence name; value: N occurence in quartet analyses																			In -> Undefined, OUT -> defined
									\$href_parameter_setting_of_option->{k}		,	# normalization setup (0-> lowest of all 6 spin scores, 1-> lowest of best spins)														In -> Ddefined, OUT -> Unchanged
		) ;
	}
	#############################
	
	#################################################################################################################### END OVERALL SPLIT SCORE ASSIGNMENT
	
	
	
	################################################################################################ START Triangle Quartets
	
	################################
	# definition of examined quartet topologies
	my 			@quartets_all	= sort keys %spin_value_of_quartet_topology ;
	unless ( 	@quartets_all == 3 ){ die "BUG-ERROR: Number of best supported quartet trees unequal 3 in subroutine &perform_quartet_analyses!\n\tPlease, report BUG to system developer!\n" }
	################################
	
	
	
	################################
	# Assign variabels to single quartet split score for print OUT
	for ( 0 .. 2 ){
		
		$href_grafik_colours->{$quartets_all[$_]}	=	$href_grafik_colours->{$_};
		#print "\ngr ", $quartets_all[$_], " color ", $href_grafik_colours->{$quartets_all[$_]}
	}
	
	my	%q_number_of_top ;
			$q_number_of_top{$quartets_all[0]}	= 'Q1' ;
			$q_number_of_top{$quartets_all[1]}	= 'Q2' ;
			$q_number_of_top{$quartets_all[2]}	= 'Q3' ;
	
	my	(
			@lol_data_ter		,	# array data matrix for triangle print out list1: integer, quartet number; list2: 0, 1 or 2 -> topology number; value: split score
			%name_of_repeat		,	# key: integer, quartet number, value: result string for triangle print out
			%sum_score_of_top	,	# key: integer (0, 1, or 2) of topology array; value: Total score for each quartet topology
			@lol_scores_of_top	,	# array data matrix list1: integer (0, 1, or 2) of topology array; value: list of single split scores
		) ;
	
	##############
	# for each quartet...
	my	$outcla		= $href_txt_filenames->{txtqsupport} ;
	
	if ( $href_txt_filenames->{txtqsupport} ne 'none'){
		
		open OUTcla, ">$href_name_subfolder->{txt}/$outcla"	or die "OUTFILE-ERROR: Cannot open info outfile ", $href_name_subfolder->{txt}, "/", $outcla, "in subroutine &perform_quartet_analyses!\n"
	}
	
	my		$repeat		=	0;
	
	for my $q_number ( 1 .. @$aref_quartets ){
		
		unless ( $rejected_quartet_of_number{$q_number} ){
			
			if ( $href_txt_filenames->{txtqsupport} ne 'none'){ print OUTcla	"Split Analysis q" , $q_number, ":\t" }
			
			my @values ;
			
			##############
			# for each possible quartet topology
			for my $topology ( 0 .. 2 ){
				
				$lol_data_ter[$repeat][$topology]		=	$spin_value_of_quartet_topology{$quartets_all[$topology]}[$repeat];	# Assign split score for triangle print out
				push @values							,	$spin_value_of_quartet_topology{$quartets_all[$topology]}[$repeat];	# Assign split score for triangle result string
				$sum_score_of_top{$topology}			+=	$spin_value_of_quartet_topology{$quartets_all[$topology]}[$repeat];	# Assign split score for mean split score calculation
				push @{$lol_scores_of_top[$topology]}	,	$spin_value_of_quartet_topology{$quartets_all[$topology]}[$repeat];	# Assign split score for medain split score calculation
			}
			##############
			
			##############
			# Generate triangle result info string for each quartet run
			($name_of_repeat{$repeat}			=	@$aref_quartets[$q_number-1])	=~	s/:/,/g ;
			$name_of_repeat{$repeat}			=	"(".$name_of_repeat{$repeat}.")\nSupport Q1: ".$values[0]."\nSupport Q2: ".$values[1]."\nSupport Q3: ".$values[2]."\n".$href_grafik_colours->{qdots} ;
			##############
			
			##############
			# Result print of each quartet split support for different subclade relationships
			if ( $href_txt_filenames->{txtqsupport} ne 'none'){
				
				print OUTcla 	$aref_quartets->[$repeat], "\tSplit Support ",
								"\tQ1\t", $quartets_all[0] ,":\t", $values[0],
								"\tQ2\t", $quartets_all[1] ,":\t", $values[1],
								"\tQ3\t", $quartets_all[2] ,":\t", $values[2], "\n"
			}
			
			$repeat++
		}
		else { print OUTcla	"Split Analysis q", $q_number, ":\tQuartet sequence length below allowed threshold -> rejected\n" }
	}
	if ( $href_txt_filenames->{txtqsupport} ne 'none'){ close OUTcla }
	################################
	
	
	
	################################ Mean Calculation
	# Caculate Mean Score of single quartet values for each of the three quartet topologies
	my %mean_of_top ;	# key: toponumber (0, 1, or 2); value: mean split score
	print	OUTinfo	"\n\n\tOverall Split Signal (Mean):";
	for ( 0 .. 2 ){
		
		$mean_of_top{$_}	=	$sum_score_of_top{$_}/$N_final_quartets	;
		
		print OUTinfo	"\n\tQ", $_+1 ," ", $quartets_all[$_], ":\t", $mean_of_top{$_}
	}
	################################
	
	
	
	################################ Median Calculation
	# Caculate Median Score of single quartet values for each of the three quartet topologies
	my %median_of_top ;	# key: toponumber (0, 1, or 2); value: median split score
	print	OUTinfo	"\n\n\tOverall Split Signal (Median):";
	
	my $med_total	;
	for ( 0 .. 2 ){
		
		my	@scores_t			=		@{$lol_scores_of_top[$_]}		;
			$median_of_top{$_}	=		&median( \@scores_t )			;
			$med_total			+=		$median_of_top{$_}				;
			
			@scores_t			=		()								;
	}
	
	for ( 0 .. 2 ){
		
		$median_of_top{$_}	=	$median_of_top{$_}	/ $med_total ;
		print OUTinfo	"\n\tQ", $_+1 ," ", $quartets_all[$_], ":\t", $median_of_top{$_}
	}
	
		@lol_scores_of_top	=	();
	################################
	
	
	################################################################################################ START Triangle Quartets
	
	################################ Print Triangle of Quartets
	# print triangle for all quartets
	unless ( $href_txt_filenames->{triall} eq 'none' ){
		
		&print_svg_triangle_graphic(
								
								\@lol_data_ter						,	# list1: data run (quartet) number; list 2: quartet topology; value: split	 support found for given quartet and topology	In -> Defined, OUT -> Unchanged
								\$href_txt_filenames->{triall}		,	# defined output filename																									In -> Defined, OUT -> Unchanged
								\@quartets_all						,	# list of sorted clade topologies																							In -> Defined, OUT -> Unchanged
								\%name_of_repeat					,	# key: integer, quartet number, value: result string for triangle print out													In -> Defined, OUT -> Unchanged
								\$href_graphic_title->{triall}		,	# key: graphic outfilename; value: graphic title																			In -> Defined, OUT -> Unchanged
								\%mean_of_top						,	# key: toponumber (0, 1, or 2); value: mean split score																		In -> Defined, OUT -> Unchanged
								\%median_of_top						,	# key: toponumber (0, 1, or 2); value: median split score																	In -> Defined, OUT -> Unchanged
								\%$href_graphic_setup				,	# key: svgoption; value: option setup																						In -> Defined, OUT -> Unchanged
								\$href_name_subfolder->{svg}		,	# key: subfolder code; value: subfolder name																				In -> Defined, OUT -> Unchanged
								\'all_quartets'						,	# code for median score print out																							In -> Defined, OUT -> Unchanged
								\'Mean'								,	# code for triangle overall mean or median print out																		In -> Defined, OUT -> Unchanged
								\$href_grafik_colours->{triangle}	,	# triangle background color																									In -> Defined, OUT -> Unchanged
								\$href_grafik_colours->{back}		,	# background color svg																										IN -> Defined, OUT -> Unchanged
								\$href_grafik_colours->{header}		,	# header color																												IN -> Defined, OUT -> Unchanged
		)
	}
		@lol_data_ter	= 	()	;
	################################
	
	################################################################################################ END Triangle Quartets
	
	
	
	################################################################################################ START Triangle TAXA MEAN
	
	################################
	# Color adjustment for different
	# circle print of subclades in triangle output for taxa
	my @colors				=	( $href_grafik_colours->{clan1}, $href_grafik_colours->{clan2}, $href_grafik_colours->{clan3}, $href_grafik_colours->{clan4} ) ;
	my %color_of_subclade	;
	my $cladecounter		=	0 ;
	for my $subclade ( sort keys %$href_hol_taxa_of_subgroup ){
		
		$color_of_subclade{$subclade}	=	$colors[$cladecounter];
		$cladecounter++
	}
	################################
	
	
	
	################################
	# print triangle of single taxa mean values
	my		$outtax		= $href_txt_filenames->{txttsupport} ;
	my		$outmed		= $href_txt_filenames->{txttsupmedi} ;
	
	if ( $href_txt_filenames->{txttsupport} ne 'none'){
		
		open	OUTtax,		">$href_name_subfolder->{txt}/$outtax" or die "OUTFILE-ERROR: Cannot open info outfile ", $href_name_subfolder->{txt}, "/", $outtax , "in subroutine &perform_quartet_analyses!\n" ;
	}
	
	if ( $href_txt_filenames->{txttsupmedi} ne 'none'){
		
		open	OUTmed,		">$href_name_subfolder->{txt}/$outmed" or die "OUTFILE-ERROR: Cannot open info outfile ", $href_name_subfolder->{txt}, "/", $outmed , "in subroutine &perform_quartet_analyses!\n" ;
	}
	
		$repeat					=	0	;
	my	@lol_data_ter_med 		=	()	;
	my	%name_of_repeat_median	=	()	;
			%name_of_repeat		=	()	;
	
	for my $taxon ( sort keys %subclade_of_taxon ){
		
		if ( $seen_taxa{$taxon} ){
			
			################
			# For each of the three possible quartet topologies
			# mean and median score of obtained single split support values
			my %val_of_top ;
			my %med_of_top ;
			
			for my $topology ( 0 .. 2 ){
				
				my	@tax_values			=	@{$hol_val_of_top_of_tax{$quartets_all[$topology]}{$taxon}} ;
					@tax_values			=	sort {$a<=>$b} @tax_values ;
				
				$med_of_top{$topology}	=	&median( \@tax_values ) ;
				$val_of_top{$topology}	=	0 ;
				
				for my $val ( @tax_values ){ $val_of_top{$topology} += $val }
				
				$val_of_top{$topology}	=	$val_of_top{$topology}/@tax_values ;
			}
			################
			
			################
			# For each of the three possible quartet topologies
			# calulation of split proportion compared to the two alternatives (q1 support = q1 / q1+q2+q3)
			
			for my $topology ( 0 .. 2 ){
				
				$lol_data_ter[$repeat][$topology]				=	$val_of_top{$topology} / ( $val_of_top{0} + $val_of_top{1} + $val_of_top{2} ) ;
				$lol_data_ter_med[$repeat][$topology]			=	$med_of_top{$topology} / ( $med_of_top{0} + $med_of_top{1} + $med_of_top{2} ) ;
				
				#print "\n", $lol_data_ter_med[$repeat][$topology]
			}
			#print "\n", ($lol_data_ter_med[$repeat][0] + $lol_data_ter_med[$repeat][1] + $lol_data_ter_med[$repeat][2] );
			
			$name_of_repeat_median{$repeat}	=	"Clan: ".$subclade_of_taxon{$taxon}." Seq: ".$taxon."\nMedian Support Q1: ".$lol_data_ter_med[$repeat][0]."\nMedian Support Q2: ".$lol_data_ter_med[$repeat][1]."\nMedian Support Q3: ".$lol_data_ter_med[$repeat][2]."\n".$color_of_subclade{$subclade_of_taxon{$taxon}} ;
			$name_of_repeat{$repeat}		=	"Clan: ".$subclade_of_taxon{$taxon}." Seq: ".$taxon."\nMean Support Q1: ".$val_of_top{0}."\nMean Support Q2: ".$val_of_top{1}."\nMean Support Q3: ".$val_of_top{2}."\n".$color_of_subclade{$subclade_of_taxon{$taxon}} ;
			
			unless ( $href_txt_filenames->{txttsupport} eq 'none'){
				
				print OUTtax	"Split Analysis s", $repeat+1, ":\t", $taxon ,"\tClan:\t", $subclade_of_taxon{$taxon} , "\tMean Split Support ",
								"Q1\t", $quartets_all[0], ":\t" , $val_of_top{0} , "\tQ2:\t", $quartets_all[1], "\t" , $val_of_top{1} , "\tQ3:\t", $quartets_all[2], "\t" , $val_of_top{2}, "\n" ;
			}
			
			unless ( $href_txt_filenames->{txttsupmedi} eq 'none'){
				
				print OUTmed	"Split Analysis s", $repeat+1, ":\t", $taxon ,"\tClan:\t", $subclade_of_taxon{$taxon} , "\tMedian Split Support ",
								"Q1\t", $quartets_all[0], ":\t" , $lol_data_ter_med[$repeat][0] , "\tQ2:\t", $quartets_all[1], "\t" , $lol_data_ter_med[$repeat][1] , "\tQ3:\t", $quartets_all[2], "\t" , $lol_data_ter_med[$repeat][2], "\n" ;
			}
			$repeat++;
			################
		}
	}
	if ( $href_txt_filenames->{txttsupport} ne 'none'){ close OUTtax }
	if ( $href_txt_filenames->{txttsupmedi} ne 'none'){ close OUTmed }
	
	unless ( $href_txt_filenames->{tritax} eq 'none' ){
		
		&print_svg_triangle_graphic(
								
								\@lol_data_ter						,	# list1: data run (quartet) number; list 2: quartet topology; value: mean split support found for given taxon, quartet and topology		In -> Defined, OUT -> Unchanged
								\$href_txt_filenames->{tritax}		,	# defined output filename																												In -> Defined, OUT -> Unchanged
								\@quartets_all						,	# list of sorted clade topologies																										In -> Defined, OUT -> Unchanged
								\%name_of_repeat					,	# key: integer, quartet number, value: result string for triangle print out																In -> Defined, OUT -> Unchanged
								\$href_graphic_title->{tritax}		,	# key: graphic outfilename; value: graphic title																						In -> Defined, OUT -> Unchanged
								\%mean_of_top						,	# key: toponumber (0, 1, or 2); value: mean split score																					In -> Defined, OUT -> Unchanged
								\%median_of_top						,	# key: toponumber (0, 1, or 2); value: median split score																				In -> Defined, OUT -> Unchanged
								\%$href_graphic_setup				,	# key: svgoption; value: option setup																									In -> Defined, OUT -> Unchanged
								\$href_name_subfolder->{svg}		,	# key: subfolder code; value: subfolder name																							In -> Defined, OUT -> Unchanged
								\'taxa'								,	# code for no median score print out																									In -> Defined, OUT -> Unchanged
								\'Mean'								,	# code for triangle overall mean or median print out																					In -> Defined, OUT -> Unchanged
								\$href_grafik_colours->{triangle}	,	# triangle background color																												In -> Defined, OUT -> Unchanged
								\$href_grafik_colours->{back}		,	# background color svg																													IN -> Defined, OUT -> Unchanged
								\$href_grafik_colours->{header}		,	# header color																															IN -> Defined, OUT -> Unchanged
		)
	}
		@lol_data_ter	= 	()	;
	
	unless ( $href_txt_filenames->{trimed} eq 'none' ){
		
		&print_svg_triangle_graphic(
								
								\@lol_data_ter_med					,	# list1: data run (quartet) number; list 2: quartet topology; value: median split support found for given taxon, quartet and topology	In -> Defined, OUT -> Unchanged
								\$href_txt_filenames->{trimed}		,	# defined output filename																												In -> Defined, OUT -> Unchanged
								\@quartets_all						,	# list of sorted clade topologies																										In -> Defined, OUT -> Unchanged
								\%name_of_repeat_median				,	# key: quartet number (0 to @quartets-1); value: quartet taxa (tax1:tax2:tax3:tax4)														In -> Defined, OUT -> Unchanged
								\$href_graphic_title->{trimed}		,	# key: graphic outfilename; value: graphic title																						In -> Defined, OUT -> Unchanged
								\%median_of_top						,	# key: toponumber (0, 1, or 2); value: median split score																				In -> Defined, OUT -> Unchanged
								\%mean_of_top						,	# key: toponumber (0, 1, or 2); value: mean split score																					In -> Defined, OUT -> Unchanged
								\%$href_graphic_setup				,	# key: svgoption; value: option setup																									In -> Defined, OUT -> Unchanged
								\$href_name_subfolder->{svg}		,	# key: subfolder code; value: subfolder name																							In -> Defined, OUT -> Unchanged
								\'taxa'								,	# code for no median score print out																									In -> Defined, OUT -> Unchanged
								\'Median'							,	# code for triangle overall mean or median print out																					In -> Defined, OUT -> Unchanged
								\$href_grafik_colours->{triangle}	,	# triangle background color																												In -> Defined, OUT -> Unchanged
								\$href_grafik_colours->{back}		,	# background color svg																													IN -> Defined, OUT -> Unchanged
								\$href_grafik_colours->{header}		,	# header color																															IN -> Defined, OUT -> Unchanged
		)
	}
		@lol_data_ter	= 	()	;
	################################
	
	################################################################################################ END Triangle TAXA MEAN
	
	
	
	################################################################################################ START Bar Chart of Topology Order Distribution
	
	################################ all spin topologies
	# print svg bar plot output of each spin subclade topology
	unless ( $href_txt_filenames->{barDs} eq 'none' ){
		
		if ( $href_parameter_setting_of_option->{r} == 1 ){
			
			my	@order_number	=	qw/1 2 3 4 5 6/ ;
			my	@all_six_trees	=	sort {$sum_value_of_spin_dir{$b}	<=> $sum_value_of_spin_dir{$a}}	keys %sum_value_of_spin_dir ;
			
			&print_svg_barchart_graphic(
								
								\%$href_grafik_colours				,	# key: best spin topologies; value: assigned svg color			In -> Defined, OUT -> Unchanged
								\@all_six_trees						,	# list of all clade topologies									In -> Defined, OUT -> Unchanged
								\%hoh_counter_order_top				,	# key1: quartet tree; key2: order position; value: N occurence	In -> Defined, OUT -> Unchanged
								\$N_final_quartets					,	# total number of quartets										In -> Defined, OUT -> Unchanged
								\$href_graphic_title->{barDs}		,	# svg graphic title												In -> Defined, OUT -> Unchanged
								\$href_txt_filenames->{barDs}		,	# svg outfile name												In -> Defined, OUT -> Unchanged
								\@order_number						,	# order numbers (1 best to 6 lowest support)					In -> Defined, OUT -> Unchanged
								\%$href_graphic_setup				,	# key: svgoption; value: option setup							In -> Defined, OUT -> Unchanged
								\$href_name_subfolder->{svg}		,	# key: subfolder code; value: subfolder name					In -> Defined, OUT -> Unchanged
								\$href_grafik_colours->{back}		,	# background color svg											In -> Defined, OUT -> Unchanged
								\$href_grafik_colours->{header}		,	# header color													In -> Defined, OUT -> Unchanged
								\$href_grafik_colours->{barleft}	,	# barplot color of lower supported spin trees					In -> Defined, OUT -> Unchanged
								\$href_grafik_colours->{header}		,	# header color													In -> Defined, OUT -> Unchanged
			);
			
			@all_six_trees	=	()
		}
	}
	################################
	
	
	################################ best or summarized topologies
	# print svg bar plot output of best subclade topologies
	unless ( $href_txt_filenames->{barDb} eq 'none' ){
		
		my	%best_score_of_topo	; # key: best spin topology; value: total score of topology
		for my $t ( @quartets_all ){ $best_score_of_topo{$t} = $sum_value_of_spin_dir{$t} }
		my	@best_spin_trees	=	sort { $best_score_of_topo{$b}<=>$best_score_of_topo{$a} } keys %best_score_of_topo ;
		my	@order_number		=	qw/1 2 3/ ;
		
		my %hoh_counter_order_besttop;
		for my $t ( keys %hol_taxvalues_of_subclade ){
			
			$hoh_counter_order_besttop{$t}{1}	=	exists ( $hol_taxvalues_of_subclade{$t}{0} ) ? ( @{$hol_taxvalues_of_subclade{$t}{0}} / 2 ) : ( 0 ) ;
			$hoh_counter_order_besttop{$t}{2}	=	exists ( $hol_taxvalues_of_subclade{$t}{1} ) ? ( @{$hol_taxvalues_of_subclade{$t}{1}} / 2 ) : ( 0 ) ;
			$hoh_counter_order_besttop{$t}{3}	=	exists ( $hol_taxvalues_of_subclade{$t}{2} ) ? ( @{$hol_taxvalues_of_subclade{$t}{2}} / 2 ) : ( 0 ) ;
		}
		
		&print_svg_barchart_graphic(
								
								\%$href_grafik_colours				,	# key: best spin topologies; value: assigned svg color				In -> Defined, OUT -> Unchanged
								\@best_spin_trees					,	# list of best clade topologies										In -> Defined, OUT -> Unchanged
								\%hoh_counter_order_besttop			,	# key1: quartet tree; key2: order position; value: N occurence		In -> Defined, OUT -> Unchanged
								\$N_final_quartets					,	# total number of quartets											In -> Defined, OUT -> Unchanged
								\$href_graphic_title->{barDb}		,	# svg graphic title													In -> Defined, OUT -> Unchanged
								\$href_txt_filenames->{barDb}		,	# svg outfile name													In -> Defined, OUT -> Unchanged
								\@order_number						,	# order numbers (1 best to 6 lowest support)						In -> Defined, OUT -> Unchanged
								\%$href_graphic_setup				,	# key: svgoption; value: option setup								In -> Defined, OUT -> Unchanged
								\$href_name_subfolder->{svg}		,	# key: subfolder code; value: subfolder name						In -> Defined, OUT -> Unchanged
								\$href_grafik_colours->{back}		,	# background color svg												In -> Defined, OUT -> Unchanged
								\$href_grafik_colours->{header}		,	# header color														In -> Defined, OUT -> Unchanged
								\$href_grafik_colours->{barleft}	,	# barplot color of lower supported spin trees						In -> Defined, OUT -> Unchanged
		);
		
		@best_spin_trees			=	() ;
		%hoh_counter_order_besttop	=	() ;
	}
	################################
	
	################################################################################################ END Bar Chart of Topology Order Distribution
	
	
	
	################################################################################################ START Split Support Box Plot
	
	################################ Mean Split support
	# print svg mean split support graphic
	my	%sorted_value_of_top ;
	my	%N_of_top ;
	my	@sorted_topo_number_by_value	= sort {$mean_of_top{$b}<=>$mean_of_top{$a}} keys %mean_of_top ;
	my	@sorted_topologies_by_value		= ( $quartets_all[$sorted_topo_number_by_value[0]], $quartets_all[$sorted_topo_number_by_value[1]], $quartets_all[$sorted_topo_number_by_value[2]] ) ;
	
	for (0 .. 2){
		
		$sorted_value_of_top{$sorted_topologies_by_value[$_]}	= $mean_of_top{$sorted_topo_number_by_value[$_]} ;
		$N_of_top{$sorted_topologies_by_value[$_]}				= $mean_of_top{$sorted_topo_number_by_value[$_]} ;
	}
	
	unless ( $href_txt_filenames->{splimean} eq 'none' ){
		
		&print_svg_split_graphic(
								
								\%$href_grafik_colours					,	# key: topology support (0,1, or 2); value: assigned svg color								In -> Defined, OUT -> Unchanged
								\@sorted_topologies_by_value			,	# list of sorted clade topologies															In -> Defined, OUT -> Unchanged
								\$href_txt_filenames->{splimean}		,	# Output svg name																			In -> Defined, OUT -> Unchanged
								\$href_graphic_title->{splimean}		,	# Graphic header text																		In -> Defined, OUT -> Unchanged
								\%sorted_value_of_top					,	# key: toponumber (0, 1, or 2); value: mean split score										In -> Defined, OUT -> Unchanged
								\%$href_graphic_setup					,	# key: svgoption; value: option setup														In -> Defined, OUT -> Unchanged
								\%N_of_top								,	# key: toponumber (0, 1, or 2); value: mean split score, only important in	split Number	In -> Defined, OUT -> Unchanged
								\'Mean Split'							,	# value code of figure legend print out														In -> Defined, OUT -> Unchanged
								\%q_number_of_top						,	# key alphabetically sorted topology, value: Q code											In -> Defined, OUT -> Unchanged
								\$href_name_subfolder->{svg}			,	# key: subfolder code; value: subfolder name												In -> Defined, OUT -> Unchanged
								\$href_parameter_setting_of_option->{r}	,	# value (0 -> no directive multiple clan analysis, 1 -> spin directive)						In -> Defined, OUT -> Unchanged
								\$href_grafik_colours->{back}			,	# background color svg																		IN -> Defined, OUT -> Unchanged
								\$href_grafik_colours->{header}			,	# header color																				IN -> Defined, OUT -> Unchanged
							) ;
	}
	################################
	
	
	
	################################ Median Split support
	# print svg median split support graphic
	%sorted_value_of_top			=	() ;
	@sorted_topo_number_by_value	=	sort {$median_of_top{$b}<=>$median_of_top{$a}} keys %median_of_top ;
	@sorted_topologies_by_value		= ( $quartets_all[$sorted_topo_number_by_value[0]], $quartets_all[$sorted_topo_number_by_value[1]], $quartets_all[$sorted_topo_number_by_value[2]] ) ;
	
	for (0 .. 2){
		
		$sorted_value_of_top{$sorted_topologies_by_value[$_]}	= $median_of_top{$sorted_topo_number_by_value[$_]} ;
		$N_of_top{$sorted_topologies_by_value[$_]}						= $median_of_top{$sorted_topo_number_by_value[$_]}
	}
	
	unless ( $href_txt_filenames->{splimedi} eq 'none' ){
		
		&print_svg_split_graphic(
								
								\%$href_grafik_colours					,	# key: topology support (0,1, or 2); value: assigned svg color								In -> Defined, OUT -> Unchanged
								\@sorted_topologies_by_value			,	# list of sorted clade topologies															In -> Defined, OUT -> Unchanged
								\$href_txt_filenames->{splimedi}		,	# Output svg name																			In -> Defined, OUT -> Unchanged
								\$href_graphic_title->{splimedi}		,	# Graphic header text																		In -> Defined, OUT -> Unchanged
								\%sorted_value_of_top					,	# key: toponumber (0, 1, or 2); value: mean split score										In -> Defined, OUT -> Unchanged
								\%$href_graphic_setup					,	# key: svgoption; value: option setup														In -> Defined, OUT -> Unchanged
								\%N_of_top								,	# key: toponumber (0, 1, or 2); value: mean split score, only important in	split Number	In -> Defined, OUT -> Unchanged
								\'Median Split'							,	# value code of figure legend print out														In -> Defined, OUT -> Unchanged
								\%q_number_of_top						,	# key alphabetically sorted topology, value: Q code											In -> Defined, OUT -> Unchanged
								\$href_name_subfolder->{svg}			,	# key: subfolder code; value: subfolder name												In -> Defined, OUT -> Unchanged
								\$href_parameter_setting_of_option->{r}	,	# value (0 -> no directive multiple clan analysis, 1 -> spin directive)						In -> Defined, OUT -> Unchanged
								\$href_grafik_colours->{back}			,	# background color svg																		IN -> Defined, OUT -> Unchanged
								\$href_grafik_colours->{header}			,	# header color																				IN -> Defined, OUT -> Unchanged
							 ) ;
	}
	################################
	
	################################################################################################ END Split Support Box Plot
	
	
	
	################################
	# print txt info about resolved quartets
	unless ( ( @$aref_quartets == 1 ) || ( $href_txt_filenames->{txtDb} eq 'none' ) ){
		
		my	%best_score_of_topo	; # key: best spin topology; value: total score of topology
		for my $t ( @quartets_all ){ $best_score_of_topo{$t} = $sum_value_of_spin_dir{$t} }
		my	@sorted_best_trees	=	sort { $best_score_of_topo{$b}<=>$best_score_of_topo{$a} } keys %best_score_of_topo ;
		
		&print_txt_quartet_results(
								
								\%hol_taxvalues_of_subclade			,	# key1: subclade, key2: support code (0,1), value: list quartet taxa, split support value		In -> Defined, OUT -> Unchanged
								\$href_txt_filenames->{txtDb}		,	# text outfile name																				In -> Defined, OUT -> Unchanged
								\$href_name_subfolder->{txt}		,	# key: subfolder code; value: subfolder name													In -> Defined, OUT -> Unchanged
								\@sorted_best_trees					,	# list of best topologies, sorted by scoring													In -> Defined, OUT -> Unchanged
		) ;
		
		@sorted_best_trees	=	() ;
	}
	################################
	
	if ( ( $href_parameter_setting_of_option->{N} == 0 ) && ( $href_txt_filenames->{txtrejected} ne 'none' ) ){ close	OUTrej	}
	
	unless ( $href_txt_filenames->{txtrejected} eq 'none' ){ print "\n\tPrint TXT ", $href_txt_filenames->{txtrejected} }
	unless ( $href_txt_filenames->{txttsupport} eq 'none' ){ print "\n\tPrint TXT ", $href_txt_filenames->{txttsupport} }
	unless ( $href_txt_filenames->{txtqsupport} eq 'none' ){ print "\n\tPrint TXT ", $href_txt_filenames->{txtqsupport} }
	unless ( $href_txt_filenames->{txtaddprint} eq 'none' ){ print "\n\tPrint TXT ", $href_txt_filenames->{txtaddprint} }
	unless ( $href_txt_filenames->{txtssupport} eq 'none' ){ print "\n\tPrint TXT ", $href_txt_filenames->{txtssupport} }
	#################################################################################################################### END Result Print OUT
}

sub non_directive_clan_evaluation{
	
	my	$aref_split_topologies					=	$_[0]	;	#	list of split topologies (1,2,3)																										In -> Defined, OUT -> Unchanged
	my	$href_spin_value_of_quartet_topology	=	$_[1]	;	#	key: quartet_topology; value: list of spin_values																						In -> Defined, OUT -> Changed
	my	$href_hoh_seen_quartet					=	$_[2]	;	#	key1: support number (0,1,2 -> 0 element: best supported quartet ), key2: topology; value: count number									In -> Undefined, OUT -> defined
	my	$href_hol_taxvalues_of_subclade			=	$_[3]	;	#	key1: subclade topology; key2: support code (0->best, 1-> 2nd best); value1: sorted taxonnames (comma separated); value2: split score	In -> Undefined, OUT -> defined
	my	$href_hol_val_of_top_of_tax				=	$_[4]	;	#	key1: subclade topology; key2: taxonname; value: list of taxon observed split scores													In -> Undefined, OUT -> defined
	my	$href_rejected_quartet_of_number		=	$_[5]	;	#	key1: quartet number; value: counter of rejection (0 or 1)																				In -> Defined, OUT -> Unchanged
	my	$sref_out_symetric						=	$_[6]	;	#	name summarized split support, former M-ICE symetric outfile																			In -> Defined, OUT -> Unchanged
	my	$sref_msa_infile_name					=	$_[7]	;	#	msa infile name																															In -> Defined, OUT -> Unchanged
	my	$sref_name_folder						=	$_[8]	;	#	key: subfolder code; value: subfolder name																								In -> Defined, OUT -> Unchanged
	my	$href_seen_taxa							=	$_[9]	;	#	key: sequence name; value: N occurence in quartet analyses																				In -> Undefined, OUT -> defined
	my	$sref_norm_level						= $_[10]; # if 1 -> subsequent normalisation in respect of the lowest topology score (set to zero)															In -> Defined, OUT -> Unchanged
	
	
	print "\n\t---------------------------------------------------\n\n\tSplit Result(s) of Single Quartets...\n" ;
	
	
	#############################
	# (1) Identification of highest split score due to different spin directions of the same quartet tree
	# (2) Subsequent normalisation in respect of the lowest topology score (set to zero)
	# (3) Calculation posterior probability for each split signal n (Signal=(n/signal q1+q2+3)
	# (4) Counting of best, second best and lowest supported quartet tree
	# (5) Assign best, second best and lowest support to taxa of each quartet analysis
	# (6) 4BASAL print out ONLY if Outfilename is not 'none' (see header outfile declaration) !!!!!!
	my %best_value_of_quartet;
	my $counter_quartet = 0 ;
	
	for my $quartet_number ( 1 .. @$aref_split_topologies ){ #print "\nq number: ", $quartet_number;
		
		unless ( $href_rejected_quartet_of_number->{$quartet_number} ){
			
			my @taxa	=	split ":", $aref_split_topologies->[$quartet_number-1] ;
			
			my (
					$total_score ,
					%highest_spin_of_topology ,
					%seen_topology ,
			) ;
			
			# (1)
			for my $topology_dir ( sort keys %$href_spin_value_of_quartet_topology ){ #print "\n", $topology_dir ;
				
				(	my	$top		=	$topology_dir )	=~	s/\),\(/\)::\(/ ;
					my	@sisters	=	split "::", $top ;
						@sisters	=	sort @sisters ;
						$top		=	join ",", @sisters ; #print "\n", $top ;
				
				unless ( $seen_topology{$top} ){#print "\n unseen:", $top , "\tscore: ", $href_spin_value_of_quartet_topology->{$topology_dir}[$counter_quartet];
					
					$seen_topology{$top}++ ;
					
					$highest_spin_of_topology{$top} = $href_spin_value_of_quartet_topology->{$topology_dir}[$counter_quartet] ; #print "\t spin: ", $highest_spin_of_topology{$top} ;
				}
				
				elsif ( $highest_spin_of_topology{$top} < $href_spin_value_of_quartet_topology->{$topology_dir}[$counter_quartet] ){#print "\n seen:", $top ;
					
					$highest_spin_of_topology{$top} = $href_spin_value_of_quartet_topology->{$topology_dir}[$counter_quartet] ; #print "\t spin: ", $highest_spin_of_topology{$top} ;
				}
			}
			
			# (2) subsequent normalization in respect of the lowest best toopology score (set to zero) if -k option equal 0
			my @sorted_spin_topos	= sort {$highest_spin_of_topology{$b}	<=> $highest_spin_of_topology{$a}}	keys %highest_spin_of_topology ;
			if ( $$sref_norm_level == 0 ){
				
				if ( $highest_spin_of_topology{$sorted_spin_topos[2]} > 0 ){
					
					$highest_spin_of_topology{$sorted_spin_topos[0]} -= abs($highest_spin_of_topology{$sorted_spin_topos[2]}) ;
					$highest_spin_of_topology{$sorted_spin_topos[1]} -= abs($highest_spin_of_topology{$sorted_spin_topos[2]}) ;
					$highest_spin_of_topology{$sorted_spin_topos[2]} -= abs($highest_spin_of_topology{$sorted_spin_topos[2]})
				}
				
				elsif ( $highest_spin_of_topology{$sorted_spin_topos[2]} < 0 ){
					
					$highest_spin_of_topology{$sorted_spin_topos[0]} += abs($highest_spin_of_topology{$sorted_spin_topos[2]}) ;
					$highest_spin_of_topology{$sorted_spin_topos[1]} += abs($highest_spin_of_topology{$sorted_spin_topos[2]}) ;
					$highest_spin_of_topology{$sorted_spin_topos[2]} += abs($highest_spin_of_topology{$sorted_spin_topos[2]})
				}
			}
			
			#exit;
			# (3)
			$total_score = $highest_spin_of_topology{$sorted_spin_topos[0]} + $highest_spin_of_topology{$sorted_spin_topos[1]} + $highest_spin_of_topology{$sorted_spin_topos[2]} ;
			unless ( $total_score > 0 ){ die "BUG-ERROR: Cannot calculate normalised split scorees. Total number of split score equal zero in subroutine &non_directive_clan_evaluation!\n\tPlease, report BUG to system developer!\n" }
			
			for my $top ( keys %highest_spin_of_topology ){
				
				push @{$best_value_of_quartet{$top}}, ( $highest_spin_of_topology{$top} / $total_score ) ;
				
				for my $taxon ( @taxa ){ push @{$href_hol_val_of_top_of_tax->{$top}{$taxon}}, $best_value_of_quartet{$top}[$counter_quartet]; $href_seen_taxa->{$taxon}++ }
			}
			
			# (4)
			$href_hoh_seen_quartet->{0}{$sorted_spin_topos[0]}++ ;
			$href_hoh_seen_quartet->{1}{$sorted_spin_topos[1]}++ ;
			$href_hoh_seen_quartet->{2}{$sorted_spin_topos[2]}++ ;
			
			for my $t ( @sorted_spin_topos ){ print "\n\tq", $quartet_number ,": Best score ", $t , ":\t", $best_value_of_quartet{$t}[$counter_quartet] }
			print "\n" ;
			
			# (5)
			push @{$href_hol_taxvalues_of_subclade->{$sorted_spin_topos[0]}{0}} , ( $aref_split_topologies->[$quartet_number-1], $best_value_of_quartet{$sorted_spin_topos[0]}[$counter_quartet] ) ;
			push @{$href_hol_taxvalues_of_subclade->{$sorted_spin_topos[1]}{1}} , ( $aref_split_topologies->[$quartet_number-1], $best_value_of_quartet{$sorted_spin_topos[1]}[$counter_quartet] ) ;
			push @{$href_hol_taxvalues_of_subclade->{$sorted_spin_topos[2]}{2}} , ( $aref_split_topologies->[$quartet_number-1], $best_value_of_quartet{$sorted_spin_topos[2]}[$counter_quartet] ) ;
			
			
			# (6)
			unless ( $$sref_out_symetric eq 'none'){
				
				open	OUTpipe,	">>$$sref_name_folder/$$sref_out_symetric" or die "OUTFILE-ERROR: Cannot open info outfile ", $$sref_name_folder, "/", $$sref_out_symetric ," in subroutine &non_directive_clan_evaluation!\n" ;
				print	OUTpipe		$$sref_msa_infile_name, "\t", $quartet_number, ":\t", $aref_split_topologies->[$quartet_number-1],"\t",
												$sorted_spin_topos[0], "\t", $best_value_of_quartet{$sorted_spin_topos[0]}[$counter_quartet], "\t",
												$sorted_spin_topos[1], "\t", $best_value_of_quartet{$sorted_spin_topos[1]}[$counter_quartet], "\t",
												$sorted_spin_topos[2], "\t", $best_value_of_quartet{$sorted_spin_topos[2]}[$counter_quartet], "\n",
				;
				close	OUTpipe
			}
			
			
			$counter_quartet++
		}
		
		else { print "\n\tq", $quartet_number ," rejected!\n" }
	}
	#############################
	
	
	
	#############################
	# For each quartet tree delete lower spin scores
	%$href_spin_value_of_quartet_topology	= %best_value_of_quartet ;
	%best_value_of_quartet					= () ;
	
	print "\n\t---------------------------------------------------\n" ;
	#############################
}

sub spin_directive_clan_evaluation{
	
	my	$aref_split_topologies					=	$_[0]	;	#	list of split topologies (1,2,3)																										In -> Defined, OUT -> Unchanged
	my	$href_sum_value_of_spin_dir				=	$_[1]	;	#	key: spin directive clan topology; value: total score after single quartet analyses														In -> Defined, OUT -> Unchanged
	my	$href_spin_value_of_quartet_topology	=	$_[2]	;	#	key1: threshold code, key2: quartet topo, value: quartet weight																			In -> Defined, OUT -> Changed
	my	$href_hoh_seen_quartet					=	$_[3]	;	#	key1: support number (0,1,2 -> 0 element: best supported quartet ), key2: topology; value: count number									In -> Undefined, OUT -> defined
	my	$href_hol_taxvalues_of_subclade			=	$_[4]	;	#	key1: subclade topology; key2: support code (0->best, 1-> 2nd best); value1: sorted taxonnames (comma separated); value2: split score	In -> Undefined, OUT -> defined
	my	$href_hol_val_of_top_of_tax				=	$_[5]	;	#	key1: subclade topology; key2: taxonname; value: list of taxon observed split scores													In -> Undefined, OUT -> defined
	my	$href_rejected_quartet_of_number		=	$_[6]	;	#	key1: quartet number; value: counter of rejection (0 or 1)																				In -> Defined, OUT -> Unchanged
	my	$sref_out_symetric						=	$_[7]	;	#	name summarized split support, former M-ICE symetric outfile																			In -> Defined, OUT -> Unchanged
	my	$sref_msa_infile_name					=	$_[8]	;	#	msa infile name																															In -> Defined, OUT -> Unchanged
	my	$sref_name_folder						=	$_[9]	;	#	key: subfolder code; value: subfolder name																								In -> Defined, OUT -> Unchanged
	my	$href_seen_taxa							=	$_[10]	;	#	key: sequence name; value: N occurence in quartet analyses																				In -> Undefined, OUT -> defined
	my	$sref_norm_level						=	$_[11]	;	#	if 1 -> subsequent normalisation in respect of the lowest topology score (set to zero)													In -> Defined, OUT -> Unchanged
	
	print "\n\t---------------------------------------------------\n\n\tSplit Result(s) of Single Quartets...\n" ;
	
	
	#############################
	# Sort highest split signal found for each quartet in descending order
	my @sorted_spin_topos	= sort {$href_sum_value_of_spin_dir->{$b}	<=> $href_sum_value_of_spin_dir->{$a}}	keys %$href_sum_value_of_spin_dir ;
	#############################
	
	
	
	#############################
	# Identify highest spin signal of each of the three spin quartet topologies
	# delete spin scores of other spin quartet topologies
	my	%seen_topology	;
	
	for my $spin_topo ( @sorted_spin_topos ){
		
		(	my	$top	= $spin_topo )		=~	s/\),\(/\)::\(/ ; #print "\n", $top ;
			my	@sisters					=	split "::", $top ;
			
		unless (	( $seen_topology{$sisters[0]} ) &&
					( $seen_topology{$sisters[1]} ) ){
			
			$seen_topology{$sisters[0]}++ ; #print "\nseen: ", $spin_topo;
			$seen_topology{$sisters[1]}++ ;
		}
		else { delete $href_spin_value_of_quartet_topology->{$spin_topo} }
	}
	#############################
	
	
	
	#############################
	# (1) subsequent normalization in respect of the lowest best toopology score (set to zero) if -k option equal 0
	# (2) Calculation posterior probability for each split signal n (Signal=(n/signal q1+q2+3)
	# (3) Counting of best, second best and lowest supported quartet tree
	# (4) Assign best, second best and lowest support to taxa of each quartet analysis
	my $counter_quartet = 0 ;
	
	for my $quartet_number ( 1 .. @$aref_split_topologies ){
		
		unless ( $href_rejected_quartet_of_number->{$quartet_number} ){
			
			my @taxa	=	split ":", $aref_split_topologies->[$quartet_number-1] ;
			
			my (
					$total_score ,
					%highest_spin_of_topology ,
			) ;
			
			# (1)
			if ( $$sref_norm_level == 0 ){
				
				my %value_of_top ;
				for my $topology_dir ( keys %$href_spin_value_of_quartet_topology ){ $value_of_top{$topology_dir} = $href_spin_value_of_quartet_topology->{$topology_dir}[$counter_quartet] }
				my @sorted_spin_topos	= sort {$value_of_top{$b}	<=> $value_of_top{$a}}	keys %value_of_top ;
				
				if ( $value_of_top{$sorted_spin_topos[2]} > 0 ){
					
					$href_spin_value_of_quartet_topology->{$sorted_spin_topos[0]}[$counter_quartet] -= abs($value_of_top{$sorted_spin_topos[2]}) ;
					$href_spin_value_of_quartet_topology->{$sorted_spin_topos[1]}[$counter_quartet] -= abs($value_of_top{$sorted_spin_topos[2]}) ;
					$href_spin_value_of_quartet_topology->{$sorted_spin_topos[2]}[$counter_quartet] -= abs($value_of_top{$sorted_spin_topos[2]})
				}
				
				elsif ( $value_of_top{$sorted_spin_topos[2]} < 0 ){
					
					$href_spin_value_of_quartet_topology->{$sorted_spin_topos[0]}[$counter_quartet] += abs($value_of_top{$sorted_spin_topos[2]}) ;
					$href_spin_value_of_quartet_topology->{$sorted_spin_topos[1]}[$counter_quartet] += abs($value_of_top{$sorted_spin_topos[2]}) ;
					$href_spin_value_of_quartet_topology->{$sorted_spin_topos[2]}[$counter_quartet] += abs($value_of_top{$sorted_spin_topos[2]})
				}
			}
			
			# (2)
			for my $topology_dir ( keys %$href_spin_value_of_quartet_topology ){ $total_score += $href_spin_value_of_quartet_topology->{$topology_dir}[$counter_quartet] } #print "\nscore total: ", $total_score;
			unless ( $total_score > 0 ){ die "BUG-ERROR: Cannot calculate normalised split scorees. Total number of split score equal zero in subroutine &spin_directive_clan_evaluation!\n\tPlease, report BUG to system developer!\n" }
			
			for my $topology_dir ( keys %$href_spin_value_of_quartet_topology ){ #print "\nquartet: ", $topology_dir, "\t", $href_spin_value_of_quartet_topology->{$topology_dir}[$counter_quartet] ;
				
				$href_spin_value_of_quartet_topology->{$topology_dir}[$counter_quartet] = (	$href_spin_value_of_quartet_topology->{$topology_dir}[$counter_quartet] / $total_score ) ; #print "\n norm: ", $href_spin_value_of_quartet_topology->{$topology_dir}[$counter_quartet];
				$highest_spin_of_topology{$topology_dir}																=		$href_spin_value_of_quartet_topology->{$topology_dir}[$counter_quartet] ; #print "\n highest ", $href_spin_value_of_quartet_topology->{$topology_dir}[$counter_quartet];
				
				for my $t ( @taxa ){ push @{$href_hol_val_of_top_of_tax->{$topology_dir}{$t}}, $highest_spin_of_topology{$topology_dir}; $href_seen_taxa->{$t}++ }
			}
			
			# (3)
			my @sorted_spin_topos	= sort {$highest_spin_of_topology{$b}	<=> $highest_spin_of_topology{$a}}	keys %highest_spin_of_topology ;
			$href_hoh_seen_quartet->{0}{$sorted_spin_topos[0]}++ ;
			$href_hoh_seen_quartet->{1}{$sorted_spin_topos[1]}++ ;
			$href_hoh_seen_quartet->{2}{$sorted_spin_topos[2]}++ ;
			
			for my $t ( @sorted_spin_topos ){ print "\n\tq", $quartet_number ,": Best score ", $t , ":\t", $highest_spin_of_topology{$t} }
			print "\n" ;
			
			# (4)
			push @{$href_hol_taxvalues_of_subclade->{$sorted_spin_topos[0]}{0}} , ( $aref_split_topologies->[$quartet_number-1], $highest_spin_of_topology{$sorted_spin_topos[0]} ) ;
			push @{$href_hol_taxvalues_of_subclade->{$sorted_spin_topos[1]}{1}} , ( $aref_split_topologies->[$quartet_number-1], $highest_spin_of_topology{$sorted_spin_topos[1]} ) ;
			push @{$href_hol_taxvalues_of_subclade->{$sorted_spin_topos[2]}{2}} , ( $aref_split_topologies->[$quartet_number-1], $highest_spin_of_topology{$sorted_spin_topos[2]} ) ;
			
			# (5)
			unless ( $$sref_out_symetric eq 'none'){
				
				open	OUTpipe,	">>$$sref_name_folder/$$sref_out_symetric" or die "OUTFILE-ERROR: Cannot open info outfile ", $$sref_name_folder, "/", $$sref_out_symetric ," in subroutine &spin_directive_clan_evaluation!\n" ;
				print	OUTpipe		$$sref_msa_infile_name, "\tq", $quartet_number, ":\t", $aref_split_topologies->[$quartet_number-1],"\t",
									$sorted_spin_topos[0], "\t", $href_spin_value_of_quartet_topology->{$sorted_spin_topos[0]}[$counter_quartet], "\t",
									$sorted_spin_topos[1], "\t", $href_spin_value_of_quartet_topology->{$sorted_spin_topos[1]}[$counter_quartet], "\t",
									$sorted_spin_topos[2], "\t", $href_spin_value_of_quartet_topology->{$sorted_spin_topos[2]}[$counter_quartet], "\n",
				;
				close	OUTpipe
			}
			
			$counter_quartet++
		}
		
		else { print "\n\tq", $quartet_number ," rejected!" }
	}
	#############################
	
	#############################
	# print total spin score info
	#for my $t ( @sorted_spin_topos ){print "\nt", $t, "\t", $href_sum_value_of_spin_dir->{$t}}
	print "\n\t---------------------------------------------------\n" ;
	#############################
}

sub recode_pattern{

	# &recode_pattern
	# my $recoded_pattern = &recode_pattern(\@pattern, \$state_handling, \state_sequence) ;

	my $aref_pattern_states		= $_[0] ; #
	my $sref_recoding_class_key	= $_[1] ; #
	my $sref_sequence_state		= $_[2] ; #

	#######################
	# If msa pattern should be recoded after state classes (like purine or pyrimidine if nucleotide data)
	# substitute single character states of list @$aref_pattern_states to corresponding class states
	if ( $$sref_recoding_class_key == 1 ){

		if ( $$sref_sequence_state eq 'nuc' ){

			map { s/A|G/P/; $_ } @$aref_pattern_states;	# substitute A and G to P (Purine)
			map { s/C|T/Y/; $_ } @$aref_pattern_states	# substitute C and T to Y (Pyrimidine)
		}

		elsif ( $$sref_sequence_state eq 'aa' ){

			map { s/A|W|M|I|L|F|P/H/; $_			} @$aref_pattern_states;	# substitute A|W|M|I|L|F|P to H (Hydrophobic)
			map { s/C|G|T|N|Y|R|S|K|D|V|E|Q/Y/; $_	} @$aref_pattern_states;	# substitute C|G|T|N|Y|R|S|K|D|V|E|Q to Y (Hydrophilic)
		}

		else{ die "\n\tBUG-ERROR: cannot assign data type in subroutine &recode_pattern!\n\tPlease, report BUG to system developer!\n" }
	}
	#######################

	##############################################
	# START quartet pattern recoding

	#######################
	# definition of single quartet pattern codes
	# key state pattern; value: quartet pattern code
	my %code_of_pattern = (
							'XXYY' => 'A',
							'XYXY' => 'B',
							'XYYX' => 'C',
							'XXYW' => 'D',
							'XYWW' => 'E',
							'XYWY' => 'F',
							'XYXW' => 'G',
							'XYYW' => 'H',
							'XYWX' => 'J',
							'XYYY' => 'K',
							'XYXX' => 'L',
							'XXYX' => 'M',
							'XXXY' => 'N',
							'XYWZ' => 'V',
							'XXXX' => 'I'
						) ;
	#######################

	#######################
	# START quartet pattern recoding via if queries
	if ( $aref_pattern_states->[0] eq $aref_pattern_states->[1] ){

		if ( $aref_pattern_states->[1] eq $aref_pattern_states->[2] ){

			if ( $aref_pattern_states->[2] eq $aref_pattern_states->[3] )	{ return $code_of_pattern{XXXX} }	# I
			else 																													{ return $code_of_pattern{XXXY} }	# N
		}

		elsif( $aref_pattern_states->[1] eq $aref_pattern_states->[3] )	{ return $code_of_pattern{XXYX} }	# M
		elsif( $aref_pattern_states->[2] eq $aref_pattern_states->[3] )	{ return $code_of_pattern{XXYY} }	# A
		else																														{ return $code_of_pattern{XXYW} }	# D
	}

	elsif( $aref_pattern_states->[0] eq $aref_pattern_states->[2] ){

		if		( $aref_pattern_states->[2] eq $aref_pattern_states->[3] )	{ return $code_of_pattern{XYXX} }	# L
		elsif	( $aref_pattern_states->[1] eq $aref_pattern_states->[3] )	{ return $code_of_pattern{XYXY} }	# B
		else																															{ return $code_of_pattern{XYXW} }	# G
	}

	elsif( $aref_pattern_states->[0] eq $aref_pattern_states->[3] ){

		if		( $aref_pattern_states->[1] eq $aref_pattern_states->[2] )	{ return $code_of_pattern{XYYX} }	# C
		else																															{ return $code_of_pattern{XYWX} }	# J
	}

	elsif( $aref_pattern_states->[1] eq $aref_pattern_states->[2] ){

		if		( $aref_pattern_states->[2] eq $aref_pattern_states->[3] )	{ return $code_of_pattern{XYYY} }	# K
		else																															{ return $code_of_pattern{XYYW} }	# H
	}

	elsif( $aref_pattern_states->[2] eq $aref_pattern_states->[3] )		{ return $code_of_pattern{XYWW} }	# E

	elsif( $aref_pattern_states->[1] eq $aref_pattern_states->[3] )		{ return $code_of_pattern{XYWY} }	# F

	else 																															{ return $code_of_pattern{XYWZ} }	# V
	# END quartet pattern recoding via if queries
	#######################
}

sub create_msa{
	
	my $href_msasequence_of_taxon	=	$_[0] ;	#
	my $sref_ml_method_infile		=	$_[1] ;	#
	
	#######################
	# Print OUT of masked/unmasked sequences in phylip format for ML pattern estimation
					$$sref_ml_method_infile = "Penguin_quartet_msa.phy" ;
	open OUT,		">$$sref_ml_method_infile" or die "OUTFILE-ERROR: Cannot open MSA outfile ", $$sref_ml_method_infile, " in subroutine &create_msa!\n";
	
	my @taxa	=	sort keys %$href_msasequence_of_taxon ;
	my $N_taxa	=	@taxa ;
	my $L_seq	=	length $href_msasequence_of_taxon->{$taxa[0]} ;
	
	print OUT "\t", $N_taxa, "\t", $L_seq, "\n" ;
	for my $taxon (  @taxa ){ print OUT $taxon, "   ", $href_msasequence_of_taxon->{$taxon}, "\n" }
	close OUT;
	#######################
}

sub write_p4_script{

	# &write_p4_script( \%$href_masked_sequence_of_taxon, \$ml_infile_name, \$split_topology_1, \$split_topology_2, \$split_topology_3, \$p4_infile_name ) ;

	my	$sref_msa_infile	=	$_[0]	;	# p4 alignment infile				-> IN (not changed)
	my	$sref_topo_1		=	$_[1]	;	# quartet topology (tax1,tax2)		-> IN (not changed)
	my	$sref_topo_2		=	$_[2]	;	# quartet topology (tax1,tax3)		-> IN (not changed)
	my	$sref_topo_3		=	$_[3]	;	# quartet topology (tax1,tax4)		-> IN (not changed)
	my	$sref_gamma			=	$_[4]	;	# gamma start value					-> IN (not changed)
	my	$sref_pinv			=	$_[5]	;	# pinv start value					-> IN (not changed)
	my	$sref_model			=	$_[6]	;	# subst. model						-> IN (not changed)
	my	$sref_p4_infile		=	$_[7]	;	# name of the p4 infile				-> OUT (defined)
	my	$sref_seq_type		=	$_[8]	;	# sequence type						-> IN (not changed)



	##############################################
	# Change rooted quartet tree to unrooted tree
	# ((A,B),(C,D)) -> (A,B,(C,D)) should make P4 analysis a little bit faster
	(	my	$p4_tree1	=	$$sref_topo_1	)	=~	s/\(\(/\(/	;	$p4_tree1	=~	s/\),/,/	;
	(	my	$p4_tree2	=	$$sref_topo_2	)	=~	s/\(\(/\(/	;	$p4_tree2	=~	s/\),/,/	;
	(	my	$p4_tree3	=	$$sref_topo_3	)	=~	s/\(\(/\(/	;	$p4_tree3	=~	s/\),/,/	;
	##############################################



	##############################################
	# ML model parameters
	# später ersetzen durch string referenzen, da diese werte übers terminal vom user definiert werden

	my	$ncat	= 4		;

	my %p4_parameters = (

		'gamma_cat'	=>	$ncat											,
		'pinv'		=>	$$sref_pinv										,
		'pinvfree'	=>	1												,
		'alpha'		=>	$$sref_gamma									,
		'asrv'		=>	"    t.newGdasrv(free=1, val=".$$sref_gamma.")"	,
		'model'		=>	$$sref_model									,
	) ;

	if ( $$sref_gamma == 100 ){

		$p4_parameters{asrv}		=	()		;
		$p4_parameters{pinv}		=	'0.0'	;
		$p4_parameters{pinvfree}	=	'0'		;
		$p4_parameters{gamma_cat}	=	1		;
	}
	##############################################



	##############################################
	# Definement python code ML model parameters
	# JC needs another input code as e.g. GTR
	my	(
			$newcomp	,
			$rmatrix	,
	) ;

	if		( $$sref_seq_type	eq	'nuc' ){

		if		( $p4_parameters{model}	eq	'JC' ){

			$newcomp	=	"    t.newComp(free=0, spec=\'equal\')" ;
			$rmatrix	=	"    t.newRMatrix(free=0, spec=\'ones\')" ;
		}

		elsif	( $p4_parameters{model}	eq	'F81' ){

			$newcomp	=	"    t.newComp(free=1, spec=\'empirical\')" ;
			$rmatrix	=	"    t.newRMatrix(free=0, spec=\'ones\')" ;

		}

		elsif	( $p4_parameters{model}	eq	'K2P' ){

			$newcomp	=	"    t.newComp(free=0, spec=\'equal\')" ;
			$rmatrix	=	"    t.newRMatrix(free=1, spec=\'2p\')" ;
		}

		elsif	( $p4_parameters{model}	eq 'HKY' ){

			$newcomp	=	"    t.newComp(free=1, spec=\'empirical\')" ;
			$rmatrix	=	"    t.newRMatrix(free=1, spec=\'2p\')" ;
		}

		elsif	( $p4_parameters{model}	eq	'GTR' ){

			$newcomp	=	"    t.newComp(free=1, spec=\'empirical\')" ;
			$rmatrix	=	"    t.newRMatrix(free=1, spec=\'specified\', val=[2., 3., 4., 5., 6., 7.])" ;
		}

		else{ die "\n\tBUG-ERROR: Cannot assign defined ", $$sref_seq_type ," ML model ", $p4_parameters{model}, " in subroutine &write_p4_script!\n\tPlease, report BUG to system developer!\n\t" }

		$p4_parameters{seqtype}		=	'dna'	;
		$p4_parameters{npattern}	=	'256'	;
	}

	elsif	( $$sref_seq_type	eq	'aa' ){

			$newcomp	=	"    t.newComp(free=1, spec=\'empirical\')" ;
			$rmatrix	=	"    t.newRMatrix(free=0, spec=\'".$$sref_model."\')" ;

			$p4_parameters{seqtype}		=	'protein'	;
			$p4_parameters{npattern}	=	'160000'	;
	}

	else	{ die "\n\tBUG-ERROR: Cannot assign sequence type ", $$sref_seq_type, " in subroutine &write_p4_script!\n\tPlease, report BUG to system developer!\n\t" }
	##############################################



	##############################################
	# open p4 input script
	# print script code
	$$sref_p4_infile	= "p4_temp.py" ;
	open	OUTp4,	">$$sref_p4_infile" or die "OUTFILE-ERROR: Cannot open p4 inputfile ", $$sref_p4_infile," in subroutine &write_p4_script!\n" ;

	print	OUTp4	"var.warnReadNoFile = False\n",
					"read(\'", $$sref_msa_infile  ,"\')\n",
					"a = var.alignments[0]\n" ,
					"d = Data([a])\n" ,
					"aLength = len( var.alignments[0] )\n" ,
					"print (aLength)\n" ,
					"\n" ,
					"myTreeString = \"\"\"\n" ,
					$p4_tree1, ";\n" ,
					$p4_tree2, ";\n" ,
					$p4_tree3, ";\n" ,
					"\"\"\"\n" ,
					"read(myTreeString)\n" ,
					"\n" ,
					"for t in var.trees:\n" ,
					"    \n" ,
					"    t.write(), \"\\n\"\n" ,
					"\n" ,
					"    t.taxNames = a.taxNames\n" ,
					"    t.data = d\n" ,
					"\n" ,
					$newcomp, "\n" ,
					$rmatrix, "\n" ,
					"    t.setNGammaCat(nGammaCat=",								$p4_parameters{gamma_cat}	, ")\n" ,
					$p4_parameters{asrv}, "\n",
					"    t.setPInvar(free=", $p4_parameters{pinvfree},", val=",		$p4_parameters{pinv}		, ")\n" ,
					"\n" ,
					"    t.optLogLike()\n" ,
					"    a2 = func.newEmptyAlignment(dataType=\'", $p4_parameters{seqtype} ,"\', taxNames=a.taxNames, length=", $p4_parameters{npattern}, ")\n" ,
					"\n" ,
					"    for s in a2.sequences:\n" ,
					"        s.sequence = list(s.sequence)\n" ,
					"\n" ,
					"    posn = 0\n" ,
					"    nn = list(a2.symbols)\n" ,
					"    for p0 in nn:\n" ,
					"        for p1 in nn:\n" ,
					"            for p2 in nn:\n" ,
					"                for p3 in nn:\n" ,
					"                    a2.sequences[0].sequence[posn] = p0\n" ,
					"                    a2.sequences[1].sequence[posn] = p1\n" ,
					"                    a2.sequences[2].sequence[posn] = p2\n" ,
					"                    a2.sequences[3].sequence[posn] = p3\n" ,
					"                    posn += 1\n" ,
					"\n" ,
					"    for s in a2.sequences:\n" ,
					"        s.sequence = ''.join(s.sequence)\n" ,
					"\n" ,
					"    d2 = Data([a2])\n" ,
					"    t.data = d2\n" ,
					"    t.calcLogLike()\n" ,
					"    t.getSiteLikes()\n" ,
					"    print (sum(t.siteLikes))\n" ,
					"    sDict = {}\n" ,
					"    sList = []\n" ,
					"\n" ,
					"    for posn in range(", $p4_parameters{npattern}, "):\n" ,
					"        sl = ''.join(a2.sequenceSlice(posn))\n" ,
					"        sDict[sl] = [aLength * t.siteLikes[posn], 0]\n" ,
					"        sList.append(sl)\n" ,
					"\n" ,
					"    for posn in range(aLength):\n" ,
					"        sl = ''.join(a.sequenceSlice(posn))\n" ,
					"        sDict[sl][1] += 1\n" ,
					"\n" ,
					"    t.write(), \"\\n\"\n" ,
					"    for sl in sList:\n" ,
					"        print (sl, \" %9i\" % sDict[sl][1], \" %11.2f\" % sDict[sl][0])\n" ,
	;

	close	OUTp4;
	##############################################
}

sub read_in_p4_results{
	
	my $href_hoh_found_N_of_topo_of_pattern		=	$_[0] ; # key1: topology; key2: recoded site pattern; value: N site pattern in original data	-> IN (defined)
	my $href_hoh_expected_N_of_topo_of_pattern	=	$_[1] ; # key1: topology; key2: recoded site pattern; value: expected ML N site pattern			-> IN (defined)
	my $sref_p4_outfile_name					=	$_[2] ; # name of the p4 outfile																-> IN (not changed)
	my $href_parameter_of_option				=	$_[3] ;	# key: option command, value: option value												-> IN (not changed)
	my $href_data_of_property					=	$_[4] ;	# key: property; value : property value													-> IN (not changed)
	my $sref_p4_error_message					= $_[5] ; # defined by 1 if p4 stops with an error prompt											-> IN (changed)
	
	########################################################## START EXTRACTION OF P4 RESULTS
	open IN, "<$$sref_p4_outfile_name" or die "\n\tINFILE-ERROR:Cannot read IN P4 result file ", $$sref_p4_outfile_name , "!\n" ;
	print "\n\tRead p4 result file\n";
	
	my $topology ;
	while ( my $line = <IN> ){
		#print $line;
		chomp $line ;
		
		#############################
		# Extraction of each topology with associated difference between expected and observed pattern frequencies
		if 			( $line =~ /Bad/ 								){ $$sref_p4_error_message = $line }
		elsif 	(	( $line =~ /^\(\w/ ) && ( $line !~ /\d+\);$/ )	){
			
			##############################################
			# Change unrooted quartet tree to rooted tree
			# (A,B,(C,D)) -> ((A,B),(C,D))
			(	$topology = $line )	=~ s/\s+|\n|;//g	;	#print "\tTopology:\t", $topology, "\n" ;
				$topology			=~ s/^\(/\(\(/g		;	#print "\tTopology:\t", $topology, "\n" ;
				$topology			=~ s/,\(/\),\(/g	;	#print "\tTopology:\t", $topology, "\n" ;
		}
		
		elsif	(	(	$line =~ /^\w\w\w\w\s+\d+\s+\d+/	) ||
					(	$line =~ /^\(.\w\w\w\w.,/			) ){
				
				$line			=~	s/\(|\)|,|\'//g ;
				$line			=~	s/\s+/ /g ; #print $line;
			my	@line_blocks	=	split " ", $line ;
			
			#############################
			# recode quartet pattern to M-ICE pattern code
			# best pattern match tree works currently only for nucleotide data
			# and not for RY coded sequence states
			# No error prompt implemented, yet !!!
				$line_blocks[0]		=~	tr/acgt/ACGT/ ;
			my	@site_pattern		=	split "", $line_blocks[0] ;
			my	$recoded_pattern 	=	&recode_pattern(
													
													\@site_pattern,							# list of pattern states							IN (not changed)
													\$href_parameter_of_option->{c},		# code for non-RY coded states )otherwise -> 1)		IN (not changed)
													\$href_data_of_property->{type}			# code for nucleotide sequences						IN (not changed)
										) ;
			#############################
			
			############################# (1)
			# Add observed pattern number to %hoh_found_N_of_topo_of_pattern
			$href_hoh_found_N_of_topo_of_pattern->{$topology}{$recoded_pattern} += $line_blocks[1] ;
			#############################
			
			############################# (2)
			# Add expected pattern number to %hoh_expected_N_of_topo_of_pattern
			$href_hoh_expected_N_of_topo_of_pattern->{$topology}{$recoded_pattern} += $line_blocks[2] ;
			#############################
		}
		#############################
	}
	close IN ;
	
	#############################
	# delete p4 result file
	unlink $$sref_p4_outfile_name ;
	########################################################## END EXTRACTION OF P4 RESULTS
}

sub start_phyquart{

	##########################################################
	# T11: Nap = N observed split signal - N observed signal plesiomorph (Np) - N expected convergent signal (Nk) -> each value reduced by proportion (shortest number of singelton * 4) / total number of singeltons (separately for observed and expected pattern)
	# Convergent split signal reduced of given topology reduced for each of the two others AFTER calculating the mean convergent signal value
	my	$sref_filename_in						=	$_[0]	;	#	quartet name																			In -> Defined, OUT -> Unchanged
	my	$aref_split_topologies					=	$_[1]	;	#	list of split topologies (1,2,3)														In -> Defined, OUT -> Unchanged
	my	$href_hoh_found_N_of_topo_of_pattern	=	$_[2]	;	#	key1: topology; key2: recoded site pattern; value: N site pattern in original data		In -> Defined, OUT -> Unchanged
	my	$href_hoh_expected_N_of_topo_of_pattern	=	$_[3]	;	#	key1: topology; key2: recoded site pattern; value: expected ML N site pattern			In -> Defined, OUT -> Unchanged
	my	$href_spin_value_of_quartet_topology	=	$_[4]	;	#	key1: threshold code, key2: quartet topo, value: quartet weight							In -> Undefined, OUT -> defined
	my	$sref_info_print						=	$_[5]	;	#	print code (1 -> print all) (0 - print only important output)							In -> Defined, OUT -> Unchanged
	my	$sref_quartet_taxa						=	$_[6]	;	#	string of quartet taxa for calculation print out ('t1:t2:t3:t4')						In -> Defined, OUT -> Unchanged
	my	$sref_outfile_q_calc					=	$_[7]	;	#	name single quartet calculation															In -> Defined, OUT -> Unchanged
	my	$sref_name_folder						=	$_[8]	;	#	key: subfolder code; value: subfolder name												In -> Defined, OUT -> Unchanged
	my	$sref_msa_infile_name					=	$_[9]	;	#	msa infile name																			In -> Defined, OUT -> Unchanged
	my	$href_subclade_of_taxon					=	$_[10]	;	#	key: taxon name; value: defined subclade name of taxon									In -> Defined, OUT -> Unchanged
	my	$href_sum_value_of_spin_dir				=	$_[11]	;	#	key: spin directive clan topology; value: total score after single quartet analyses		In -> Undefined, OUT -> defined
	my	$href_hoh_counter_order_top				=	$_[12]	;	#	key1: quartet tree; key2: order position; value: N occurence							IN -> Defined, OUT -> Changed
	
	my $split_topology_1	= @$aref_split_topologies[0] ;
	my $split_topology_2	= @$aref_split_topologies[1] ;
	my $split_topology_3	= @$aref_split_topologies[2] ;
	##########################################################
	
	
	
	#######################
	# open out general info file for T11 calculation print outs
	unless ( $$sref_outfile_q_calc eq 'none' ){

		my 		$info_t22_out	=	$$sref_outfile_q_calc ;
					$info_t22_out	=~	s/.txt$/_q${$sref_filename_in}.txt/ ;

		mkdir									"$$sref_name_folder/single_quartet_calculations" ;
		open 		OUTinfotxt, 	">$$sref_name_folder/single_quartet_calculations/$info_t22_out" or die "OUTFILE-ERROR: Cannot open info outfile ", $$sref_name_folder, "/single_quartet_calculations/", $info_t22_out, " in subroutine &determine_quartet_topology_11!\n";
		print 	OUTinfotxt 		"Q", $$sref_filename_in, ":\t", $$sref_quartet_taxa ,"\n"
	}
	#######################


	#############################
	# store possible pattern codes as list in:
	# set number of patterns which are not given
	# in recoded pattern string of msa
	# to zero
	for my $topology ( @$aref_split_topologies ){ for my $pattern ( qw/A B C D E F G H J K L M N I V/ ){

			unless ( $href_hoh_expected_N_of_topo_of_pattern->{$topology}{$pattern} ){ $href_hoh_expected_N_of_topo_of_pattern	->{$topology}{$pattern} = 0 }
			unless ( $href_hoh_found_N_of_topo_of_pattern		->{$topology}{$pattern} ){ $href_hoh_found_N_of_topo_of_pattern			->{$topology}{$pattern} = 0 }
		}
	}
	#############################





	#############################
	# Identification of...
	# (1)... total number of observed singleton splits
	# (2)... total number of expected singleton splits
	# (3)... lowest number of observed singleton split
	# (4)... lowest number of expected singleton split

	## (1) identic for all three topologies
	my $total_N_of_obs_singeltons	= $href_hoh_found_N_of_topo_of_pattern->{$split_topology_1}{K} + $href_hoh_found_N_of_topo_of_pattern->{$split_topology_1}{L} + $href_hoh_found_N_of_topo_of_pattern->{$split_topology_1}{M} + $href_hoh_found_N_of_topo_of_pattern->{$split_topology_1}{N} ;
	####



	## (2)
	my %total_N_of_exp_sing_of_topo ;
	$total_N_of_exp_sing_of_topo{$split_topology_1}	= $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_1}{K} + $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_1}{L} + $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_1}{M} + $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_1}{N} ;
	$total_N_of_exp_sing_of_topo{$split_topology_2}	= $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_2}{K} + $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_2}{L} + $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_2}{M} + $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_2}{N} ;
	$total_N_of_exp_sing_of_topo{$split_topology_3}	= $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_3}{K} + $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_3}{L} + $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_3}{M} + $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_3}{N} ;
	####



	## (3)
	my @singletons_observed_sorted = sort {$a<=>$b}( $href_hoh_found_N_of_topo_of_pattern->{$split_topology_1}{K}, $href_hoh_found_N_of_topo_of_pattern->{$split_topology_1}{L}, $href_hoh_found_N_of_topo_of_pattern->{$split_topology_1}{M}, $href_hoh_found_N_of_topo_of_pattern->{$split_topology_1}{N} ) ;
	#for (@singletons_observed_sorted){print "\n", $_}; exit;
	####



	## (4)
	my @singletons_expected_sorted_topo_1 = sort {$a<=>$b}( $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_1}{K}, $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_1}{L}, $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_1}{M}, $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_1}{N} ) ;
	my @singletons_expected_sorted_topo_2 = sort {$a<=>$b}( $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_2}{K}, $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_2}{L}, $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_2}{M}, $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_2}{N} ) ;
	my @singletons_expected_sorted_topo_3 = sort {$a<=>$b}( $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_3}{K}, $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_3}{L}, $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_3}{M}, $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_3}{N} ) ;
	####

	#############################



	#############################
	# Calculation of percentage of equal distributed terminal branches via shortest singleton pattern and percentage of signal reduction
	# conducted for...
	# ...(1) observed and
	# ...(2) expected pattern signal
	# --------------------
	# e.g.
	# shortest singleton number (=shortest branch proportion) : 40
	# percentage of equal distributed singleton number: 4 * 40 (4 terminal branches) = 160
	# Number of overall observed singletons: 1000
	# percentage of equal distributed signal: (160/1000) * 100 = 16 %
	# percentage of signal reduction: 100 - 16 = 84 %
	# --------------------

	## (1)
	my $perc_remaining_signal_observed ;
	if ( $total_N_of_obs_singeltons > 0 ){ $perc_remaining_signal_observed	= ( ( $singletons_observed_sorted[0] * 4 ) / $total_N_of_obs_singeltons ) * 100 } else { $perc_remaining_signal_observed = 0 }

	my $percentage_signal_red_observed	= 100 - $perc_remaining_signal_observed ;
	####



	## (2)
	my $perc_remaining_signal_exp_topo_1 ;
	if ( $total_N_of_exp_sing_of_topo{$split_topology_1} > 0 )	{ $perc_remaining_signal_exp_topo_1 = ( ( $singletons_expected_sorted_topo_1[0] * 4 ) / $total_N_of_exp_sing_of_topo{$split_topology_1} ) * 100 }
	else 																												{ $perc_remaining_signal_exp_topo_1 = 0 }

	my $perc_remaining_signal_exp_topo_2 ;
	if ( $total_N_of_exp_sing_of_topo{$split_topology_2} > 0 )	{ $perc_remaining_signal_exp_topo_2 = ( ( $singletons_expected_sorted_topo_2[0] * 4 ) / $total_N_of_exp_sing_of_topo{$split_topology_2} ) * 100 }
	else																												{ $perc_remaining_signal_exp_topo_2 = 0 }

	my $perc_remaining_signal_exp_topo_3 ;
	if ( $total_N_of_exp_sing_of_topo{$split_topology_3} > 0 )	{ $perc_remaining_signal_exp_topo_3 = ( ( $singletons_expected_sorted_topo_3[0] * 4 ) / $total_N_of_exp_sing_of_topo{$split_topology_3} ) * 100 }
	else																												{ $perc_remaining_signal_exp_topo_3 = 0 }

	my $percentage_signal_red_exp_topo1		= 100 - $perc_remaining_signal_exp_topo_1 ;
	my $percentage_signal_red_exp_topo2		= 100 - $perc_remaining_signal_exp_topo_2 ;
	my $percentage_signal_red_exp_topo3		= 100 - $perc_remaining_signal_exp_topo_3 ;
	####

	# info print
	my @topologies = ( $split_topology_1, $split_topology_2, $split_topology_3 ) ;
	unless ( $$sref_outfile_q_calc eq 'none' ){

		my @taxa_quart	=	split ":", $$sref_quartet_taxa;

		my $topo_counter = 1 ;
		for my $topo ( @topologies ){

			if ( $topo_counter == 1 ){

				print OUTinfotxt "\nN observed singletons (sorted) for all 3 topologies:" ;
				for ( @singletons_observed_sorted ){ print OUTinfotxt "\t", $_ }
				print OUTinfotxt	"\nN observed K (",$taxa_quart[0],"): ", $href_hoh_found_N_of_topo_of_pattern->{$split_topology_1}{K} ,
													"\nN observed L (",$taxa_quart[1],"): ", $href_hoh_found_N_of_topo_of_pattern->{$split_topology_1}{L} ,
													"\nN observed M (",$taxa_quart[2],"): ", $href_hoh_found_N_of_topo_of_pattern->{$split_topology_1}{M} ,
													"\nN observed N (",$taxa_quart[3],"): ", $href_hoh_found_N_of_topo_of_pattern->{$split_topology_1}{N} ,
													"\n\nN total observed singletons for all 3 topologies:\t", $total_N_of_obs_singeltons, "\n",
													"Percentage remaining signal (observed) for all 3 topologies:\t"	, $perc_remaining_signal_observed, "\n",
													"Percentage reduced signal (observed) for all 3 topologies:\t"		, $percentage_signal_red_observed, "\n\n";
			}

			print OUTinfotxt "N expected singletons (sorted) for topology ", $topo_counter, ":" ;
			if 			( $topo_counter == 1 ){ for ( @singletons_expected_sorted_topo_1 ){ print OUTinfotxt "\t", $_ } }
			elsif 	( $topo_counter == 2 ){ for ( @singletons_expected_sorted_topo_2 ){ print OUTinfotxt "\t", $_ } }
			elsif 	( $topo_counter == 3 ){ for ( @singletons_expected_sorted_topo_3 ){ print OUTinfotxt "\t", $_ } }
			print OUTinfotxt	"\nN expected K (",$taxa_quart[0],"): ", $href_hoh_expected_N_of_topo_of_pattern->{$topo}{K};
			print OUTinfotxt	"\nN expected L (",$taxa_quart[1],"): ", $href_hoh_expected_N_of_topo_of_pattern->{$topo}{L};
			print OUTinfotxt	"\nN expected M (",$taxa_quart[2],"): ", $href_hoh_expected_N_of_topo_of_pattern->{$topo}{M};
			print OUTinfotxt	"\nN expected N (",$taxa_quart[3],"): ", $href_hoh_expected_N_of_topo_of_pattern->{$topo}{N};

			print OUTinfotxt "\nN total expected singletons (sorted) for topology ", $topo_counter, ":\t", $total_N_of_exp_sing_of_topo{$topo} ;

			print OUTinfotxt "\nPercentage remaining signal (expected):\t";
			if 			( $topo_counter == 1 ){ print OUTinfotxt $perc_remaining_signal_exp_topo_1, "%" }
			elsif 	( $topo_counter == 2 ){ print OUTinfotxt $perc_remaining_signal_exp_topo_2, "%" }
			elsif 	( $topo_counter == 3 ){ print OUTinfotxt $perc_remaining_signal_exp_topo_3, "%" }

			print OUTinfotxt "\nPercentage reduced signal (expected):\t";
			if 			( $topo_counter == 1 ){ print OUTinfotxt $percentage_signal_red_exp_topo1, "%" }
			elsif 	( $topo_counter == 2 ){ print OUTinfotxt $percentage_signal_red_exp_topo2, "%" }
			elsif 	( $topo_counter == 3 ){ print OUTinfotxt $percentage_signal_red_exp_topo3, "%" }

			print OUTinfotxt "\n\n" ;

			$topo_counter++
		}
	}
	#############################



	#############################
	# Calculation of oberved, overall split supporting pattern number for each of the three possible topologies (nb)
	# Split supporting symmetric and asymmetric pattern topology 1 ((T1,T2),(T3,T4)) -> A + D + E
	# Split supporting symmetric and asymmetric pattern topology 2 ((T1,T3),(T2,T4)) -> B + F + G
	# Split supporting symmetric and asymmetric pattern topology 3 ((T1,T4),(T2,T3)) -> A + H + J
	my $spin_topology_1_nb	= 0 + $href_hoh_found_N_of_topo_of_pattern->{$split_topology_1}{A} + $href_hoh_found_N_of_topo_of_pattern->{$split_topology_1}{D} + $href_hoh_found_N_of_topo_of_pattern->{$split_topology_1}{E} ;
	my $spin_topology_2_nb	= 0 + $href_hoh_found_N_of_topo_of_pattern->{$split_topology_2}{B} + $href_hoh_found_N_of_topo_of_pattern->{$split_topology_2}{F} + $href_hoh_found_N_of_topo_of_pattern->{$split_topology_2}{G} ;
	my $spin_topology_3_nb	= 0 + $href_hoh_found_N_of_topo_of_pattern->{$split_topology_3}{C} + $href_hoh_found_N_of_topo_of_pattern->{$split_topology_3}{H} + $href_hoh_found_N_of_topo_of_pattern->{$split_topology_3}{J} ;

	unless ( $$sref_outfile_q_calc eq 'none' ){

		print OUTinfotxt	"Observed Split Signal (Nb)",
							"\nN observed pattern ", $split_topology_1, ":\t", $spin_topology_1_nb, "\tSymetric (A):", $href_hoh_found_N_of_topo_of_pattern->{$split_topology_1}{A}, "\tAsymetric (D):", $href_hoh_found_N_of_topo_of_pattern->{$split_topology_1}{D}, "\tAsymetric (E):", $href_hoh_found_N_of_topo_of_pattern->{$split_topology_1}{E},
							"\nN observed pattern ", $split_topology_2, ":\t", $spin_topology_2_nb, "\tSymetric (B):", $href_hoh_found_N_of_topo_of_pattern->{$split_topology_1}{B}, "\tAsymetric (F):", $href_hoh_found_N_of_topo_of_pattern->{$split_topology_1}{F}, "\tAsymetric (G):", $href_hoh_found_N_of_topo_of_pattern->{$split_topology_1}{G},
							"\nN observed pattern ", $split_topology_3, ":\t", $spin_topology_3_nb, "\tSymetric (C):", $href_hoh_found_N_of_topo_of_pattern->{$split_topology_1}{C}, "\tAsymetric (H):", $href_hoh_found_N_of_topo_of_pattern->{$split_topology_1}{H}, "\tAsymetric (J):", $href_hoh_found_N_of_topo_of_pattern->{$split_topology_1}{J}, "\n"
	}
	#############################



	#############################
	# Reduction of observed signal nb_red in proportion of reduced singletons
	# --------------------
	# e.g.
	# Observed signal pattern topology1 (asymetric and symetric split pattern): 1000
	# Reduction proportion ($percentage_signal_red_observed): 58 %
	# Remaining split signal: 1000 - (( 1000/100 ) * 58 ) = 420 = 42% remaining pattern = $perc_remaining_signal_observed
	my $spin_topology_1_nb_red = $spin_topology_1_nb - (( $spin_topology_1_nb / 100 ) * $percentage_signal_red_observed ) ;
	my $spin_topology_2_nb_red = $spin_topology_2_nb - (( $spin_topology_2_nb / 100 ) * $percentage_signal_red_observed ) ;
	my $spin_topology_3_nb_red = $spin_topology_3_nb - (( $spin_topology_3_nb / 100 ) * $percentage_signal_red_observed ) ;


	unless ( $$sref_outfile_q_calc eq 'none' ){

		print	OUTinfotxt	"\nPercentage remaining signal (observed) after reduction:\n",
							"topology 1 ", $split_topology_1, "\t", $spin_topology_1_nb_red, "\n",
							"topology 2 ", $split_topology_2, "\t", $spin_topology_2_nb_red, "\n",
							"topology 3 ", $split_topology_3, "\t", $spin_topology_3_nb_red, "\n\n"
	}
	#############################



	#############################
	# Assignment of observed plesiomorphic split pattern number in accordance to spin direction (only asymmetric splits) (np)
	# Each topology can be divided into two splin directions
	#--------------------------------------------------------------------------------------
	# e.g.
	# left to right (lr):
	# topology 1 ((T1,T2),(T3,T4)): (T1,T2) -> (T3,T4) , plesiomorphic pattern XX|YZ -> D
	# right to left (rl):
	# topology 1 ((T1,T2),(T3,T4)): (T3,T4) <- (T1,T2) , plesiomorphic pattern YZ|XX -> E
	#--------------------------------------------------------------------------------------
	my $spin_topology_1_lr_np	= 0 + $href_hoh_found_N_of_topo_of_pattern->{$split_topology_1}{D} ;
	my $spin_topology_1_rl_np	= 0 + $href_hoh_found_N_of_topo_of_pattern->{$split_topology_1}{E} ;
	my $spin_topology_2_lr_np	= 0 + $href_hoh_found_N_of_topo_of_pattern->{$split_topology_2}{G} ;
	my $spin_topology_2_rl_np	= 0 + $href_hoh_found_N_of_topo_of_pattern->{$split_topology_2}{F} ;
	my $spin_topology_3_lr_np	= 0 + $href_hoh_found_N_of_topo_of_pattern->{$split_topology_3}{J} ;
	my $spin_topology_3_rl_np	= 0 + $href_hoh_found_N_of_topo_of_pattern->{$split_topology_3}{H} ;

	unless ( $$sref_outfile_q_calc eq 'none' ){

		print	OUTinfotxt	"Plesiomorphic Split Signal (Np)",
							"\nN observed plesiomorphic states ", $split_topology_1, "lr:\t", $spin_topology_1_lr_np,
							"\nN observed plesiomorphic states ", $split_topology_1, "rl:\t", $spin_topology_1_rl_np,
							"\nN observed plesiomorphic states ", $split_topology_2, "lr:\t", $spin_topology_2_lr_np,
							"\nN observed plesiomorphic states ", $split_topology_2, "rl:\t", $spin_topology_2_rl_np,
							"\nN observed plesiomorphic states ", $split_topology_3, "lr:\t", $spin_topology_3_lr_np,
							"\nN observed plesiomorphic states ", $split_topology_3, "rl:\t", $spin_topology_3_rl_np, "\n"
	}
	#############################



	#############################
	# Reduction of observed plesiomorphic split pattern number in accordance to spin direction (only asymmetric splits) (np_red)
	my $spin_topology_1_lr_np_red	= $spin_topology_1_lr_np - (( $spin_topology_1_lr_np / 100 ) * $percentage_signal_red_observed ) ;
	my $spin_topology_1_rl_np_red	= $spin_topology_1_rl_np - (( $spin_topology_1_rl_np / 100 ) * $percentage_signal_red_observed ) ;
	my $spin_topology_2_lr_np_red	= $spin_topology_2_lr_np - (( $spin_topology_2_lr_np / 100 ) * $percentage_signal_red_observed ) ;
	my $spin_topology_2_rl_np_red	= $spin_topology_2_rl_np - (( $spin_topology_2_rl_np / 100 ) * $percentage_signal_red_observed ) ;
	my $spin_topology_3_lr_np_red	= $spin_topology_3_lr_np - (( $spin_topology_3_lr_np / 100 ) * $percentage_signal_red_observed ) ;
	my $spin_topology_3_rl_np_red	= $spin_topology_3_rl_np - (( $spin_topology_3_rl_np / 100 ) * $percentage_signal_red_observed ) ;


	unless ( $$sref_outfile_q_calc eq 'none' ){

		print	OUTinfotxt	"\nN observed plesiomorphic states ", $split_topology_1, "lr reduced:\t", $spin_topology_1_lr_np_red,
							"\nN observed plesiomorphic states ", $split_topology_1, "rl reduced:\t", $spin_topology_1_rl_np_red,
							"\nN observed plesiomorphic states ", $split_topology_2, "lr reduced:\t", $spin_topology_2_lr_np_red,
							"\nN observed plesiomorphic states ", $split_topology_2, "rl reduced:\t", $spin_topology_2_rl_np_red,
							"\nN observed plesiomorphic states ", $split_topology_3, "lr reduced:\t", $spin_topology_3_lr_np_red,
							"\nN observed plesiomorphic states ", $split_topology_3, "rl reduced:\t", $spin_topology_3_rl_np_red, "\n"
	}
	#############################



	#############################
	# Calculate expected number of convergent split pattern for each of the three possible quartet topologies (nk)
	#--------------------------------------------------------------------------------------
	# e.g. split supporting pattern of topology 1 ((T1,T2),(T3,T4)):
	# A, D, E
	# convergent split pattern which are supporting topology 1 in both other quartets:
	# expected number of A, D, E in topology 2 and topology 3
	# spin direction has to be considered for D and E
	# ((T1,T2),(T3,T4)) lr : convergent signal on the right side of the internal branch (XY:ZZ) -> E
	# ((T1,T2),(T3,T4)) rl : convergent signal on the left  side of the internal branch (XX:YZ) -> D
	# symmetric split numbers supporting topology 1 in other quartets are spin independent
	#--------------------------------------------------------------------------------------
	my $spin_topology_1_lr_nk = ( 0 + $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_2}{A} + $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_3}{A} + $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_2}{E} + $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_3}{E} ) ;
	my $spin_topology_1_rl_nk = ( 0 + $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_2}{A} + $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_3}{A} + $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_2}{D} + $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_3}{D} ) ;
	my $spin_topology_2_lr_nk = ( 0 + $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_1}{B} + $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_3}{B} + $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_1}{F} + $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_3}{F} ) ;
	my $spin_topology_2_rl_nk = ( 0 + $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_1}{B} + $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_3}{B} + $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_1}{G} + $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_3}{G} ) ;
	my $spin_topology_3_lr_nk = ( 0 + $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_1}{C} + $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_2}{C} + $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_1}{H} + $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_2}{H} ) ;
	my $spin_topology_3_rl_nk = ( 0 + $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_1}{C} + $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_2}{C} + $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_1}{J} + $href_hoh_expected_N_of_topo_of_pattern->{$split_topology_2}{J} ) ;

	unless ( $$sref_outfile_q_calc eq 'none' ){

		print	OUTinfotxt	"\nConvergent Split Signal (Nk)",
							"\nN expected convergent split pattern (ML)", $split_topology_1, "lr:\t", $spin_topology_1_lr_nk,
							"\nN expected convergent split pattern (ML)", $split_topology_1, "rl:\t", $spin_topology_1_rl_nk,
							"\nN expected convergent split pattern (ML)", $split_topology_2, "lr:\t", $spin_topology_2_lr_nk,
							"\nN expected convergent split pattern (ML)", $split_topology_2, "rl:\t", $spin_topology_2_rl_nk,
							"\nN expected convergent split pattern (ML)", $split_topology_3, "lr:\t", $spin_topology_3_lr_nk,
							"\nN expected convergent split pattern (ML)", $split_topology_3, "rl:\t", $spin_topology_3_rl_nk, "\n"
	}
	#############################



	#############################
	# Calculate...
	# ...(1) mean reduction proportion for each of the three topologies (inferred for each topology from the remaining two topology singletons
	# ...(2) the expected number of reduced convergent split pattern for each of the three possible quartet topologies (nk_red)
	# ...(3) the mean of convergent split signal

	## (1)
	my $mean_red_prop_for_exp_convergencies_top_1 = ( $percentage_signal_red_exp_topo2 + $percentage_signal_red_exp_topo3 ) / 2 ;
	my $mean_red_prop_for_exp_convergencies_top_2 = ( $percentage_signal_red_exp_topo1 + $percentage_signal_red_exp_topo3 ) / 2 ;
	my $mean_red_prop_for_exp_convergencies_top_3 = ( $percentage_signal_red_exp_topo1 + $percentage_signal_red_exp_topo2 ) / 2 ;
	####

	## (2)
	my $spin_topology_1_lr_nk_red = $spin_topology_1_lr_nk - (( $spin_topology_1_lr_nk / 100 ) * $mean_red_prop_for_exp_convergencies_top_1 ) ;
	my $spin_topology_1_rl_nk_red = $spin_topology_1_rl_nk - (( $spin_topology_1_rl_nk / 100 ) * $mean_red_prop_for_exp_convergencies_top_1 ) ;
	my $spin_topology_2_lr_nk_red = $spin_topology_2_lr_nk - (( $spin_topology_2_lr_nk / 100 ) * $mean_red_prop_for_exp_convergencies_top_2 ) ;
	my $spin_topology_2_rl_nk_red = $spin_topology_2_rl_nk - (( $spin_topology_2_rl_nk / 100 ) * $mean_red_prop_for_exp_convergencies_top_2 ) ;
	my $spin_topology_3_lr_nk_red = $spin_topology_3_lr_nk - (( $spin_topology_3_lr_nk / 100 ) * $mean_red_prop_for_exp_convergencies_top_3 ) ;
	my $spin_topology_3_rl_nk_red = $spin_topology_3_rl_nk - (( $spin_topology_3_rl_nk / 100 ) * $mean_red_prop_for_exp_convergencies_top_3 ) ;
	####

	## (3)
	my $spin_topology_1_lr_nk_red_mean = $spin_topology_1_lr_nk_red / 2 ;
	my $spin_topology_1_rl_nk_red_mean = $spin_topology_1_rl_nk_red / 2 ;
	my $spin_topology_2_lr_nk_red_mean = $spin_topology_2_lr_nk_red / 2 ;
	my $spin_topology_2_rl_nk_red_mean = $spin_topology_2_rl_nk_red / 2 ;
	my $spin_topology_3_lr_nk_red_mean = $spin_topology_3_lr_nk_red / 2 ;
	my $spin_topology_3_rl_nk_red_mean = $spin_topology_3_rl_nk_red / 2 ;
	####

	unless ( $$sref_outfile_q_calc eq 'none' ){

		print	OUTinfotxt	"\nMean reduction proportion for convergent signal of topology 1 "	, $split_topology_1, "\t", $mean_red_prop_for_exp_convergencies_top_1, "\t", (($percentage_signal_red_exp_topo2 + $percentage_signal_red_exp_topo3) / 2),
							"\nMean reduction proportion for convergent signal of topology 2 "	, $split_topology_2, "\t", $mean_red_prop_for_exp_convergencies_top_2, "\t", (($percentage_signal_red_exp_topo1 + $percentage_signal_red_exp_topo3) / 2),
							"\nMean reduction proportion for convergent signal of topology 3 "	, $split_topology_3, "\t", $mean_red_prop_for_exp_convergencies_top_3, "\t", (($percentage_signal_red_exp_topo1 + $percentage_signal_red_exp_topo2) / 2),
							"\n",
							"\nN expected convergent split pattern (ML) "						, $split_topology_1, "lr reduced:\t", $spin_topology_1_lr_nk_red,
							"\nN expected convergent split pattern (ML) "						, $split_topology_1, "rl reduced:\t", $spin_topology_1_rl_nk_red,
							"\nN expected convergent split pattern (ML) "						, $split_topology_2, "lr reduced:\t", $spin_topology_2_lr_nk_red,
							"\nN expected convergent split pattern (ML) "						, $split_topology_2, "rl reduced:\t", $spin_topology_2_rl_nk_red,
							"\nN expected convergent split pattern (ML) "						, $split_topology_3, "lr reduced:\t", $spin_topology_3_lr_nk_red,
							"\nN expected convergent split pattern (ML) "						, $split_topology_3, "rl reduced:\t", $spin_topology_3_rl_nk_red,
							"\n",
							"\nMean N expected convergent split pattern (ML) "					, $split_topology_1, "lr reduced:\t", $spin_topology_1_lr_nk_red_mean,
							"\nMean N expected convergent split pattern (ML) "					, $split_topology_1, "rl reduced:\t", $spin_topology_1_rl_nk_red_mean,
							"\nMean N expected convergent split pattern (ML) "					, $split_topology_2, "lr reduced:\t", $spin_topology_2_lr_nk_red_mean,
							"\nMean N expected convergent split pattern (ML) "					, $split_topology_2, "rl reduced:\t", $spin_topology_2_rl_nk_red_mean,
							"\nMean N expected convergent split pattern (ML) "					, $split_topology_3, "lr reduced:\t", $spin_topology_3_lr_nk_red_mean,
							"\nMean N expected convergent split pattern (ML) "					, $split_topology_3, "rl reduced:\t", $spin_topology_3_rl_nk_red_mean, "\n"
	}
	#############################



	#############################
	# Calculate potentially correct (apomorph) split signal for each spin direction of a given topology (nap)
	# nap = nb (observed) - np (expected) - nk (expected)
	my %spin_value_of_topology ;
	$spin_value_of_topology{$split_topology_1."::lr"}	= $spin_topology_1_nb_red - $spin_topology_1_lr_np_red - $spin_topology_1_lr_nk_red_mean ;
	$spin_value_of_topology{$split_topology_1."::rl"}	= $spin_topology_1_nb_red - $spin_topology_1_rl_np_red - $spin_topology_1_rl_nk_red_mean ;
	$spin_value_of_topology{$split_topology_2."::lr"}	= $spin_topology_2_nb_red - $spin_topology_2_lr_np_red - $spin_topology_2_lr_nk_red_mean ;
	$spin_value_of_topology{$split_topology_2."::rl"}	= $spin_topology_2_nb_red - $spin_topology_2_rl_np_red - $spin_topology_2_rl_nk_red_mean ;
	$spin_value_of_topology{$split_topology_3."::lr"}	= $spin_topology_3_nb_red - $spin_topology_3_lr_np_red - $spin_topology_3_lr_nk_red_mean ;
	$spin_value_of_topology{$split_topology_3."::rl"}	= $spin_topology_3_nb_red - $spin_topology_3_rl_np_red - $spin_topology_3_rl_nk_red_mean ;
	#############################



	#############################
	# sort topologies by number of potentially apomorph signal (nap)
	my @sorted_topologies			= sort {$spin_value_of_topology{$b}			<=> $spin_value_of_topology{$a}}			keys %spin_value_of_topology ;

	unless ( $$sref_outfile_q_calc eq 'none' ){

		print OUTinfotxt "\nPotential apomorph split signal of each spin dependend taxon topology:" ;
		for (@sorted_topologies)		{ print OUTinfotxt "\n", $_,  "\t", $spin_value_of_topology{$_} }

		if ( $spin_value_of_topology{$sorted_topologies[5]} < 0 ){

			print OUTinfotxt	"\n\nSpin value of lowest supported quartet < 0:\t", $spin_value_of_topology{$sorted_topologies[5]},
												"\nNormalized spin values (+", abs ($spin_value_of_topology{$sorted_topologies[5]}), "):\n",
												"\nPotentially apomorph split signal of each spin dependend clan topology after normalization:\n"
		}

		elsif ( $spin_value_of_topology{$sorted_topologies[5]} > 0 ){

			print OUTinfotxt	"\n\nSpin value of lowest supported quartet > 0:\t", $spin_value_of_topology{$sorted_topologies[5]},
												"\nNormalized spin values (-", abs ($spin_value_of_topology{$sorted_topologies[5]}), "):\n",
												"\nPotentially apomorph split signal of each spin dependend clan topology after normalization:\n"
		}
	}
	#############################



	#############################
	# (1) Assign taxonnames to defined subclades according to each single quartet topology and spin direction
	# e.g.:
	# $top_q = (taxA,taxB),(taxC,taxD) -> taxA,taxB,taxC,taxD
	# $topology = (subclade1,subclade2),(subclade3,subclade4)
	#
	# (2) Normalisation of single spin values due to the lowest spin value (set to zero)
	# Sum up of single spin dependend spin scores
	#
	# (3) Count number of order position of each quartet topology (positions possible)
	my $counter_position = 1 ;
	for my $clan_tree	( @sorted_topologies ){

		#### (1)
		(my $clan_tree_spin	= $clan_tree)	=~	s/(\()+|(\))+//g ;

		my	@tree_spin_parts	=	split "::", $clan_tree_spin ;
		my	@taxa							=	split ",", $tree_spin_parts[0] ;
		my	@sister_le				=	sort ( $href_subclade_of_taxon->{$taxa[0]}, $href_subclade_of_taxon->{$taxa[1]} ) ;
		my	@sister_ri				=	sort ( $href_subclade_of_taxon->{$taxa[2]}, $href_subclade_of_taxon->{$taxa[3]} ) ;
		my	$sister_le				=	join ",", @sister_le ;	#print "\n", $sister_le, "\n";
		my	$sister_ri				=	join ",", @sister_ri ;	#print "\n", $sister_ri, "\n"; exit;

		my	$clan_topology_spin ;
		if		( $tree_spin_parts[1] eq "lr" ){ $clan_topology_spin	=	"(".$sister_le."),(".$sister_ri.")" }
		elsif	( $tree_spin_parts[1] eq "rl" ){ $clan_topology_spin	=	"(".$sister_ri."),(".$sister_le.")" }
		else	{ die "\n\tBUG-ERROR: Cannot assign spin \"", $tree_spin_parts[1], "\" in PhyQuart subroutine!\n\tPlease, report BUG to system developer!\n\t" }
		####

		#### (2)
		my	$norm_spin_value ;
		if		( $spin_value_of_topology{$sorted_topologies[5]} < 	0	){ $norm_spin_value = ( $spin_value_of_topology{$clan_tree} + abs ($spin_value_of_topology{$sorted_topologies[5]}) ) ; $href_sum_value_of_spin_dir->{$clan_topology_spin} += $norm_spin_value }
		elsif	( $spin_value_of_topology{$sorted_topologies[5]} >= 0	){ $norm_spin_value = ( $spin_value_of_topology{$clan_tree} - abs ($spin_value_of_topology{$sorted_topologies[5]}) ) ; $href_sum_value_of_spin_dir->{$clan_topology_spin} += $norm_spin_value }
		else	{ die "\n\tBUG-ERROR: Cannot assign split values in PhyQuart subroutine!\n\tPlease, report BUG to system developer!\n\t" }

		push @{$href_spin_value_of_quartet_topology->{$clan_topology_spin}}, $norm_spin_value ;
		####

		#### (3)
		$href_hoh_counter_order_top->{$clan_topology_spin}{$counter_position}++ ;
		$counter_position++ ;
		####

		unless ( $$sref_outfile_q_calc eq 'none' ){

			print OUTinfotxt	$clan_tree, "\t", $norm_spin_value, "\n"
		}
	}

	close OUTinfotxt ;
}

sub median{

	my	$aref_vals	=	$_[0];
		@$aref_vals	=	sort {$a <=> $b} @$aref_vals	;
	my	$len		=	@$aref_vals	;

	#odd?
	if	( $len%2 ){ return $aref_vals->[int($len/2)] }

	#even
	else{ return ($aref_vals->[int($len/2)-1] + $aref_vals->[int($len/2)])/2 }
}

sub print_svg_triangle_graphic{

	my	$aref_lol_data					=	$_[0]	;	#	list1: data run (quartet) number; list 2: quartet topology (0,1,2) ; value: split support found for given quartet and topology	In -> Defined, OUT -> Unchanged
	my	$sref_outname					=	$_[1]	;	#	output filename																													In -> Defined, OUT -> Unchanged
	my	$aref_sorted_top				=	$_[2]	;	#	list of sorted clade topologies																									In -> Defined, OUT -> Unchanged
	my	$href_name_of_repeat			=	$_[3]	;	#	key: quartet number (0 to @quartets-1); value: quartet taxa (tax1:tax2:tax3:tax4)												In -> Defined, OUT -> Unchanged
	my	$sref_graphic_title				=	$_[4]	;	#	title of svg graphic																											In -> Defined, OUT -> Unchanged
	my	$href_mean_of_top				=	$_[5]	;	#	key: toponumber (0, 1, or 2); value: mean split score																			In -> Defined, OUT -> Unchanged
	my	$href_medi_of_top				=	$_[6]	;	#	key: toponumber (0, 1, or 2); value: median split score																			In -> Defined, OUT -> Unchanged
	my	$href_graphic_setup				=	$_[7]	;	#	key: svgoption; value: option setup																								In -> Defined, OUT -> Unchanged
	my	$sref_name_folder				=	$_[8]	;	#	key: subfolder code; value: subfolder name																						In -> Defined, OUT -> Unchanged
	my	$sref_median_print_code			=	$_[9]	;	#	code for median score print out																									In -> Defined, OUT -> Unchanged
	my	$sref_mean_median_legend		=	$_[10]	;	#	code for triangle overall mean or median print out																				In -> Defined, OUT -> Unchanged
	my	$sref_triangle_color			=	$_[11]	;	#	triangle background color																										In -> Defined, OUT -> Unchanged
	my	$sref_svg_background			=	$_[12]	;	#	background color svg																											In -> Defined, OUT -> Unchanged
	my	$sref_header_color				=	$_[13]	;	#	header color																													In -> Defined, OUT -> Unchanged

	#################################################################### START Triangle Coordinates Calculation

	#################
	# Calculate single x y coordinates of given triplet data
	# Single triplet data values must sum each to 1
	my	@lol_x_y ;					# list1: data run number; list 2: x,y coordinates for triangle output of triplet run
	my	%hol_values_of_triposition;	# key: triplet position (0,1,or 2); value: list of observed triplet scores

	DATA:
	for (0 .. @$aref_lol_data-1){

		( $lol_x_y[$_][0], $lol_x_y[$_][1] )	=	&calculate_xy	( $aref_lol_data->[$_][1], $aref_lol_data->[$_][2] ) ;

		push	@{$hol_values_of_triposition{0}},	$aref_lol_data->[$_][0] ;
		push	@{$hol_values_of_triposition{1}},	$aref_lol_data->[$_][1] ;
		push	@{$hol_values_of_triposition{2}},	$aref_lol_data->[$_][2] ;
	}
	#################

	#################
	# store mean and median values of single triplet positions
	# calculate x,y coordinates for mean output given a triangle
	my @mean_coordinates ;
	( $mean_coordinates[0], $mean_coordinates[1] )		= &calculate_xy	( $href_mean_of_top->{1}, $href_mean_of_top->{2} ) ;

	my @median_coordinates ;
	( $median_coordinates[0], $median_coordinates[1] )	= &calculate_xy	( $href_medi_of_top->{1}, $href_medi_of_top->{2} ) ;
	#################

	#################
	# Subroutine Calculation
	sub calculate_xy{

		my $b	=	$_[0] ;
		my $c	=	$_[1] ;

		my $x	=	0.5 * ( 2 * $b + $c ) ;
		my $y	=	( sqrt(3) / 2 ) * $c ;

		my @x_y	= ( $x, $y ) ;

		return @x_y ;
	}
	#################

	#################################################################### END Triangle Coordinates Calculation



	#################################################################### Triangle SVG Print OUT ####################################################################

	###################################### START SVG Print

	#################
	# Standard values of grafik output considering subtree circles, terminal branches, and internal node circle radius
	my	$stroke_width_branches				=	$href_graphic_setup->{stroke_width_branches}	;
	my	$stroke_color_branches				=	$href_graphic_setup->{stroke_color_branches}	;
	my	$axis_stroke_opacity				=	$href_graphic_setup->{axis_stroke_opacity}		;
	my	$axis_fill_opacity					=	$href_graphic_setup->{axis_fill_opacity}		;
	my	$axis_fill_color					=	$href_graphic_setup->{axis_fill_color}			;
	my	$text_cursor						=	"evt.target.ownerDocument.getElementById('spktText0').firstChild.nodeValue="	;
	my	$text_color_cursor					=	"evt.target.setAttributeNS(null, 'fill', '"		;
	my	$text_curser_fill					=	$href_graphic_setup->{text_curser_fill}			;
	my	$text_font_size						=	$href_graphic_setup->{text_font_size}			;
	my	$text_stroke						=	$href_graphic_setup->{text_stroke}				;
	my	$text_stroke_2						=	$href_graphic_setup->{text_stroke_2}			;
	my	$text_font_weight					=	$href_graphic_setup->{text_font_weight}			;
	my	$text_stroke_width					=	$href_graphic_setup->{text_stroke_width}		;
	my	$text_fill							=	$href_graphic_setup->{text_fill}				;
	my	$text_fill_2						=	$href_graphic_setup->{text_fill_2}				;
	my	$text_font_family					=	$href_graphic_setup->{text_font_family}			;
	my	$scale_color						=	$href_graphic_setup->{scale_color}				;
	my	$scale_stroke_width					=	$href_graphic_setup->{scale_stroke_width}		;
	my	$scale_std_value					=	$href_graphic_setup->{scale_std_value}			;
	my	$scale_font_size					=	$href_graphic_setup->{scale_font_size}			;
	my	$legend_font_size					=	$href_graphic_setup->{legend_font_size}			;
	my	$legend_font_weight					=	$href_graphic_setup->{legend_font_weight}		;
	my	$legend_stroke_width				=	$href_graphic_setup->{legend_stroke_width}		;
	my	$circle_radius						=	$href_graphic_setup->{circle_radius}			;
	my	$circle_radius_mean					=	$href_graphic_setup->{circle_radius_mean}		;
	my	$circle_stroke_color				=	$href_graphic_setup->{circle_stroke_color}		;
	my	$circle_stroke_color_mean			=	$href_graphic_setup->{circle_stroke_color_mean}	;
	my	$circle_stroke_width				=	$href_graphic_setup->{circle_stroke_width}		;
	my	$circle_fill_color					=	$href_graphic_setup->{circle_fill_color}		;
	my	$circle_fill_color_mean				=	$href_graphic_setup->{circle_fill_color_mean}	;
	my	$circle_curser_fill					=	$href_graphic_setup->{circle_curser_fill}		;
	my	$split_factor						=	$href_graphic_setup->{split_factor}				;
	my	$font_size							=	$href_graphic_setup->{font_size}				;
	################

	#################
	# Definement of docu dimensions
	my	$zoom_factor	=	150 ;

	my	$docu_width		=	2 * $zoom_factor ;
	my	$docu_height	=	2 * $zoom_factor ;
	my	$y_point		=	1 * $zoom_factor ;

	my	@x_y_corner1	=	( 0.5 * $zoom_factor					, $y_point + 0.5 * $zoom_factor ) ;
	my	@x_y_corner2	=	( $docu_width - ( 0.5 * $zoom_factor )	, $y_point + 0.5 * $zoom_factor ) ;
	my	@x_y_corner3	=	( 0.5 * $docu_width						, $x_y_corner1[1] - ( sqrt(3) / 2 ) * $zoom_factor ) ;
	#################

	#################
	# open OUT svg output
	open	OUTsvg, ">$$sref_name_folder/$$sref_outname" or die "\n\tOUTFILE-ERROR: Cannot open triangle graphic ", $$sref_name_folder, "/", $$sref_outname, " in subroutine &ternary!\n" ;
	#################




	#################
	# SVG Header
	print	OUTsvg	"<?xml version=\"1.0\"?>\n",
					"<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n",
					"\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n\n",
					"<svg width=\"", $docu_width*2, "\" height=\"", $docu_height, "\" fill=\"", $$sref_header_color ,"\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">\n\n",
					"\t<defs>\n",
					"\n",
					"\t</defs>\n",
					"\n";
	#################



	################################################### START triangle
	print OUTsvg	"\n\t<rect width=\"", $docu_width*2, "\" height=\"", $docu_height, "\" fill=\"", $$sref_svg_background, "\" />" ,
					"\n\t<polygon points=\""	, $x_y_corner1[0]	, ","	, $x_y_corner1[1]	," "	, $x_y_corner2[0]	, ","	, $x_y_corner2[1]	, " "	, $x_y_corner3[0]	, ","	, $x_y_corner3[1]	,
					"\" style=\"stroke-width:"	, $stroke_width_branches, ";stroke:"		, $stroke_color_branches, ";fill:", $$sref_triangle_color, "\" />" ;
	################# END triangle

	################# START Corner Legend
	# print lower left corner code
	print OUTsvg	"\n\t<text x=\""	, $x_y_corner1[0]-0.125*$zoom_factor		,"\" y=\""						, $x_y_corner1[1]+$text_font_size*0.5	,
					"\" style=\"font-size:"		, $font_size, ";font-family: ", $text_font_family, "; stroke-width:", $text_stroke_width	,";\" ",

					"onmouseover=\""	, $text_cursor ," ' ", $aref_sorted_top->[0] ," '; ", $text_color_cursor, $text_curser_fill ,"');\" ",
					"onmouseout=\""		, $text_cursor ," ' ';", $text_color_cursor, $text_fill ,"');\" fill=\"", $text_fill ,"\" > Q1 </text>" ;

	# print lower right corner code
	print OUTsvg	"\n\t<text x=\""	, $x_y_corner2[0]+0.05*$zoom_factor	,"\" y=\""						, $x_y_corner1[1]+$text_font_size*0.5	,
					"\" style=\"font-size:"		, $font_size, ";font-family: ", $text_font_family, "; stroke-width:", $text_stroke_width	,";\" ",
					"onmouseover=\""	, $text_cursor ," ' ", $aref_sorted_top->[1] ," '; ", $text_color_cursor, $text_curser_fill ,"');\" ",
					"onmouseout=\""		, $text_cursor ," ' ';", $text_color_cursor, $text_fill ,"');\" fill=\"", $text_fill ,"\" > Q2 </text>" ;

	# print upper corner code
	print OUTsvg	"\n\t<text x=\""	, $x_y_corner3[0]-$text_font_size*0.22	,"\" y=\""						, $x_y_corner3[1]-$text_font_size*0.25	,
					"\" style=\"font-size:"		, $font_size, ";font-family: ", $text_font_family, "; stroke-width:", $text_stroke_width	,";\" ",
					"onmouseover=\""	, $text_cursor ," ' ", $aref_sorted_top->[2] ," '; ", $text_color_cursor, $text_curser_fill ,"');\" ",
					"onmouseout=\""		, $text_cursor ," ' ';", $text_color_cursor, $text_fill ,"');\" fill=\"", $text_fill ,"\" > Q3 </text>\n" ;
	################################################### END Corner Legend



	################################################### START Data Print

	#################
	# Print single Data
	for (0 .. @$aref_lol_data-1){

		$lol_x_y[$_][0]	=	$x_y_corner1[0] + $lol_x_y[$_][0] * $zoom_factor;
		$lol_x_y[$_][1]	=	$x_y_corner1[1] - $lol_x_y[$_][1] * $zoom_factor;

		#print "\nData: ", $_, "\tx: ", $lol_x_y[$_][0], "\ty: ", $lol_x_y[$_][1] ;

		my @legend_line_values	=	split "\n", $href_name_of_repeat->{$_};

		print OUTsvg	"\n\t<circle cx=\""	, $lol_x_y[$_][0]		, "\" cy=\""			, $lol_x_y[$_][1]	,
						"\" r=\""			, $circle_radius		,
						"\" stroke=\""		, $legend_line_values[4]	, "\" stroke-width=\""	, $circle_stroke_width	,
						"\" onmouseover=\""	;

		for my $legendline ( 0 .. 3){

			my $line_ref	=	"spktText".$legendline ;
			print OUTsvg		"evt.target.ownerDocument.getElementById('", $line_ref ,"').firstChild.nodeValue= '"	, $legend_line_values[$legendline]	,	"'; "
		}

		print OUTsvg	$text_color_cursor, $circle_curser_fill	,	"');\" onmouseout=\""	;

		for my $legendline ( 0 .. 3){

			my $line_ref	=	"spktText".$legendline ;
			print OUTsvg		"evt.target.ownerDocument.getElementById('"	, $line_ref ,"').firstChild.nodeValue= ' ';"
		}

		print OUTsvg	$text_color_cursor, $legend_line_values[4]	,	"');"	,
						"\" fill=\""		, $legend_line_values[4]	, "\" />"	;
	}
	#################

	#################
	# Print median data
	if ( $$sref_median_print_code eq 'all_quartets' ){

		$median_coordinates[0]	=	$x_y_corner1[0] + $median_coordinates[0] * $zoom_factor;
		$median_coordinates[1]	=	$x_y_corner1[1] - $median_coordinates[1] * $zoom_factor;

		my @legend_data_median	=	( "Support (Overall): ", "Median Support Q1: ".$href_medi_of_top->{0}, "Median Support Q2: ".$href_medi_of_top->{1}, "Median Support Q3: ".$href_medi_of_top->{2}, "tomato" ) ;

		print OUTsvg	"\n\t<circle cx=\""	, $median_coordinates[0]			, "\" cy=\""			, $median_coordinates[1]	,
						"\" r=\""			, $circle_radius_mean		,
						"\" stroke=\""		, $legend_data_median[4]	, "\" stroke-width=\""	, $circle_stroke_width		, "\" onmouseover=\""	;

		for my $legendline ( 0 .. 3){

				my $line_ref	=	"spktText".$legendline ;
				print OUTsvg	"evt.target.ownerDocument.getElementById('", $line_ref ,"').firstChild.nodeValue= '"	, $legend_data_median[$legendline]	,	"'; "
		}

		print OUTsvg	$text_color_cursor, $circle_curser_fill	,	"');\" onmouseout=\""	;

		for my $legendline ( 0 .. 3){

			my $line_ref	=	"spktText".$legendline ;
			print OUTsvg		"evt.target.ownerDocument.getElementById('"	, $line_ref ,"').firstChild.nodeValue= ' ';"
		}

		print OUTsvg	$text_color_cursor, $legend_data_median[4]	,	"');"	,
						"\" fill=\""		, $legend_data_median[4]	, "\" />"	;
	}
	#################

	#################
	# Print mean data
	$mean_coordinates[0]	=	$x_y_corner1[0] + $mean_coordinates[0] * $zoom_factor;
	$mean_coordinates[1]	=	$x_y_corner1[1] - $mean_coordinates[1] * $zoom_factor;

	my @legend_data_mean	=	( "Support (Overall): ", $$sref_mean_median_legend." Support Q1: ".$href_mean_of_top->{0}, $$sref_mean_median_legend." Support Q2: ".$href_mean_of_top->{1}, $$sref_mean_median_legend." Support Q3: ".$href_mean_of_top->{2}, "red" ) ;

	print OUTsvg	"\n\t<circle cx=\""	, $mean_coordinates[0]			, "\" cy=\""			, $mean_coordinates[1]	,
					"\" r=\""			, $circle_radius_mean		,
					"\" stroke=\""		, $legend_data_mean[4]	, "\" stroke-width=\""	, $circle_stroke_width		, "\" onmouseover=\""	;

	for my $legendline ( 0 .. 3){

			my $line_ref	=	"spktText".$legendline ;
			print OUTsvg	"evt.target.ownerDocument.getElementById('", $line_ref ,"').firstChild.nodeValue= '"	, $legend_data_mean[$legendline]	,	"'; "
	}

	print OUTsvg	$text_color_cursor, $circle_curser_fill	,	"');\" onmouseout=\""	;

	for my $legendline ( 0 .. 3){

		my $line_ref	=	"spktText".$legendline ;
		print OUTsvg		"evt.target.ownerDocument.getElementById('"	, $line_ref ,"').firstChild.nodeValue= ' ';"
	}

	print OUTsvg	$text_color_cursor, $legend_data_mean[4]	,	"');"	,
					"\" fill=\""		, $legend_data_mean[4]	, "\" />"	;
	#################

	################################################### END Data Print

	################################################### Start Legend Print
	# Print legend
	my $linespace	=	$text_font_size ;
	my @linecolor	=	( $text_curser_fill, 'black', 'black', 'black' );
	for (0 .. 3){

		print OUTsvg	"\n\t<text id=\"spktText", $_ ,"\" x=\""	, 20	,"\" y=\""						, $x_y_corner1[1]+$linespace	,
						"\" style=\" font-size:10;font-family: Arial; stroke-width:1; fill:", $linecolor[$_], "\" >  </text>" ;

		$linespace += 12
	}
	################################################### END Legend Print

	################################################### Start Title Print
	# Print title
	print OUTsvg	"\n\t<text x=\"20\" y=\""						, $x_y_corner3[1]-($text_font_size*1.5)	,
					"\" style=\"font-size:"		, $font_size, ";font-family: ", $text_font_family, "; stroke-width:", $text_stroke_width	,"; fill=", $text_fill ,"\" > ", $$sref_graphic_title ," </text>" ;
	################################################### End Title Print



	#################
	# print out svg end line and close outfile
	print OUTsvg				"\n</svg>\n" ;
	close OUTsvg ;
	#################

	#################
	# Print terminal info
	print "\n\tPrint SVG ", $$sref_outname ;
	#################
}

sub print_svg_split_graphic{

	my	$href_colours_grafik		=	$_[0]		;	#	key: topology support (0,1, or 2); value: assigned svg color					In -> Defined, OUT -> Unchanged
	my	$aref_quartets_all			=	$_[1]		;	#	list of sorted clade topologies																				In -> Defined, OUT -> Unchanged
	my	$sref_outfilename				=	$_[2]		;	#	outfilename splitgraphic																							In -> Defined, OUT -> Unchanged
	my	$sref_graphic_title			=	$_[3]		;	#	title of svg graphic																									In -> Defined, OUT -> Unchanged
	my	$href_value_of_top			=	$_[4]		;	#	key: toponumber (0, 1, or 2); value: mean or median split score				In -> Defined, OUT -> Unchanged
	my	$href_graphic_setup			=	$_[5]		;	#	key: svgoption; value: option setup																		In -> Defined, OUT -> Unchanged
	my	$href_N_of_top					=	$_[6]		;	#	key: toponumber (0, 1, or 2); value: N resolution											In -> Defined, OUT -> Unchanged
	my	$sref_legend_text				=	$_[7]		;	#	value code of figure legend print out																	In -> Defined, OUT -> Unchanged
	my	$href_q_number_of_top		=	$_[8]		;	#	key alphabetically sorted topology, value: Q code											In -> Defined, OUT -> Unchanged
	my	$sref_name_folder				=	$_[9]		;	#	key: subfolder code; value: subfolder name														In -> Defined, OUT -> Unchanged
	my	$sref_r_option					=	$_[10]	;	#	value (0 -> no directive multiple clan analysis, 1 -> spin directive)	In -> Defined, OUT -> Unchanged
	my	$sref_svg_background		=	$_[11]	;	#	background color svg																									In -> Defined, OUT -> Unchanged
	my	$sref_header_color			=	$_[12]	;	#	header color																													In -> Defined, OUT -> Unchanged


	##############################
	# Assign Grafik colours to possible quartets
	my	%color_of_topology ;	# key: quartet topology; value: grafik colour
		$color_of_topology{$aref_quartets_all->[0]}	= $href_colours_grafik->{$aref_quartets_all->[0]} ;
		$color_of_topology{$aref_quartets_all->[1]}	= $href_colours_grafik->{$aref_quartets_all->[1]} ;
		$color_of_topology{$aref_quartets_all->[2]}	= $href_colours_grafik->{$aref_quartets_all->[2]} ;


	##############################



	################################
	# identify clade names of highest supported split for print out
	( my $clades_best_supported	= $aref_quartets_all->[0] ) =~ s/\(|\)//g ;
	my @cladetaxa_best_support	= split ",", $clades_best_supported ;
	################################



	################################
	# Definement of vector coordinates for single graphic lines
	my @coordinates_terminal_branch_ul		= ( 50, 75, 150, 175 ) ;																																																				# Definement upper, left terminal branch ( x1, y1, x2, y2 coordinates )
	my @coordinates_internal_branch_lr		= ( $coordinates_terminal_branch_ul[2], $coordinates_terminal_branch_ul[3], $coordinates_terminal_branch_ul[2]+100*$href_value_of_top->{$aref_quartets_all->[0]}	, $coordinates_terminal_branch_ul[3] 							) ;	# Definement internal branch upper horizontal ( x1, y1, x2, y2 coordinates )
	my @coordinates_internal_branch_vu		= ( $coordinates_terminal_branch_ul[2], $coordinates_terminal_branch_ul[3], $coordinates_terminal_branch_ul[2]								, $coordinates_terminal_branch_ul[3]+100*$href_value_of_top->{$aref_quartets_all->[1]}	) ;	# Definement internal branch upper left vertical ( x1, y1, x2, y2 coordinates )
	my @coordinates_internal_branch_rl		= ( $coordinates_terminal_branch_ul[2], $coordinates_internal_branch_vu[3], $coordinates_internal_branch_lr[2]								, $coordinates_internal_branch_vu[3] 													) ;	# Definement internal branch lower horizontal ( x1, y1, x2, y2 coordinates )
	my @coordinates_terminal_branch_ur		= ( $coordinates_internal_branch_lr[2], $coordinates_terminal_branch_ul[3], $coordinates_internal_branch_lr[2]+100						, $coordinates_terminal_branch_ul[1] 													) ;	# Definement upper, right terminal branch ( x1, y1, x2, y2 coordinates )
	my @coordinates_terminal_branch_ll		= ( $coordinates_internal_branch_vu[2], $coordinates_internal_branch_vu[3], $coordinates_terminal_branch_ul[0]								, $coordinates_internal_branch_vu[3]+100												) ;	# Definement lower, left terminal branch ( x1, y1, x2, y2 coordinates )
	my @coordinates_terminal_branch_lr		= ( $coordinates_internal_branch_rl[2], $coordinates_internal_branch_rl[3], $coordinates_terminal_branch_ur[2]								, $coordinates_terminal_branch_ll[3]													) ;	# Definement lower, right terminal branch ( x1, y1, x2, y2 coordinates )
	my @coordinates_internal_branch_uv		= ( $coordinates_internal_branch_lr[2], $coordinates_internal_branch_lr[3], $coordinates_internal_branch_rl[2]								, $coordinates_internal_branch_rl[3]													) ;	# Definement internal branch upper right vertical ( x1, y1, x2, y2 coordinates )
	my @coordinates_internal_branch_vl		= ( $coordinates_internal_branch_vu[2], $coordinates_internal_branch_vu[3], $coordinates_internal_branch_vu[2]								, $coordinates_internal_branch_vu[3]+100*$href_value_of_top->{$aref_quartets_all->[2]}	) ;	# Definement internal branch lower left vertical ( x1, y1, x2, y2 coordinates )
	my @coordinates_internal_branch_ho		= ( $coordinates_internal_branch_vl[2], $coordinates_internal_branch_vl[3], $coordinates_internal_branch_lr[2]								, $coordinates_internal_branch_vl[3]													) ; # Definement internal branch upper horizontal 2 ( x1, y1, x2, y2 coordinates )
	my @coordinates_internal_branch_lv		= ( $coordinates_internal_branch_rl[2], $coordinates_internal_branch_rl[3], $coordinates_internal_branch_ho[2]								, $coordinates_internal_branch_ho[3]													) ; # Definement internal branch lower right vertical ( x1, y1, x2, y2 coordinates )

	my @coordinates_taxongroup_1			= ( 35,										70										) ; # Taxon Text upper left terminal branch (I)
	my @coordinates_taxongroup_2			= ( 35,										$coordinates_terminal_branch_ll[3]+12	) ; # Taxon Text lower left terminal branch (II)
	my @coordinates_taxongroup_3			= ( $coordinates_terminal_branch_ur[2]-15,	70										) ; # Taxon Text upper right terminal branch (III)
	my @coordinates_taxongroup_4			= ( $coordinates_terminal_branch_ur[2]-15,	$coordinates_terminal_branch_ll[3]+12	) ; # Taxon Text lower right terminal branch (IV)

	my @coordinates_taxonnames_1			= ( 0,	$coordinates_taxongroup_2[1]+20 ) ; # Subclade print group I
	my @coordinates_taxonnames_2			= ( 0,	$coordinates_taxonnames_1[1]+12 ) ; # Subclade print group II
	my @coordinates_taxonnames_3			= ( 0,	$coordinates_taxonnames_2[1]+12 ) ; # Subclade print group III
	################################


	######################################
	# Standard values of grafik output considering subtree circles, terminal branches, and internal node circle radius
	my $stroke_width_branches				=	$href_graphic_setup->{stroke_width_branches}	;
	my $stroke_color_branches				=	$href_graphic_setup->{stroke_color_branches}	;
	my $axis_stroke_opacity					=	$href_graphic_setup->{axis_stroke_opacity}		;
	my $axis_fill_opacity						=	$href_graphic_setup->{axis_fill_opacity}			;
	my $axis_fill_color							=	$href_graphic_setup->{axis_fill_color}				;
	my $text_font_size							=	$href_graphic_setup->{text_font_size}					;
	my $text_stroke									=	$href_graphic_setup->{text_stroke}						;
	my $text_stroke_2								=	$href_graphic_setup->{text_stroke_2}					;
	my $text_font_weight						=	$href_graphic_setup->{text_font_weight}				;
	my $text_stroke_width						=	$href_graphic_setup->{text_stroke_width}			;
	my $text_fill										=	$href_graphic_setup->{text_fill}							;
	my $text_fill_2									=	$href_graphic_setup->{text_fill_2}						;
	my $text_font_family						=	$href_graphic_setup->{text_font_family}				;
	my $text_curser_fill						=	$href_graphic_setup->{text_curser_fill}				;
	my $scale_color									=	$href_graphic_setup->{scale_color}						;
	my $scale_stroke_width					=	$href_graphic_setup->{scale_stroke_width}			;
	my $scale_std_value							=	$href_graphic_setup->{scale_std_value}				;
	my $scale_font_size							=	$href_graphic_setup->{scale_font_size}				;
	my $legend_font_size						=	$href_graphic_setup->{legend_font_size}				;
	my $legend_font_weight					=	$href_graphic_setup->{legend_font_weight}			;
	my $legend_stroke_width					=	$href_graphic_setup->{legend_stroke_width}		;
	my $split_factor								=	$href_graphic_setup->{split_factor}						;
	my $font_size										=	$href_graphic_setup->{font_size}							;
	######################################


	######################################
	# Standard values of grafik output considering subtree circles, terminal branches, and internal node circle radius
	#my $stroke_width_branches				=	1	;
	#my $stroke_color_branches				=	"black"	;
	#my $axis_stroke_opacity					= "1.000000"	;
	#my $axis_fill_opacity					= "0.000000"	;
	#my $axis_fill_color						= "#000000"		;
	#my $text_font_size						=	30	;
	#my $text_stroke							=	"black"	;
	#my $text_stroke_2						=	"orange"	;
	#my $text_font_weight					=	"normal"	;	# bold
	#my $text_stroke_width					=	0	;
	#my $text_fill							=	"black"	;
	#my $text_fill_2							=	"orange"	;
	#my $text_font_family					=	"Arial"	;
	#my $text_curser_fill					=	"grey"		;
	#my $scale_color							=	"black"	;
	#my $scale_stroke_width					=	1	;
	#my $scale_std_value						=	0.1	;
	#my $scale_font_size						=	12	;
	#my $legend_font_size					=	12	;
	#my $legend_font_weight					=	"normal"	;
	#my $legend_stroke_width					=	0	;
	#my $split_factor						=	200	;
	#my $font_size							=	12			;
	######################################


	#################
	# Definement of docu dimensions
	my	$zoom_factor	=	100 ;

	my	$docu_width		=	8 * $zoom_factor ;
	my	$docu_height	=	4 * $zoom_factor ;
	#################
	#my $docu_width	= 1000 ;
	#my $docu_height	= $coordinates_taxonnames_3[1]+12 ;



	######################################
	# Print OUT vector graphic
	open	OUTsvg, ">$$sref_name_folder/$$sref_outfilename" or die "\n\tOUTFILE-ERROR: Cannot open split graphic ", $$sref_name_folder, "/", $$sref_outfilename, " in subroutine &print_mice_split_graphic!\n" ;

	########### SVG Header
	print	OUTsvg	"<?xml version=\"1.0\"?>\n",
					"<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n",
					"\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n\n",
					"<svg width=\"", $docu_width, "\" height=\"", $docu_height, "\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">\n\n" ;

	if ( $$sref_r_option == 1 ){

		print	OUTsvg	"\t<defs>\n",
						"<marker id=\"arrow\" markerWidth=\"10\" markerHeight=\"10\" refX=\"0\" refY=\"2.5\" orient=\"auto\" markerUnits=\"strokeWidth\">\n",
						"<path d=\"M0,0 L0,5 L9,2.5 z\" fill=\"grey\" />" ,
						"</marker>",
						"\t</defs>\n\n"
	}
	else{

		print OUTsvg	"\t<defs>\n",
						"\n",
						"\t</defs>\n\n"
	}

	################################################### Start Title Print
	# Print title
	print OUTsvg	"\n\t<rect width=\"", $docu_width, "\" height=\"", $docu_height, "\" fill=\"", $$sref_svg_background, "\" />" ,
					"\n\t<text x=\"35\" y=\"35\" style=\"font-size:12;font-family: Arial; stroke-width:1; fill:", $$sref_header_color, "\" > ", $$sref_graphic_title ," </text>" ;
	################################################### End Title Print


	########### SVG Terminal Branches
	# upper left:
	print OUTsvg	"\n\t<line x1=\"", $coordinates_terminal_branch_ul[0], "\" y1=\"", $coordinates_terminal_branch_ul[1], "\" x2=\"", $coordinates_terminal_branch_ul[2], "\" y2=\"", $coordinates_terminal_branch_ul[3], "\" style=\"stroke-width:", $stroke_width_branches, ";stroke:", $stroke_color_branches, ";fill:", $axis_fill_color, ";stroke-opacity:",$axis_stroke_opacity ,";fill-opacity:", $axis_fill_opacity, "\"/>" ;

	# upper right:
	print OUTsvg	"\n\t<line x1=\"", $coordinates_terminal_branch_ur[0], "\" y1=\"", $coordinates_terminal_branch_ur[1], "\" x2=\"", $coordinates_terminal_branch_ur[2], "\" y2=\"", $coordinates_terminal_branch_ur[3], "\" style=\"stroke-width:", $stroke_width_branches, ";stroke:", $stroke_color_branches, ";fill:", $axis_fill_color, ";stroke-opacity:",$axis_stroke_opacity ,";fill-opacity:", $axis_fill_opacity, "\"/>" ;

	# lower left:
	print OUTsvg	"\n\t<line x1=\"", $coordinates_terminal_branch_ll[0], "\" y1=\"", $coordinates_terminal_branch_ll[1], "\" x2=\"", $coordinates_terminal_branch_ll[2], "\" y2=\"", $coordinates_terminal_branch_ll[3], "\" style=\"stroke-width:", $stroke_width_branches, ";stroke:", $stroke_color_branches, ";fill:", $axis_fill_color, ";stroke-opacity:",$axis_stroke_opacity ,";fill-opacity:", $axis_fill_opacity, "\"/>" ;

	# lower right:
	print OUTsvg	"\n\t<line x1=\"", $coordinates_terminal_branch_lr[0], "\" y1=\"", $coordinates_terminal_branch_lr[1], "\" x2=\"", $coordinates_terminal_branch_lr[2], "\" y2=\"", $coordinates_terminal_branch_lr[3], "\" style=\"stroke-width:", $stroke_width_branches, ";stroke:", $stroke_color_branches, ";fill:", $axis_fill_color, ";stroke-opacity:",$axis_stroke_opacity ,";fill-opacity:", $axis_fill_opacity, "\"/>" ;


	########### SVG internal horizontal Branches
	# internal horizontal branch up
	print OUTsvg	"\n\t<line x1=\"", $coordinates_internal_branch_lr[0], "\" y1=\"", $coordinates_internal_branch_lr[1], "\" x2=\"", $coordinates_internal_branch_lr[2], "\" y2=\"", $coordinates_internal_branch_lr[3], "\" style=\"stroke-width:", $stroke_width_branches, ";fill:", $axis_fill_color, ";stroke-opacity:",$axis_stroke_opacity ,";fill-opacity:", $axis_fill_opacity, "\" onmouseover=\"evt.target.ownerDocument.getElementById('spktText0').firstChild.nodeValue= '", $href_q_number_of_top->{$aref_quartets_all->[0]} ,": ", $aref_quartets_all->[0] ,"'; evt.target.ownerDocument.getElementById('spktText1').firstChild.nodeValue= '", $$sref_legend_text ," Support: ", $href_N_of_top->{$aref_quartets_all->[0]}, "'; evt.target.setAttributeNS(null, 'stroke', 'grey');\" onmouseout=\"evt.target.ownerDocument.getElementById('spktText0').firstChild.nodeValue= ' ';evt.target.ownerDocument.getElementById('spktText1').firstChild.nodeValue= ' '; evt.target.setAttributeNS(null, 'stroke', '", $color_of_topology{$aref_quartets_all->[0]} ,"');\" stroke=\"", $color_of_topology{$aref_quartets_all->[0]}, "\" />" ;

	# internal horizontal branch middle
	print OUTsvg	"\n\t<line x1=\"", $coordinates_internal_branch_rl[0], "\" y1=\"", $coordinates_internal_branch_rl[1], "\" x2=\"", $coordinates_internal_branch_rl[2], "\" y2=\"", $coordinates_internal_branch_rl[3], "\" style=\"stroke-width:", $stroke_width_branches, ";fill:", $axis_fill_color, ";stroke-opacity:",$axis_stroke_opacity ,";fill-opacity:", $axis_fill_opacity, "\" onmouseover=\"evt.target.ownerDocument.getElementById('spktText0').firstChild.nodeValue= '", $href_q_number_of_top->{$aref_quartets_all->[0]} ,": ", $aref_quartets_all->[0] ,"'; evt.target.ownerDocument.getElementById('spktText1').firstChild.nodeValue= '", $$sref_legend_text ," Support: ", $href_N_of_top->{$aref_quartets_all->[0]} ,"'; evt.target.setAttributeNS(null, 'stroke', 'grey');\" onmouseout=\"evt.target.ownerDocument.getElementById('spktText0').firstChild.nodeValue= ' ';evt.target.ownerDocument.getElementById('spktText1').firstChild.nodeValue= ' '; evt.target.setAttributeNS(null, 'stroke', '", $color_of_topology{$aref_quartets_all->[0]} ,"');\" stroke=\"", $color_of_topology{$aref_quartets_all->[0]}, "\" />" ;

	# internal horizontal branch bottom
	print OUTsvg	"\n\t<line x1=\"", $coordinates_internal_branch_ho[0], "\" y1=\"", $coordinates_internal_branch_ho[1], "\" x2=\"", $coordinates_internal_branch_ho[2], "\" y2=\"", $coordinates_internal_branch_ho[3], "\" style=\"stroke-width:", $stroke_width_branches, ";fill:", $axis_fill_color, ";stroke-opacity:",$axis_stroke_opacity ,";fill-opacity:", $axis_fill_opacity, "\" onmouseover=\"evt.target.ownerDocument.getElementById('spktText0').firstChild.nodeValue= '", $href_q_number_of_top->{$aref_quartets_all->[0]} ,": ", $aref_quartets_all->[0] ,"'; evt.target.ownerDocument.getElementById('spktText1').firstChild.nodeValue= '", $$sref_legend_text ," Support: ", $href_N_of_top->{$aref_quartets_all->[0]} ,"'; evt.target.setAttributeNS(null, 'stroke', 'grey');\" onmouseout=\"evt.target.ownerDocument.getElementById('spktText0').firstChild.nodeValue= ' ';evt.target.ownerDocument.getElementById('spktText1').firstChild.nodeValue= ' '; evt.target.setAttributeNS(null, 'stroke', '", $color_of_topology{$aref_quartets_all->[0]} ,"');\" stroke=\"", $color_of_topology{$aref_quartets_all->[0]}, "\" />" ;


	########### SVG internal vertical upper Branches
	# internal vertical, upper branch left
	print OUTsvg	"\n\t<line x1=\"", $coordinates_internal_branch_vu[0], "\" y1=\"", $coordinates_internal_branch_vu[1], "\" x2=\"", $coordinates_internal_branch_vu[2], "\" y2=\"", $coordinates_internal_branch_vu[3], "\" style=\"stroke-width:", $stroke_width_branches, ";fill:", $axis_fill_color, ";stroke-opacity:",$axis_stroke_opacity ,";fill-opacity:", $axis_fill_opacity, "\" onmouseover=\"evt.target.ownerDocument.getElementById('spktText0').firstChild.nodeValue= '", $href_q_number_of_top->{$aref_quartets_all->[1]} ,": ", $aref_quartets_all->[1] ,"'; evt.target.ownerDocument.getElementById('spktText1').firstChild.nodeValue= '", $$sref_legend_text ," Support: ", $href_N_of_top->{$aref_quartets_all->[1]} ,"'; evt.target.setAttributeNS(null, 'stroke', 'grey');\" onmouseout=\"evt.target.ownerDocument.getElementById('spktText0').firstChild.nodeValue= ' ';evt.target.ownerDocument.getElementById('spktText1').firstChild.nodeValue= ' '; evt.target.setAttributeNS(null, 'stroke', '", $color_of_topology{$aref_quartets_all->[1]} ,"');\" stroke=\"", $color_of_topology{$aref_quartets_all->[1]}, "\" />" ;

	# internal vertical, upper branch right
	print OUTsvg	"\n\t<line x1=\"", $coordinates_internal_branch_uv[0], "\" y1=\"", $coordinates_internal_branch_uv[1], "\" x2=\"", $coordinates_internal_branch_uv[2], "\" y2=\"", $coordinates_internal_branch_uv[3], "\" style=\"stroke-width:", $stroke_width_branches, ";fill:", $axis_fill_color, ";stroke-opacity:",$axis_stroke_opacity ,";fill-opacity:", $axis_fill_opacity, "\" onmouseover=\"evt.target.ownerDocument.getElementById('spktText0').firstChild.nodeValue= '", $href_q_number_of_top->{$aref_quartets_all->[1]} ,": ", $aref_quartets_all->[1] ,"'; evt.target.ownerDocument.getElementById('spktText1').firstChild.nodeValue= '", $$sref_legend_text ," Support: ", $href_N_of_top->{$aref_quartets_all->[1]} ,"'; evt.target.setAttributeNS(null, 'stroke', 'grey');\" onmouseout=\"evt.target.ownerDocument.getElementById('spktText0').firstChild.nodeValue= ' ';evt.target.ownerDocument.getElementById('spktText1').firstChild.nodeValue= ' '; evt.target.setAttributeNS(null, 'stroke', '", $color_of_topology{$aref_quartets_all->[1]} ,"');\" stroke=\"", $color_of_topology{$aref_quartets_all->[1]}, "\" />" ;


	########### SVG internal vertical lower Branches
	# internal vertical, lower branch left
	print OUTsvg	"\n\t<line x1=\"", $coordinates_internal_branch_vl[0], "\" y1=\"", $coordinates_internal_branch_vl[1], "\" x2=\"", $coordinates_internal_branch_vl[2], "\" y2=\"", $coordinates_internal_branch_vl[3], "\" style=\"stroke-width:", $stroke_width_branches, ";fill:", $axis_fill_color, ";stroke-opacity:",$axis_stroke_opacity ,";fill-opacity:", $axis_fill_opacity, "\" onmouseover=\"evt.target.ownerDocument.getElementById('spktText0').firstChild.nodeValue= '", $href_q_number_of_top->{$aref_quartets_all->[2]} ,": ", $aref_quartets_all->[2] ,"'; evt.target.ownerDocument.getElementById('spktText1').firstChild.nodeValue= '", $$sref_legend_text ," Support: ", $href_N_of_top->{$aref_quartets_all->[2]} ,"'; evt.target.setAttributeNS(null, 'stroke', 'grey');\" onmouseout=\"evt.target.ownerDocument.getElementById('spktText0').firstChild.nodeValue= ' ';evt.target.ownerDocument.getElementById('spktText1').firstChild.nodeValue= ' '; evt.target.setAttributeNS(null, 'stroke', '", $color_of_topology{$aref_quartets_all->[2]} ,"');\" stroke=\"", $color_of_topology{$aref_quartets_all->[2]}, "\" />" ;

	# internal vertical, lower branch right
	print OUTsvg	"\n\t<line x1=\"", $coordinates_internal_branch_lv[0], "\" y1=\"", $coordinates_internal_branch_lv[1], "\" x2=\"", $coordinates_internal_branch_lv[2], "\" y2=\"", $coordinates_internal_branch_lv[3], "\" style=\"stroke-width:", $stroke_width_branches, ";fill:", $axis_fill_color, ";stroke-opacity:",$axis_stroke_opacity ,";fill-opacity:", $axis_fill_opacity, "\" onmouseover=\"evt.target.ownerDocument.getElementById('spktText0').firstChild.nodeValue= '", $href_q_number_of_top->{$aref_quartets_all->[2]} ,": ", $aref_quartets_all->[2] ,"'; evt.target.ownerDocument.getElementById('spktText1').firstChild.nodeValue= '", $$sref_legend_text ," Support: ", $href_N_of_top->{$aref_quartets_all->[2]} ,"'; evt.target.setAttributeNS(null, 'stroke', 'grey');\" onmouseout=\"evt.target.ownerDocument.getElementById('spktText0').firstChild.nodeValue= ' ';evt.target.ownerDocument.getElementById('spktText1').firstChild.nodeValue= ' '; evt.target.setAttributeNS(null, 'stroke', '", $color_of_topology{$aref_quartets_all->[2]} ,"');\" stroke=\"", $color_of_topology{$aref_quartets_all->[2]}, "\" />" ;


	########### SVG Terminal Node Definitions
	# Taxon I
	print OUTsvg "\t<text x=\"", $coordinates_taxongroup_1[0] ,"\" y=\"", $coordinates_taxongroup_1[1] ,"\" style=\"font-size:"	, $font_size, ";font-family: ", $text_font_family, "; stroke-width:", $text_stroke_width	,";fill:black\" > ", $cladetaxa_best_support[0] ," </text>\n" ;

	# Taxon II
	print OUTsvg "\t<text x=\"", $coordinates_taxongroup_2[0] ,"\" y=\"", $coordinates_taxongroup_2[1] ,"\" style=\"font-size:"	, $font_size, ";font-family: ", $text_font_family, "; stroke-width:", $text_stroke_width	,";fill:black\" > ", $cladetaxa_best_support[1] ," </text>\n" ;

	# Taxon III
	print OUTsvg "\t<text x=\"", $coordinates_taxongroup_3[0] ,"\" y=\"", $coordinates_taxongroup_3[1] ,"\" style=\"font-size:"	, $font_size, ";font-family: ", $text_font_family, "; stroke-width:", $text_stroke_width	,";fill:black\" > ", $cladetaxa_best_support[2] ," </text>\n" ;

	# Taxon IV
	print OUTsvg "\t<text x=\"", $coordinates_taxongroup_4[0] ,"\" y=\"", $coordinates_taxongroup_4[1] ,"\" style=\"font-size:"	, $font_size, ";font-family: ", $text_font_family, "; stroke-width:", $text_stroke_width	,";fill:black\" > ", $cladetaxa_best_support[3] ," </text>\n" ;


	########### Spin Arrow
	if ( $$sref_r_option == 1 ){

		my $x2 = $coordinates_internal_branch_lr[2]-9 ;
		print OUTsvg "\n\t<line x1=\"", $coordinates_internal_branch_lr[0], "\" y1=\"", $coordinates_terminal_branch_ul[1], "\" x2=\"", $x2 , "\" y2=\"", $coordinates_terminal_branch_ul[1], , "\" style=\"stroke-width:", $stroke_width_branches, ";stroke:", $stroke_color_branches, ";fill:", $axis_fill_color, ";stroke-opacity:",$axis_stroke_opacity ,";fill-opacity:", $axis_fill_opacity, "\" marker-end=\"url(#arrow)\"/>\n" }
	###########


	########### SVG Subclade Definitions
	# Subclade I
	#print OUTsvg "\t<text x=\"", $coordinates_taxonnames_1[0] ,"\" y=\"", $coordinates_taxonnames_1[1] ,"\" style=\"font-size:10;font-family: Arial;fill:", $color_of_topology{$sorted_topologies_by_value[0]} ,"\" > ", $sorted_topologies_by_value[0]," Mean Weight: ", $mean_score_of_top{$sorted_topologies_by_value[0]}/100 ," </text>\n" ;

	# Subclade II
	#print OUTsvg "\t<text x=\"", $coordinates_taxonnames_2[0] ,"\" y=\"", $coordinates_taxonnames_2[1] ,"\" style=\"font-size:10;font-family: Arial;fill:", $color_of_topology{$sorted_topologies_by_value[1]} ,"\" > ", $sorted_topologies_by_value[1]," Mean Weight: ", $mean_score_of_top{$sorted_topologies_by_value[1]}/100 ," </text>\n" ;

	# Subclade III
	#print OUTsvg "\t<text x=\"", $coordinates_taxonnames_3[0] ,"\" y=\"", $coordinates_taxonnames_3[1] ,"\" style=\"font-size:10;font-family: Arial;fill:", $color_of_topology{$sorted_topologies_by_value[2]} ,"\" > ", $sorted_topologies_by_value[2]," Mean Weight: ", $mean_score_of_top{$sorted_topologies_by_value[2]}/100 ," </text>\n" ;


	################################################### Start Legend Print
	# Print legend
	my $linespace	=	$text_font_size ;
	my @linecolor	=	( $text_curser_fill, 'black' );

	print OUTsvg	"\n\t<text id=\"spktText0\" x=\"35\" y=\""	, $coordinates_taxonnames_1[1]	, "\" style=\" font-size:10;font-family: Arial; stroke-width:1; fill:", $linecolor[0], "\" >  </text>" ;
	print OUTsvg	"\n\t<text id=\"spktText1\" x=\"35\" y=\""	, $coordinates_taxonnames_2[1]	, "\" style=\" font-size:10;font-family: Arial; stroke-width:1; fill:", $linecolor[1], "\" >  </text>" ;
	################################################### END Legend Print


	##############################
	# print out svg end line and close outfile
	print OUTsvg				"\n</svg>\n" ;
	close OUTsvg ;
	##############################


	##############################
	# Print terminal info
	print "\n\tPrint SVG ", $$sref_outfilename ;
	##############################


	######################################
}

sub print_svg_barchart_graphic{

	my	$href_colours_grafik	=	$_[0]	;	#	key: topology support (0,1, or 2); value: assigned svg color	In -> Defined, OUT -> Unchanged
	my	$aref_trees_six			=	$_[1]	;	#	list of sorted clade topologies									In -> Defined, OUT -> Unchanged
	my	$href_hoh_seen_quartet	=	$_[2]	;	#	key: sublade topology; value: N occurence						In -> Defined, OUT -> Unchanged
	my	$sref_Nquartets			=	$_[3]	;	#	total number of quartets										In -> Defined, OUT -> Unchanged
	my	$sref_graphic_title		=	$_[4]	;	#	svg graphic title												In -> Defined, OUT -> Unchanged
	my	$sref_outfilename_bar	=	$_[5]	;	#	svg outfile name												In -> Defined, OUT -> Unchanged
	my	$aref_score_number		=	$_[6]	;	#	list number of each quartet topology (0,1,2)					In -> Defined, OUT -> Unchanged
	my	$href_graphic_setup		=	$_[7]	;	#	key: svgoption; value: option setup								In -> Defined, OUT -> Unchanged
	my	$sref_name_folder		=	$_[8]	;	#	key: subfolder code; value: subfolder name						In -> Defined, OUT -> Unchanged
	my	$sref_svg_background	=	$_[9]	;	#	background color svg											In -> Defined, OUT -> Unchanged
	my	$sref_header_color		=	$_[10]	;	#	header color													In -> Defined, OUT -> Unchanged
	my	$sref_bar_color			=	$_[11]	;	#	barplot color of lower supported spin trees						In -> Defined, OUT -> Unchanged

	my	$stroke_width_branches				=	$href_graphic_setup->{stroke_width_branches}	;
	my	$bar_stroke_width					=	$href_graphic_setup->{bar_stroke_width}			;
	my	$text_stroke_width					=	$href_graphic_setup->{text_stroke_width}		;
	my	$text_font_family					=	$href_graphic_setup->{text_font_family}			;
	my	$legend_stroke_width				=	$href_graphic_setup->{legend_stroke_width}		;
	my	$raster_stroke_width				=	$href_graphic_setup->{raster_stroke_width}		;
	my	$axis_stroke_color					=	$href_graphic_setup->{axis_stroke_color}		;
	my	$axis_fill_color					=	$href_graphic_setup->{axis_fill_color}			;
	my	$axis_stroke_opacity				=	$href_graphic_setup->{axis_stroke_opacity}		;
	my	$axis_fill_opacity					=	$href_graphic_setup->{axis_fill_opacity}		;
	my	$raster_stroke_color				=	$href_graphic_setup->{raster_stroke_color}		;
	my	$raster_stroke_opacity				=	$href_graphic_setup->{raster_stroke_opacity}	;
	my	$raster_fill_color					=	$href_graphic_setup->{raster_fill_color}		;
	my	$raster_fill_opacity				=	$href_graphic_setup->{raster_fill_opacity}		;
	my	$bar_stroke_color					=	$href_graphic_setup->{bar_stroke_color}			;
	my	$bar_fill_color						=	$href_graphic_setup->{bar_fill_color}			;
	my	$bar_stroke_opacity					=	$href_graphic_setup->{bar_stroke_opacity}		;
	my	$bar_stroke_fill					=	$href_graphic_setup->{bar_stroke_fill}			;
	my	$text_curser_fill					=	$href_graphic_setup->{text_curser_fill}			;
	my	$font_size							=	$href_graphic_setup->{font_size}				;

	my	$Nbar_plots							=	@$aref_score_number ;
	my	$N_trees							=	@$aref_trees_six ;

	my $N_quartets			= $$sref_Nquartets	;																# Total Number of analysed quartets
	my $x_value				= 50		;																		# start point x value
	my $y_value				= 320	;																			# start point y value
	my $width_bar_plot		= 100	;																			# width length of single bar plot
	my $width_bar_space		= 60	;																			# width length of bar plot space
	my $high_bar_plot		= 200	;																			# maximum high of each bar plot (100%)
	my $axis_space_end		= 10	;																			# axis space at the end of x and y axis
	my $bar_text_distance	= 10	;																			# distance number of observed quartets <-> top of assigned bar plot
	my $length_x_line_yaxis	= 20	;																			# length legend x axis on y axes
	my $legend_descr_space	= 20	;																			# space between single legend descriptions
	my $stroke_dasharray	= 10	;																			# define length stroke dashs
	my $stroke_dash_space	= 2		;																			# define space between stroke dashs
	my $axis_x_start		= $x_value + 40	;																	# x start point x axis
	my $axis_y_start		= $y_value - 20	;																	# y start point y axis
	my $legend_y_line_le	= $axis_x_start  - ($length_x_line_yaxis / 2) ;										# length legend line left to the y axis
	my $legend_y_line_ri	= $axis_x_start  + ($length_x_line_yaxis / 2) ;										# length legend line right to the y axis
	my $length_x_axis		= $axis_x_start  + $Nbar_plots * $width_bar_plot + $Nbar_plots * $width_bar_space + $axis_space_end ;		# length x axis
	my $length_y_axis		= $axis_y_start  - $high_bar_plot - $axis_space_end ;								# length y axis
	my $bar_plot_start_x	= $axis_x_start  + $width_bar_space ;												# x startpoint first bar plot
	my $y_value_header		= $axis_y_start  - $high_bar_plot - $axis_space_end - 20 ;							# y start point header
	my $x_value_header		= $length_x_axis / 2 ;																# x start point header
	my $y_legend_x_axis		= $y_value ;																		# y start point legend x axis
	my $y_legend_descript	= $axis_y_start  - $high_bar_plot ;													# y start point legend description
	my $x_legend_descript	= $length_x_axis ;																	# y start point legend description
	my $x_legend_y_axis		= $x_value ;																		# y start point legend description
	my $y_legend_y_axis		= $length_y_axis - 10 ;																# y start point legend description
	my $docu_width			= $x_legend_descript + 50 ;															# document width
	my $docu_height			= 400 * $N_trees;																	# document height



	##############################
	# Assign Grafik colours to possible quartets
	my	%colour_of_quartet ;	# key: quartet topology; value: grafik colour
	my	$colour_other_trees = $$sref_bar_color ; # light blue

	for my $t ( @$aref_trees_six ){

		if		( $href_colours_grafik->{$t} )	{ $colour_of_quartet{$t}	= $href_colours_grafik->{$t} }
		else 									{ $colour_of_quartet{$t}	= $colour_other_trees }
	}
	##############################



	##############################
	my %legend_of_order		= (
								1	=>	'Best'		,
								2	=>	'2nd'		,
								3	=>	'3rd'		,
								4	=>	'4th'		,
								5	=>	'5th'		,
								6	=>	'Lowest'	,
	) ;
	##############################



	##############################
	# Titel of each barplot
	my $header		=	$$sref_graphic_title ;
	my $bartitle	=	"N Split Score Ranking in Single Quartet Analysis" ;
	##############################



	##############################
	# Open out svg file and svg header print
	my	$outfile			=	$$sref_outfilename_bar ;

	open	OUTsvg, ">$$sref_name_folder/$outfile" or die "\n\tOUTFILE-ERROR: Cannot open split graphic ", $$sref_name_folder, "/", $outfile, " in subroutine &print_mice_reconstruction_bar_svg!\n" ;
	print	OUTsvg	"<?xml version=\"1.0\"?>\n",
					"<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n",
					"\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n\n",
					"<svg width=\"", $docu_width, "\" height=\"", $docu_height, "\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">\n\n",
					"\t<defs>\n",
					"\n",
					"\t</defs>\n",
					"\n",
	;
	print OUTsvg	"\n\t<rect width=\"", $docu_width, "\" height=\"", $docu_height, "\" fill=\"", $$sref_svg_background, "\" />" ;
	print OUTsvg	"\n\t<text x=\"", 0.5*$x_value ,"\" y=\"12\" style=\"font-size:12;font-family: Arial; stroke-width:1; fill:", $$sref_header_color, ";\" > ", $$sref_graphic_title ," </text>" ;
	##############################



	##############################
	# Generate bar charts of each topology given the number of identified order positions
	for my $topology ( @$aref_trees_six ){

		##############################
		# Define x axis start points of each bar plot
		# Store N of occurence of each quartet topology as hasvalue of start point x axis
		# Store quartet topology as hashvalue of start point x axis
		# Store y start point for text above bar plot of start point x axis
		my (
				%bar_plot_y_of_x	,	# keys: x value start of each bar plot (e.g. 60) ; values: y value end of each bar plot (e.g. 180)
				%topo_of_x_start	,	# keys: x value start of each bar plot ; values: quartet topology
				%valu_of_x_start	,	# keys: quartet topology ; values Number of quartet occurency
				%colour_of_x_start	,	# keys: x value start of each bar plot ; value: barplot colour
				%N_text_start		,	# keys: x value start of each bar plot ; values: y value for N occurence print
		) ;

		for my $order_number ( @$aref_score_number ){

			if ( $href_hoh_seen_quartet->{$topology}{$order_number} ){

				$bar_plot_y_of_x{$bar_plot_start_x}	= $axis_y_start - (( $href_hoh_seen_quartet->{$topology}{$order_number} / $$sref_Nquartets	) * $high_bar_plot) ;		# y value of each bar plot (%)
				$valu_of_x_start{$bar_plot_start_x}	= $href_hoh_seen_quartet->{$topology}{$order_number} ;																	# assign number of oberved quartet topologies to bar start point
			}

			else{

				$bar_plot_y_of_x{$bar_plot_start_x}	= $axis_y_start ;	# y value of each bar plot (%)
				$valu_of_x_start{$bar_plot_start_x}	= 0 ;				# assign number of oberved quartet topologies to bar start point
			}

			if		( $colour_of_quartet{$topology} )	{ $colour_of_x_start{$bar_plot_start_x}	= $colour_of_quartet{$topology}	}
			else										{ $colour_of_x_start{$bar_plot_start_x}	= $colour_other_trees			} # assign bar plot colour for given topology

			$N_text_start	{$bar_plot_start_x}		= $bar_plot_y_of_x{$bar_plot_start_x} - $bar_text_distance ;	# assign y start point for text above bar plot
			$bar_plot_start_x += ( $width_bar_plot + $width_bar_space ) ;											# Define new x star point for next quartet
		}
		##############################


		##############################
		# print out y-axis
		print OUTsvg	"\n\t<line x1=\"", $axis_x_start, "\" y1=\"", $axis_y_start, "\" x2=\"", $axis_x_start, "\" y2=\"", $length_y_axis, "\" style=\"stroke:", $axis_stroke_color, ";stroke-width:", $stroke_width_branches, ";\" />" ;
		##############################


		##############################
		# print out x-axis
		print OUTsvg	"\n\t<line x1=\"", $axis_x_start, "\" y1=\"", $axis_y_start, "\" x2=\"", $length_x_axis, "\" y2=\"", $axis_y_start, "\" style=\"stroke:", $axis_stroke_color, ";stroke-width:", $stroke_width_branches, ";\" />" ;
		##############################


		##############################
		# print out legend lines and values y axis
		my $legend_value	= 0 ;
		for ( my $k = $axis_y_start; $k >= ($axis_y_start - $high_bar_plot); $k -= ($high_bar_plot/10) ){

			## legend line
			print OUTsvg	"\n\t<line x1=\"", $legend_y_line_le, "\" y1=\"", $k, "\" x2=\"", $legend_y_line_ri, "\" y2=\"", $k, "\" style=\"stroke-width:", $legend_stroke_width, ";stroke:", $axis_stroke_color, ";fill:", $axis_fill_color, ";stroke-opacity:",$axis_stroke_opacity ,";fill-opacity:", $axis_fill_opacity, "\"/>" ;

			## x raster line
			print OUTsvg	"\n\t<line x1=\"", $axis_x_start, "\" y1=\"", $k, "\" x2=\"", $length_x_axis, "\" y2=\"", $k, "\" style=\"stroke-width:", $raster_stroke_width, ";stroke:", $raster_stroke_color, ";fill:", $raster_fill_color, ";stroke-opacity:",$raster_stroke_opacity ,";fill-opacity:", $raster_fill_opacity, ";stroke-dasharray:", $stroke_dasharray, " ", $stroke_dash_space, "\"/>" ;

			## legend number
			print OUTsvg	"\n\t<text x=\"", $x_value, "\" y=\"", $k ,"\" style=\"font-size:"	, $font_size, ";font-family: ", $text_font_family, "; stroke-width:", $text_stroke_width	,";\"  > ", $legend_value ," </text>" ;

			$legend_value += 10
		}
		##############################


		##############################
		# print out legend y axis
		print OUTsvg	"\n\t<text x=\"", $x_legend_y_axis, "\" y=\"", $y_legend_y_axis ,"\" style=\"font-size:"		, $font_size, ";font-family: ", $text_font_family, "; stroke-width:", $text_stroke_width	,";\"  > % Occurence </text>" ;
		##############################


		##############################
		# print out bar plots
		# print out number of assigned quartet occurence (above bar)
		# print out legend x axis (below bar and x axis)
		my $counter 				= 1 ;
		#my @bar_colors = ( "#FF9999", "#99FF99", "#99CCFF" ) ;
		for my $x ( sort {$a<=>$b} keys %bar_plot_y_of_x ){

			## bar plots
			my $height_bar_plot	= $axis_y_start - $bar_plot_y_of_x{$x} ;
			print OUTsvg	"\n\t<rect x=\"", $x ,"\" y=\"", $bar_plot_y_of_x{$x} , "\" width=\"", $width_bar_plot, "\" height=\"", $height_bar_plot,"\" style=\"stroke-width:", $bar_stroke_width ,";stroke:", $bar_stroke_color, ";fill:", $colour_of_x_start{$x}, ";stroke-opacity:", $bar_stroke_opacity, ";fill-opacity:", $bar_stroke_opacity, "\" />" ;

			## value
			print OUTsvg	"\n\t<text x=\"", $x ,"\" y=\"", $N_text_start{$x} ,"\" style=\"font-size:"		, $font_size, ";font-family: ", $text_font_family, "; stroke-width:", $text_stroke_width	,";\" > N = ", $valu_of_x_start{$x} ," </text>" ;

			## X legend
			my $legend_text = $legend_of_order{$counter} ;
			print OUTsvg	"\n\t<text x=\"", $x, "\" y=\"", $y_legend_x_axis, "\" style=\"font-size:"		, $font_size, ";font-family: ", $text_font_family, "; stroke-width:", $text_stroke_width	,";\" > ", $legend_text ," </text>" ;

			$counter++ ;
		}
		##############################


		################################################### Start Legend Print
		# Print legend
		#my $linespace	=	12 ;
		my @linecolor	=	( $text_curser_fill, 'black' );

		print OUTsvg	"\n\t<text id=\"spktText0\" x=\"", $x_value ,"\" y=\""	, $y_legend_x_axis+24	, "\" style=\" font-size:10;font-family: Arial; stroke-width:1; fill:", $linecolor[0], "\" >  </text>" ;
		#print OUTsvg	"\n\t<text id=\"spktText1\" x=\"", $x_value ,"\" y=\""	, $y_legend_x_axis+24	, "\" style=\" font-size:10;font-family: Arial; stroke-width:1; fill:", $linecolor[1], "\" >  </text>" ;
		################################################### END Legend Print


		##############################
		# print svg header
		print OUTsvg	"\n\t<text x=\"", $x_value_header ,"\" y=\"", $y_value_header, "\" style=\"font-size:"	, $font_size, ";font-family: ", $text_font_family, "; stroke-width:", $text_stroke_width	,";\" > ", $bartitle ," </text>" ;
		##############################


		################################################### Start Title Print
		# Print title
		print OUTsvg	"\n\t<text x=\"", $x_value ,"\" y=\"", $y_value_header-24 ,"\" style=\"font-size:12;font-family: Arial; stroke-width:1; fill:black\" > ", $topology ," </text>" ;
		################################################### End Title Print


		$y_value			= $y_value	+ 400;																		# start point y value
		$axis_y_start		= $y_value - 20	;																	# y start point y axis
		$length_y_axis		= $axis_y_start  - $high_bar_plot - $axis_space_end ;								# length y axis
		$y_value_header		= $axis_y_start  - $high_bar_plot - $axis_space_end - 20 ;							# y start point header
		$y_legend_x_axis	= $y_value ;																		# y start point legend x axis
		$y_legend_descript	= $axis_y_start  - $high_bar_plot ;													# y start point legend description
		$y_legend_y_axis	= $length_y_axis - 10 ;																# y start point legend description
		$bar_plot_start_x	= $axis_x_start  + $width_bar_space ;												# x startpoint first bar plot
	}

	##############################
	# print out svg end line and close outfile
	print OUTsvg				"\n</svg>\n" ;
	close OUTsvg ;
	##############################


	##############################
	# Print terminal info
	print "\n\tPrint SVG ", $outfile ;
	##############################
}

sub print_txt_quartet_results{

	my	$href_hol_taxvalues_of_subclade	=	$_[0]	;	#	key1: subclade, key2: support code (0,1), value: list quartet taxa, split support value 	In -> Defined, OUT -> Unchanged
	my	$sref_outfile_split_txt			=	$_[1]	;	#	name of the txt split support outputfile													In -> Defined, OUT -> Unchanged
	my	$sref_name_folder				=	$_[2]	;	#	key: subfolder code; value: subfolder name													In -> Defined, OUT -> Unchanged
	my	$aref_quartets					=	$_[3]	;	#	list of best topologies, sorted by scoring													In -> Defined, OUT -> Unchanged

	######################################
	# Print OUT result txt info about resolved quartets
	my	$outfile		=	$$sref_outfile_split_txt	;

	print "\n\tPrint TXT ", $outfile ;

	open	OUTtxt, ">$$sref_name_folder/$outfile" or die "\n\tOUTFILE-ERROR: Cannot open result file ", $$sref_name_folder, "/", $outfile, " in subroutine &print_txt_quartet_results!\n" ;


	for my $subclade ( @$aref_quartets ){ #print "\nsubclade: ", $subclade;

		print OUTtxt	"Quartet topology ", $subclade ,":\n",
						"--------------------------------------------------------------------------\n",
						"Best Split Score\t\t2nd Best Split Score\t\tLowest Split Score\n"
		;


		my	@quartet_value_0	=	exists ( $href_hol_taxvalues_of_subclade->{$subclade}{0} ) ? @{$href_hol_taxvalues_of_subclade->{$subclade}{0}} : () ;
		my	@quartet_value_1	=	exists ( $href_hol_taxvalues_of_subclade->{$subclade}{1} ) ? @{$href_hol_taxvalues_of_subclade->{$subclade}{1}} : () ;
		my	@quartet_value_2	=	exists ( $href_hol_taxvalues_of_subclade->{$subclade}{2} ) ? @{$href_hol_taxvalues_of_subclade->{$subclade}{2}} : () ;


		my $N_prints ;
		if		( ( @quartet_value_0 >= @quartet_value_1 ) && ( @quartet_value_0 >= @quartet_value_2 ) ){ $N_prints = ( @quartet_value_0 / 2 ) }
		elsif	( ( @quartet_value_1 >= @quartet_value_2 ) && ( @quartet_value_1 >= @quartet_value_0 ) ){ $N_prints = ( @quartet_value_1 / 2 ) }
		elsif	( ( @quartet_value_2 >= @quartet_value_1 ) && ( @quartet_value_2 >= @quartet_value_0 ) ){ $N_prints = ( @quartet_value_2 / 2 ) }
		else	{ die "\n\tBUG-ERROR: Cannot assign split values in subroutine &start_phyquart!\n\tPlease, report BUG to system developer!\n\t" }

		for my $number ( 0 .. $N_prints-1 ){

			if ( @quartet_value_0 ){ my $q = shift @quartet_value_0 ; my $v = shift @quartet_value_0 ; print OUTtxt $q, "\t", $v, "\t"	} else { print OUTtxt "\t\t"	}
			if ( @quartet_value_1 ){ my $q = shift @quartet_value_1 ; my $v = shift @quartet_value_1 ; print OUTtxt $q, "\t", $v, "\t"	} else { print OUTtxt "\t\t"	}
			if ( @quartet_value_2 ){ my $q = shift @quartet_value_2 ; my $v = shift @quartet_value_2 ; print OUTtxt $q, "\t", $v, "\n"	} else { print OUTtxt "\t\n"	}
		}

		print OUTtxt "\n"
	}

	close OUTtxt
	######################################
}
