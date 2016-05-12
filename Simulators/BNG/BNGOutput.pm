package BNGModel;
# BNGOutput is part of the BNGModel package. This file contains Output commands
# including:  writeXML, writeSBML, writeSSC, writeMfile, writeMexfile, ...


# pragmas
use strict;
use warnings;



###
###
###


# write Model in XML format
# $err = $model->writeXML({opt=>val,..}) 
sub writeXML
{
    use strict;
    use warnings;

    my $model       = shift @_;
	my $user_params = @_ ? shift @_ : {};

	my %params = (
        'evaluate_expressions' => 1,
		'format'               => 'xml',
        'include_model'        => 1,
        'include_network'      => 0,
        'overwrite'            => 1,
	);

    # copy user_params into pass_params structures
	while ( my ($key,$val) = each %$user_params )
	{   $params{$key} = $val;	}

    # writeFile will generate the output
    return $model->writeFile( \%params );
}


# generate XML string representing the BNGL model
sub toXML
{
	my $model = shift @_;
	my $user_params = @_ ? shift @_ : {};

    # default parameters
    my %params = (
        'evaluate_expressions' => 1,
    );

    # add user parameters
    while ( my ($key,$val) = each %$user_params )
    {   $params{$key} = $val;   }

	return '' if $BNGModel::NO_EXEC;

    # get BNG version	
	my $version = BNGversion();

    # get mopdel name
	my $model_name = $model->Name;

    # define size of indent
	my $indent = "    ";

    # are we evaluating expressions?
    my $evaluate_expressions = $params{'evaluate_expressions'};

    # Begin writing XML #
	# HEADER
	my $xml =  "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
              ."<!-- Created by BioNetGen $version  -->\n"
              ."<sbml xmlns=\"http://www.sbml.org/sbml/level3\" level=\"3\" version=\"1\">\n"
              ."  <model id=\"$model_name\">\n";


	# Parameters
	$xml .= $indent . "<ListOfParameters>\n";
	my $indent2 = "  " . $indent;
	my $plist   = $model->ParamList;
	foreach my $param ( @{$plist->Array} )
    {
		my $value;
		my $type;
		my $do_print = 0;
		if ( $param->Type =~ /^Constant/ )
		{
			$value = ($evaluate_expressions) ? sprintf "%.8g", $param->evaluate([], $plist) : $param->toString($plist);
			$type  = ($evaluate_expressions) ? "Constant" : $param->Type;
			$do_print = 1;
		}
		next unless $do_print;
		$xml .= sprintf( "$indent2<Parameter id=\"%s\"", $param->Name );
		$xml .= " type=\"$type\"";
		$xml .= " value=\"$value\"";
		$xml .= "/>\n";
	}
	$xml .= $indent . "</ListOfParameters>\n";

	# Molecule Types
	$xml .= $model->MoleculeTypesList->toXML($indent);

	# Compartments
	$xml .= $model->CompartmentList->toXML($indent);

	# Species
	if (@{$model->Concentrations}){
	    $xml .= $model->SpeciesList->toXML($indent,$model->Concentrations);
	} else {
	    $xml .= $model->SpeciesList->toXML($indent);
	}

	# Reaction rules
	my $string = $indent . "<ListOfReactionRules>\n";
	$indent2 = "  " . $indent;
	my $rindex  = 1;
	foreach my $rset ( @{$model->RxnRules} )
    {
		foreach my $rr ( @$rset )
        {
			$string .= $rr->toXML( $indent2, $rindex, $plist );
			++$rindex;
		}
	}
	$string .= $indent . "</ListOfReactionRules>\n";
	$xml .= $string;

	# Observables
	$string  = $indent . "<ListOfObservables>\n";
	$indent2 = "  " . $indent;
	my $oindex  = 1;
	foreach my $obs ( @{$model->Observables} )
    {
		$string .= $obs->toXML( $indent2, $oindex );
		++$oindex;
	}
	$string .= $indent . "</ListOfObservables>\n";
	$xml .= $string;

	# Functions
	$xml .= $indent . "<ListOfFunctions>\n";
	$indent2 = "  " . $indent;
	foreach my $param ( @{$plist->Array} )
    {
		next unless ( $param->Type eq "Function" );
		$xml .= $param->Ref->toXML( $plist, $indent2 );
	}
	$xml .= $indent . "</ListOfFunctions>\n";

	# FOOTER
	$xml .=  "  </model>\n"
            ."</sbml>\n";

	return $xml;
}



###
###
###



# write reaction network to SBML level 2 format
sub writeSBML
{
	my $model  = shift @_;
	my $params = @_ ? shift @_ : {};

	return '' if $BNGModel::NO_EXEC;

    # nothing to do unless a reactions are defined
	unless ( defined $model->RxnList  and  @{$model->RxnList->Array} )
	{   return "writeSBML(): No reactions in current model--nothing to do.";   }

    # get parameter list
	my $plist = $model->ParamList;

    # get model name
	my $model_name = $model->Name;

	# Strip prefixed path
	my $prefix = ( defined $params->{prefix} ) ? $params->{prefix} : $model->getOutputPrefix();
	my $suffix = ( defined $params->{suffix} ) ? $params->{suffix} : '';
	unless ( $suffix eq '' )
    {   $prefix .= "_${suffix}";   }

    # define file name
	my $file = "${prefix}.xml";

    # open file
    my $SBML;
	open( $SBML, '>', $file )  or die "Couldn't open $file: $!\n";

    # get BNG version
	my $version = BNGversion();


	# 0. HEADER
	print $SBML <<"EOF";
<?xml version="1.0" encoding="UTF-8"?>
<!-- Created by BioNetGen $version  -->
<sbml xmlns="http://www.sbml.org/sbml/level2" level="2" version="1">
  <model id="$model_name">
EOF


	# 1. Compartments (currently one dimensionsless compartment)
	print $SBML <<"EOF";
    <listOfCompartments>
      <compartment id="cell" size="1"/>
    </listOfCompartments>
EOF


	# 2. Species
	print $SBML "    <listOfSpecies>\n";

	my $use_array = @{$model->Concentrations} ? 1 : 0;
	foreach my $spec ( @{$model->SpeciesList->Array} )
	{
		my $conc;
		if ($use_array)
        {   $conc = $model->Concentrations->[ $spec->Index - 1 ];   }
		else
        {   $conc = $spec->Concentration;   }

        # If concentration is a parameter name, then evaluate the parameter
		unless ( isReal($conc) )
        {   $conc = $plist->evaluate($conc, []);   }

		printf $SBML "      <species id=\"S%d\" compartment=\"%s\" initialConcentration=\"%.8g\"",
		                                                                $spec->Index, "cell", $conc;

		if ( $spec->SpeciesGraph->Fixed )
        {   printf $SBML " boundaryCondition=\"true\"";   }

		printf $SBML " name=\"%s\"", $spec->SpeciesGraph->StringExact;
		print $SBML "/>\n";
	}
	print $SBML "    </listOfSpecies>\n";


	# 3. Parameters
	# A. Rate constants
	print $SBML "    <listOfParameters>\n";
	print $SBML "      <!-- Independent variables -->\n";
	foreach my $param ( @{$plist->Array} )
	{
	    next unless ( $param->Type eq 'Constant' );
		printf $SBML "      <parameter id=\"%s\" value=\"%.8g\"/>\n", $param->Name, $param->evaluate([], $plist);
	}
	print $SBML "      <!-- Dependent variables -->\n";
	foreach my $param ( @{$plist->Array} )
	{
	    next unless ( $param->Type eq 'ConstantExpression' );	
		printf $SBML "      <parameter id=\"%s\" constant=\"false\"/>\n", $param->Name;
	}

	# B. Observables
	if ( @{$model->Observables} )
	{
		print $SBML "      <!-- Observables -->\n";
	}
	foreach my $obs ( @{$model->Observables} )
	{
		printf $SBML "      <parameter id=\"%s\" constant=\"false\"/>\n", "Group_" . $obs->Name;
	}
	print $SBML "    </listOfParameters>\n";


	# 4. 'Rules' (for observables)
	print $SBML "    <listOfRules>\n";
	print $SBML "      <!-- Dependent variables -->\n";
	foreach my $param ( @{$plist->Array} )
    {
		next if ( $param->Expr->Type eq 'NUM' );
		printf $SBML "      <assignmentRule variable=\"%s\">\n", $param->Name;
        #print  $SBML "        <notes>\n";
        #print  $SBML "          <xhtml:p>\n";
        #printf $SBML "            %s=%s\n", $param->Name,$param->toString($plist);
        #print  $SBML "          </xhtml:p>\n";
        #print  $SBML "        </notes>\n";
		printf $SBML $param->toMathMLString( $plist, "        " );
		print $SBML "      </assignmentRule>\n";
	}
	if ( @{$model->Observables} )
    {
		print $SBML "      <!-- Observables -->\n";
		foreach my $obs ( @{$model->Observables} )
        {
			printf $SBML "      <assignmentRule variable=\"%s\">\n",
			  "Group_" . $obs->Name;
			my ( $ostring, $err ) = $obs->toMathMLString();
			if ($err) { return $err; }
			foreach my $line ( split "\n", $ostring )
            {
				print $SBML "          $line\n";
			}
			print $SBML "      </assignmentRule>\n";
		}
	}
	print $SBML "    </listOfRules>\n";


	# 5. Reactions
	print $SBML "    <listOfReactions>\n";
	my $index = 0;
	foreach my $rxn ( @{$model->RxnList->Array} )
    {
		++$index;
		printf $SBML "      <reaction id=\"R%d\" reversible=\"false\">\n", $index;

		#Get indices of reactants
		my @rindices = ();
		foreach my $spec ( @{$rxn->Reactants} )
        {
			push @rindices, $spec->Index;
		}
		@rindices = sort { $a <=> $b } @rindices;

		#Get indices of products
		my @pindices = ();
		foreach my $spec ( @{$rxn->Products} )
        {
			push @pindices, $spec->Index;
		}
		@pindices = sort { $a <=> $b } @pindices;

		print $SBML "        <listOfReactants>\n";
		foreach my $i (@rindices)
        {
			printf $SBML "          <speciesReference species=\"S%d\"/>\n", $i;
		}
		print $SBML "        </listOfReactants>\n";

		print $SBML "        <listOfProducts>\n";
		foreach my $i (@pindices)
        {
			printf $SBML "          <speciesReference species=\"S%d\"/>\n", $i;
		}
		print $SBML "        </listOfProducts>\n";

		print $SBML "        <kineticLaw>\n";
		my ( $rstring, $err ) = $rxn->RateLaw->toMathMLString( \@rindices, \@pindices, $rxn->StatFactor );
		if ($err) { return $err; }

		foreach my $line ( split "\n", $rstring )
        {
			print $SBML "          $line\n";
		}
		print $SBML "        </kineticLaw>\n";

		print $SBML "      </reaction>\n";
	}
	print $SBML "    </listOfReactions>\n";

	# 6. FOOTER
	print $SBML <<"EOF";
  </model>
</sbml>
EOF

    close $SBML;
	print "Wrote SBML to $file.\n";
	return;
}



###
###
###



sub writeSSC
{
	my $model  = shift @_;
	my $params = @_ ? shift @_ : {};
	return '' if $BNGModel::NO_EXEC;

    # get model name
	my $model_name = $model->Name;

	my $prefix = ( defined $params->{prefix} ) ? $params->{prefix} : $model->getOutputPrefix();
	my $suffix = ( defined $params->{suffix} ) ? $params->{suffix} : undef;
	if ( $suffix ) {
		$prefix .= "_${suffix}";
	}
	my $file = "${prefix}.rxn";
	open( SSCfile, ">$file" ) || die "Couldn't open $file: $!\n";
	my $version = BNGversion();
	print SSCfile
	  "--# SSC-file for model $model_name created by BioNetGen $version\n";
	print "Writing SSC translator .rxn file.....";

	#-- Compartment default for SSC ---- look more into it
	printf SSCfile
	  "region World \n  box width 1 height 1 depth 1\nsubvolume edge 1";

	# --This part correspond to seed specie
	print SSCfile "\n\n";
	print SSCfile "--# Initial molecules and their concentrations\n";
	my $sp_string = $model->SpeciesList->writeSSC( $model->Concentrations,
                                                   $model->ParamList       );
	print SSCfile $sp_string;

	# --This part in SSC corrsponds to Observables
	if ( @{$model->Observables })
    {
		print SSCfile"\n\n--# reads observables";
		print SSCfile "\n";
		foreach my $obs ( @{$model->Observables} )
        {
			my $ob_string = $obs->toStringSSC();
			if ( $ob_string =~ /\?/ )
            {
				print STDOUT " \n WARNING: SSC does not implement ?. The observable has been commented. Please see .rxn file for more details \n";
				print STDOUT "\n See Observable\n", $obs->toString();
				$ob_string = "\n" . "--#" . "record " . $ob_string;
				print SSCfile $ob_string;
			}    #putting this string as a comment and carrying on
			else
            {
				print SSCfile "\nrecord ", $ob_string;
			}
		}
	}

	# --Reaction rules
	print SSCfile" \n\n--# reaction rules\n";
	foreach my $rset ( @{$model->RxnRules} )
    {
		my $id = 0;
		my $rreverse = ( $#$rset > 0 ) ? $rset->[1] : "";
		( my $reac1, my $errorSSC ) = $rset->[0]->toStringSSC($rreverse);
		if ( $errorSSC == 1 )
        {
			print STDOUT "\nSee rule in .rxn \n",
			  $rset->[0]->toString($rreverse);
			$reac1 = "--#" . $reac1;
		}
		print SSCfile $reac1;
		print SSCfile "\n";
		if ($rreverse)
        {
			( my $reac2, my $errorSSC ) = $rset->[1]->toStringSSC($rreverse);
			if ( $errorSSC == 1 ) { $reac2 = "--#" . $reac2; }
			print SSCfile $reac2;
			print SSCfile "\n";
		}
	}
	print "\nWritten SSC file\n";
	return ();
}



# This subroutine writes a file which contains the information corresponding to the parameter block in BNG
sub writeSSCcfg
{
	my $model  = shift @_;
	my $params = @_ ? shift @_ : {};
	return '' if $BNGModel::NO_EXEC;

    # get model name
	my $model_name = $model->Name;

	# Strip prefixed path
	my $prefix = ( defined $params->{prefix} ) ? $params->{prefix} : $model->getOutputPrefix();
	my $suffix = ( defined $params->{suffix} ) ? $params->{suffix} : '';
	if ( $suffix ne '' )
    {
		$prefix .= "_${suffix}";
	}

	my $file    = "${prefix}.cfg";
	my $version = BNGversion();

	open( SSCcfgfile, ">$file" ) || die "Couldn't open $file: $!\n";
	print STDOUT "\n Writting SSC cfg file \n";
	print SSCcfgfile "# SSC cfg file for model $model_name created by BioNetGen $version\n";
	print SSCcfgfile $model->ParamList->writeSSCcfg();

	return;
}



###
###
###



# Write model to a MATLAB M-file
sub writeMfile
{	
	my $model = shift @_;
	my $params = @_ ? shift @_ : {};

    # a place to hold errors
    my $err;

    # nothing to do if NO_EXEC is true
	return '' if $BNGModel::NO_EXEC;

    # nothing to do if there are no reactions
	unless ( $model->RxnList )
	{
	    return ( "writeMfile() has nothing to do: no reactions in current model.\n"
	            ."  Did you remember to call generate_network() before attempting to\n"
	            ."  write network output?");
	}

    # get reference to parameter list
	my $plist = $model->ParamList;

	
	# get model name
	my $model_name = $model->Name;

    # Build output file name
	# ..use prefix if defined, otherwise use model name
	my $prefix = ( defined $params->{prefix} ) ? $params->{prefix} : $model->getOutputPrefix();
	# ..add suffix, if any
	my $suffix = ( defined $params->{suffix} ) ? $params->{suffix} : undef;
	if ( $suffix )
	{   $prefix .= "_${suffix}";   }

    # split prefix into volume, path and filebase
    my ($vol, $path, $filebase) = File::Spec->splitpath($prefix);
    
	# define m-script file name
    my $mscript_filebase = "${filebase}";
    my $mscript_filename = "${mscript_filebase}.m";
	my $mscript_path     = File::Spec->catpath($vol,$path,$mscript_filename);
    my $mscript_filebase_caps = uc $mscript_filebase;

    # configure options (see Matlab documentation of functions ODESET and ODE15S)
    my $odeset_abstol = 1e-4;
    if ( exists $params->{'atol'} )
    {   $odeset_abstol = $params->{'atol'};  }
    
    my $odeset_reltol = 1e-8;
    if ( exists $params->{'rtol'} )
    {   $odeset_reltol = $params->{'rtol'};  } 

    my $odeset_stats = 'off';
    if ( exists $params->{'stats'} )
    {   $odeset_stats = $params->{'stats'};  } 

    my $odeset_bdf = 'off';
    if ( exists $params->{'bdf'} )
    {   $odeset_bdf = $params->{'bdf'};  }

    my $odeset_maxorder = 5;
    if ( exists $params->{'maxOrder'} )
    {   $odeset_maxorder = $params->{'maxOrder'};  } 

    # time options for mscript
    my $t_start = 0;
    if ( exists $params->{'t_start'} )
    {   $t_start = $params->{'t_start'};  }  

    my $t_end = 10;
    if ( exists $params->{'t_end'} )
    {   $t_end = $params->{'t_end'};  } 

    my $n_steps = 20;
    if ( exists $params->{'n_steps'} )
    {   $n_steps = $params->{'n_steps'};  } 

    # configure time step dependent options
    my $odeset_maxstep = undef;
    if ( exists $params->{'max_step'} )
    {   $odeset_maxstep = $params->{'max_step'};  }     
    
    # construct ODESET function call
    my $mscript_call_odeset;
    if ( defined $odeset_maxstep )
    {
        $mscript_call_odeset = "opts = odeset( 'RelTol',   $odeset_reltol,   ...\n"
                              ."               'AbsTol',   $odeset_abstol,   ...\n"
                              ."               'Stats',    '$odeset_stats',  ...\n"
                              ."               'BDF',      '$odeset_bdf',    ...\n"
                              ."               'MaxOrder', $odeset_maxorder, ...\n"
                              ."               'MaxStep',  $odeset_maxstep    );\n";
    }
    else
    {
        $mscript_call_odeset = "opts = odeset( 'RelTol',   $odeset_reltol,   ...\n"
                              ."               'AbsTol',   $odeset_abstol,   ...\n"
                              ."               'Stats',    '$odeset_stats',  ...\n"
                              ."               'BDF',      '$odeset_bdf',    ...\n"
                              ."               'MaxOrder', $odeset_maxorder   );\n";    
    }

    # Index parameters associated with Constants, ConstantExpressions and Observables
    ($err) = $plist->indexParams();
    if ($err) { return $err };

    # and retrieve a string of expression definitions
    my $n_parameters = $plist->countType( 'Constant' );
    my $n_expressions = $plist->countType( 'ConstantExpression' ) + $n_parameters;
    (my $calc_expressions_string, $err) = $plist->getMatlabExpressionDefs();    
    if ($err) { return $err };

    # get list of parameter names and defintions for matlab
	my $mscript_param_names;
	my $mscript_param_values;
	($mscript_param_names, $mscript_param_values, $err) = $plist->getMatlabConstantNames();
    if ($err) { return $err };

    # get number of species
    my $n_species = scalar @{$model->SpeciesList->Array};
     
	# retrieve a string of observable definitions
    my $n_observables = scalar @{$model->Observables};
    my $calc_observables_string;
    ($calc_observables_string, $err) = $plist->getMatlabObservableDefs();
    if ($err) { return $err };
    
    # get list of observable names for matlab
	my $mscript_observable_names;
	($mscript_observable_names, $err) = $plist->getMatlabObservableNames();
    if ($err) { return $err };
    
    # Construct user-defined functions
    my $user_fcn_declarations = '';
    my $user_fcn_definitions  = '';
	foreach my $param ( @{ $model->ParamList->Array } )
	{
		if ( $param->Type eq 'Function' )
		{
		    # get reference to the actual Function
		    my $fcn = $param->Ref;
		    
		    # don't write function if it depends on a local observable evaluation (this is useless
		    #   since CVode can't do local evaluations)
		    next if ( $fcn->checkLocalDependency($plist) );
		    		    
		    # get function definition			    
		    my $fcn_defn = $fcn->toMatlabString( $plist, {fcn_mode=>'define', indent=>''} );

		    # add definition to the user_fcn_definitions string
		    $user_fcn_definitions .= $fcn_defn . "\n";
        }
	}
	
    # index reactions
    ($err) = $model->RxnList->updateIndex( $plist );
    if ($err) { return $err };

	# retrieve a string of reaction rate definitions
	my $n_reactions = scalar @{$model->RxnList->Array};
    my $calc_ratelaws_string;
    ($calc_ratelaws_string, $err) = $model->RxnList->getMatlabRateDefs( $plist );
    if ($err) { return $err };
    

    # get stoichiometry matrix (sparse encoding in a hashmap)
	my $stoich_hash = {};
	($err) = $model->RxnList->calcStoichMatrix( $stoich_hash );

	# retrieve a string of species deriv definitions
    my $calc_derivs_string;
    ($calc_derivs_string, $err) = $model->SpeciesList->toMatlabString( $model->RxnList, $stoich_hash, $plist );
    if ($err) { return $err };   	


    # get list of species names and initial value expressions for matlab
	my $mscript_species_names;
	my $mscript_species_init;
	($mscript_species_names, $mscript_species_init, $err) = $model->SpeciesList->getMatlabSpeciesNames( $model );
    if ($err) { return $err }; 


    ## Set up MATLAB Plot
    # fontsizes
    my $title_fontsize = 14;
    my $axislabel_fontsize = 12;
    my $legend_fontsize = 10;

    # generate code snippets for plotting observables or species
    my $mscript_plot_labels;
    my $mscript_make_plot;

    # get ylabel (either Number of Concentration)
    my $ylabel;
    if ( $model->SubstanceUnits eq 'Number' )
    {   $ylabel = 'number';   }
    elsif ( $model->SubstanceUnits eq 'Concentration' )
    {   $ylabel = 'concentration';   }
    else
    {   return "writeMfile(): I could not identify model substance units!";   }

    
    if ( @{$model->Observables} )
    {   # plot observables
        $mscript_plot_labels = "    observable_labels = { $mscript_observable_names };\n";
        
        $mscript_make_plot = "    plot(timepoints,observables_out);\n"
                            ."    title('${mscript_filebase} observables','fontSize',${title_fontsize},'Interpreter','none');\n"
                            ."    axis([${t_start} timepoints(end) 0 inf]);\n"
                            ."    legend(observable_labels,'fontSize',${legend_fontsize},'Interpreter','none');\n"
                            ."    xlabel('time','fontSize',${axislabel_fontsize},'Interpreter','none');\n"
                            ."    ylabel('${ylabel}','fontSize',${axislabel_fontsize},'Interpreter','none');\n";
    
    }
    else
    {   # plot species
        $mscript_plot_labels = "    species_labels = { $mscript_species_names };\n";
    
        $mscript_make_plot = "    plot(timepoints,species_out);\n"
                            ."    title('${mscript_filebase} species','fontSize',${title_fontsize},'Interpreter','none');\n"
                            ."    axis([${t_start} timepoints(end) 0 inf]);\n"
                            ."    legend(species_labels,'fontSize',${legend_fontsize},'Interpreter','none');\n"
                            ."    xlabel('time','fontSize',${axislabel_fontsize},'Interpreter','none');\n"
                            ."    ylabel('${ylabel}','fontSize',${axislabel_fontsize},'Interpreter','none');\n";
    }
    


    # open Mexfile and begin printing...
	open( Mscript, ">$mscript_path" ) || die "Couldn't open $mscript_path: $!\n";
    print Mscript <<"EOF";
function [err, timepoints, species_out, observables_out ] = ${mscript_filebase}( timepoints, species_init, parameters, suppress_plot )
%${mscript_filebase_caps} Integrate reaction network and plot observables.
%   Integrates the reaction network corresponding to the BioNetGen model
%   '${model_name}' and then (optionally) plots the observable trajectories,
%   or species trajectories if no observables are defined. Trajectories are
%   generated using either default or user-defined parameters and initial
%   species values. Integration is performed by the MATLAB stiff solver
%   'ode15s'. ${mscript_filebase_caps} returns an error value, a vector of timepoints,
%   species trajectories, and observable trajectories.
%   
%   [err, timepoints, species_out, observables_out]
%        = $mscript_filebase( timepoints, species_init, parameters, suppress_plot )
%
%   INPUTS:
%   -------
%   species_init    : row vector of $n_species initial species populations.
%   timepoints      : column vector of time points returned by integrator.
%   parameters      : row vector of $n_parameters model parameters.
%   suppress_plot   : 0 if a plot is desired (default), 1 if plot is suppressed.
%
%   Note: to specify default value for an input argument, pass the empty array.
%
%   OUTPUTS:
%   --------
%   err             : 0 if the integrator exits without error, non-zero otherwise.
%   timepoints      : a row vector of timepoints returned by the integrator.
%   species_out     : array of species population trajectories
%                        (columns correspond to species, rows correspond to time).
%   observables_out : array of observable trajectories
%                        (columns correspond to observables, rows correspond to time).
%
%   QUESTIONS about the BNG Mfile generator?  Email justinshogg\@gmail.com



%% Process input arguments

% define any missing arguments
if ( nargin < 1 )
    timepoints = [];
end

if ( nargin < 2 )
    species_init = [];
end

if ( nargin < 3 )
    parameters = [];
end

if ( nargin < 4 )
    suppress_plot = 0;
end


% initialize outputs (to avoid error msgs if script terminates early
err = 0;
species_out     = [];
observables_out = [];


% setup default parameters, if necessary
if ( isempty(parameters) )
   parameters = [ $mscript_param_values ];
end
% check that parameters has proper dimensions
if (  size(parameters,1) ~= 1  |  size(parameters,2) ~= $n_parameters  )
    fprintf( 1, 'Error: size of parameter argument is invalid! Correct size = [1 $n_parameters].\\n' );
    err = 1;
    return;
end

% setup default initial values, if necessary
if ( isempty(species_init) )
   species_init = initialize_species( parameters );
end
% check that species_init has proper dimensions
if (  size(species_init,1) ~= 1  |  size(species_init,2) ~= $n_species  )
    fprintf( 1, 'Error: size of species_init argument is invalid! Correct size = [1 $n_species].\\n' );
    err = 1;
    return;
end

% setup default timepoints, if necessary
if ( isempty(timepoints) )
   timepoints = linspace($t_start,$t_end,$n_steps+1)';
end
% check that timepoints has proper dimensions
if (  size(timepoints,1) < 2  |  size(timepoints,2) ~= 1  )
    fprintf( 1, 'Error: size of timepoints argument is invalid! Correct size = [t 1], t>1.\\n' );
    err = 1;
    return;
end

% setup default suppress_plot, if necessary
if ( isempty(suppress_plot) )
   suppress_plot = 0;
end
% check that suppress_plot has proper dimensions
if ( size(suppress_plot,1) ~= 1  |  size(suppress_plot,2) ~= 1 )
    fprintf( 1, 'Error: suppress_plots argument should be a scalar!\\n' );
    err = 1;
    return;
end

% define parameter labels (this is for the user's reference!)
param_labels = { $mscript_param_names };



%% Integrate Network Model
 
% calculate expressions
[expressions] = calc_expressions( parameters );

% set ODE integrator options
$mscript_call_odeset

% define derivative function
rhs_fcn = @(t,y)( calc_species_deriv( t, y, expressions ) );

% simulate model system (stiff integrator)
try 
    [timepoints, species_out] = ode15s( rhs_fcn, timepoints, species_init', opts );
catch
    err = 1;
    fprintf( 1, 'Error: some problem encounteredwhile integrating ODE network!\\n' );
    return;
end

% calculate observables
observables_out = zeros( length(timepoints), $n_observables );
for t = 1 : length(timepoints)
    observables_out(t,:) = calc_observables( species_out(t,:), expressions );
end


%% Plot Output, if desired

if ( ~suppress_plot )
    
    % define plot labels
$mscript_plot_labels
    % construct figure
$mscript_make_plot
end


%~~~~~~~~~~~~~~~~~~~~~%
% END of main script! %
%~~~~~~~~~~~~~~~~~~~~~%


% initialize species function
function [species_init] = initialize_species( params )

    species_init = zeros(1,$n_species);
$mscript_species_init
end


% user-defined functions
$user_fcn_definitions


% Calculate expressions
function [ expressions ] = calc_expressions ( parameters )

    expressions = zeros(1,$n_expressions);
$calc_expressions_string   
end



% Calculate observables
function [ observables ] = calc_observables ( species, expressions )

    observables = zeros(1,$n_observables);
$calc_observables_string
end


% Calculate ratelaws
function [ ratelaws ] = calc_ratelaws ( species, expressions, observables )

    ratelaws = zeros(1,$n_observables);
$calc_ratelaws_string
end

% Calculate species derivates
function [ Dspecies ] = calc_species_deriv ( time, species, expressions )
    
    % initialize derivative vector
    Dspecies = zeros($n_species,1);
    
    % update observables
    [ observables ] = calc_observables( species, expressions );
    
    % update ratelaws
    [ ratelaws ] = calc_ratelaws( species, expressions, observables );
                        
    % calculate derivatives
$calc_derivs_string
end


end
EOF

	close(Mscript);
	print "Wrote M-file script $mscript_path.\n";
	return ();	
}



###
###
###






sub writeMexfile
{
	my $model = shift;
	my $params = (@_) ? shift : {};

    # a place to hold errors
    my $err;

    # nothing to do if NO_EXEC is true
	return ('') if $BNGModel::NO_EXEC;

    # nothing to do if there are no reactions
	unless ( $model->RxnList )
	{
	    return ( "writeMexfile() has nothing to do: no reactions in current model.\n"
	            ."  Did you remember to call generate_network() before attempting to\n"
	            ."  write network output?");
	}

    # get reference to parameter list
	my $plist = $model->ParamList;
	
	# get model name
	my $model_name = $model->Name;
    
	# Strip prefixed path
	my $prefix = ( defined $params->{prefix} ) ? $params->{prefix} : $model->getOutputPrefix();
	my $suffix = ( defined $params->{suffix} ) ? $params->{suffix} : undef;
	if ( $suffix )
	{   $prefix .= "_${suffix}";   }
	
    # split prefix into volume, path and filebase
    my ($vol, $path, $filebase) = File::Spec->splitpath($prefix);

	# define mexfile name
    my $mex_filebase = "${filebase}_cvode";
	my $mex_filename = "${mex_filebase}.c";
    my $mex_path     = File::Spec->catpath($vol,$path,$mex_filename);
    
	# define m-script files name
    my $mscript_filebase = "${filebase}";
    my $mscript_filename = "${mscript_filebase}.m";
	my $mscript_path     = File::Spec->catpath($vol,$path,$mscript_filename);
    my $mscript_filebase_caps = uc $mscript_filebase;

    # configure options
    my $cvode_abstol = 1e-6;
    if ( exists $params->{'atol'} )
    {   $cvode_abstol = $params->{'atol'};  }
    
    my $cvode_reltol = 1e-8;
    if ( exists $params->{'rtol'} )
    {   $cvode_reltol = $params->{'rtol'};  }    

    my $cvode_max_num_steps = 2000;
    if ( exists $params->{'max_num_steps'} )
    {   $cvode_max_num_steps = $params->{'max_num_steps'};  }  

    my $cvode_max_err_test_fails = 7;
    if ( exists $params->{'max_err_test_fails'} )
    {   $cvode_max_err_test_fails = $params->{'max_err_test_fails'};  }  

    my $cvode_max_conv_fails = 10;
    if ( exists $params->{'max_conv_fails'} )
    {   $cvode_max_conv_fails = $params->{'max_conv_fails'};  }  

    my $cvode_max_step = '0.0';
    if ( exists $params->{'max_step'} )
    {   $cvode_max_step = $params->{'max_step'};  }

    # Stiff = CV_BDF,CV_NEWTON (Default); Non-stiff = CV_ADAMS,CV_FUNCTIONAL
    my $cvode_linear_multistep = 'CV_BDF';
    my $cvode_nonlinear_solver = 'CV_NEWTON';
    if ( exists $params->{'stiff'} )
    {   
        # if stiff is FALSE, then change to CV_ADAMS and CV_FUNCTIONAL
        unless ( $params->{'stiff'} )
        {
            $cvode_linear_multistep = 'CV_ADAMS';    
            $cvode_nonlinear_solver = 'CV_FUNCTIONAL';
        }
    }

    # set sparse option (only permitted with CV_NEWTON)
    my $cvode_linear_solver;
    if ( ($cvode_nonlinear_solver eq 'CV_NEWTON')  and  ($params->{'sparse'}) )
    {
        $cvode_linear_solver =     "flag = CVSpgmr(cvode_mem, PREC_NONE, 0);\n"
                              ."    if (check_flag(&flag, \"CVSpgmr\", 1))";
    }
    else
    {
        $cvode_linear_solver =     "flag = CVDense(cvode_mem, __N_SPECIES__);\n"
                              ."    if (check_flag(&flag, \"CVDense\", 1))";
    }

    # time options for mscript
    my $t_start = 0;
    if ( exists $params->{'t_start'} )
    {   $t_start = $params->{'t_start'};  }  

    my $t_end = 10;
    if ( exists $params->{'t_end'} )
    {   $t_end = $params->{'t_end'};  } 

    my $n_steps = 20;
    if ( exists $params->{'n_steps'} )
    {   $n_steps = $params->{'n_steps'};  } 

    # code snippet for cleaning up dynamic memory before exiting CVODE-MEX
    my $cvode_cleanup_memory =     "{                                  \n"
                              ."        N_VDestroy_Serial(expressions);\n"
                              ."        N_VDestroy_Serial(observables);\n"
                              ."        N_VDestroy_Serial(ratelaws);   \n"
                              ."        N_VDestroy_Serial(species);    \n"
                              ."        CVodeFree(&cvode_mem);         \n"
                              ."        return_status[0] = 1;          \n"
                              ."        return;                        \n"
                              ."    }                                  ";

    # Index parameters associated with Constants, ConstantExpressions and Observables
    ($err) = $plist->indexParams();
    if ($err) { return $err };

    # and retrieve a string of expression definitions
    my $n_parameters = $plist->countType( 'Constant' );
    my $n_expressions = $plist->countType( 'ConstantExpression' ) + $n_parameters;
    (my $calc_expressions_string, $err) = $plist->getCVodeExpressionDefs();    
    if ($err) { return $err };

    # get list of parameter names and defintions for matlab
	my $mscript_param_names;
	my $mscript_param_values;
	($mscript_param_names, $mscript_param_values, $err) = $plist->getMatlabConstantNames();
    if ($err) { return $err };


    # generate CVode references for species
    # (Do this now, because we need references to CVode species for Observable definitions and Rxn Rates)
    my $n_species = scalar @{$model->SpeciesList->Array};
    
     
	# retrieve a string of observable definitions
    my $n_observables = scalar @{$model->Observables};
    my $calc_observables_string;
    ($calc_observables_string, $err) = $plist->getCVodeObservableDefs();
    if ($err) { return $err };    
    
    # get list of observable names for matlab
	my $mscript_observable_names;
	($mscript_observable_names, $err) = $plist->getMatlabObservableNames();
    if ($err) { return $err };
    
    # Construct user-defined functions
    my $user_fcn_declarations = '';
    my $user_fcn_definitions = '';
	foreach my $param ( @{ $model->ParamList->Array } )
	{
		if ( $param->Type eq 'Function' )
		{
		    # get reference to the actual Function
		    my $fcn = $param->Ref;
		    
		    # don't write function if it depends on a local observable evaluation (this is useless
		    #   since CVode can't do local evaluations)
		    next if ( $fcn->checkLocalDependency($plist) );
		    		    
		    # get function declaration, add it to the user_fcn_declarations string
		    $user_fcn_declarations .= $fcn->toCVodeString( $plist, {fcn_mode=>'declare',indent=>''} );
		    
		    # get function definition			    
		    my $fcn_defn = $fcn->toCVodeString( $plist, {fcn_mode=>'define', indent=>''} );

		    # add definition to the user_fcn_definitions string
		    $user_fcn_definitions .= $fcn_defn . "\n";
        }
	}
	
    # index reactions
    ($err) = $model->RxnList->updateIndex( $plist );
    if ($err) { return $err };

	# retrieve a string of reaction rate definitions
	my $n_reactions = scalar @{$model->RxnList->Array};
    my $calc_ratelaws_string;
    ($calc_ratelaws_string, $err) = $model->RxnList->getCVodeRateDefs( $plist );
    if ($err) { return $err };
    

    # get stoichiometry matrix (sparse encoding in a hashmap)
	my $stoich_hash = {};
	($err) = $model->RxnList->calcStoichMatrix( $stoich_hash );

	# retrieve a string of species deriv definitions
    my $calc_derivs_string;
    ($calc_derivs_string, $err) = $model->SpeciesList->toCVodeString( $model->RxnList, $stoich_hash, $plist );
    if ($err) { return $err };   	



    # get list of species names and initial value expressions for matlab
	my $mscript_species_names;
	my $mscript_species_init;
	($mscript_species_names, $mscript_species_init, $err) = $model->SpeciesList->getMatlabSpeciesNames( $model );
    if ($err) { return $err };


    ## Set up MATLAB Plot
    # fontsizes
    my $title_fontsize = 14;
    my $axislabel_fontsize = 12;
    my $legend_fontsize = 10;

    # generate code snippets for plotting observables or species
    my $mscript_plot_labels;
    my $mscript_make_plot;

    # get ylabel (either Number of Concentration)
    my $ylabel;
    if ( $model->SubstanceUnits eq 'Number' )
    {   $ylabel = 'number';   }
    elsif ( $model->SubstanceUnits eq 'Concentration' )
    {   $ylabel = 'concentration';   }
    else
    {   return "writeMfile(): I could not identify model substance units!";   }

    
    if ( @{$model->Observables} )
    {   # plot observables
        $mscript_plot_labels = "    observable_labels = { $mscript_observable_names };\n";
        
        $mscript_make_plot = "    plot(timepoints,observables_out);\n"
                            ."    title('${mscript_filebase} observables','fontSize',${title_fontsize},'Interpreter','none');\n"
                            ."    axis([${t_start} timepoints(end) 0 inf]);\n"
                            ."    legend(observable_labels,'fontSize',${legend_fontsize},'Interpreter','none');\n"
                            ."    xlabel('time','fontSize',${axislabel_fontsize},'Interpreter','none');\n"
                            ."    ylabel('${ylabel}','fontSize',${axislabel_fontsize},'Interpreter','none');\n";
    
    }
    else
    {   # plot species
        $mscript_plot_labels = "    species_labels = { $mscript_species_names };\n";
    
        $mscript_make_plot = "    plot(timepoints,species_out);\n"
                            ."    title('${mscript_filebase} species','fontSize',${title_fontsize},'Interpreter','none');\n"
                            ."    axis([${t_start} timepoints(end) 0 inf]);\n"
                            ."    legend(species_labels,'fontSize',${legend_fontsize},'Interpreter','none');\n"
                            ."    xlabel('time','fontSize',${axislabel_fontsize},'Interpreter','none');\n"
                            ."    ylabel('${ylabel}','fontSize',${axislabel_fontsize},'Interpreter','none');\n";
    }


    # open Mexfile and begin printing...
	open( Mexfile, ">$mex_path" ) or die "Couldn't open $mex_path: $!\n";
    print Mexfile <<"EOF";
/*   
**   ${mex_filename}
**	 
**   Cvode-Mex implementation of BioNetGen model '$model_name'.
**
**   Code Adapted from templates provided by Mathworks and Sundials.
**   QUESTIONS about the code generator?  Email justinshogg\@gmail.com
**
**   Requires the CVODE libraries:  sundials_cvode and sundials_nvecserial.
**   https://computation.llnl.gov/casc/sundials/main.html
**
**-----------------------------------------------------------------------------
**
**   COMPILE in MATLAB:
**   mex -L<path_to_cvode_libraries> -I<path_to_cvode_includes>  ...
**          -lsundials_nvecserial -lsundials_cvode -lm ${mex_filename}
**
**   note1: if cvode is in your library path, you can omit path specifications.
**
**   note2: if linker complains about lib stdc++, try removing "-lstdc++"
**     from the mex configuration file "gccopts.sh".  This should be in the
**     matlab bin folder.
** 
**-----------------------------------------------------------------------------
**
**   EXECUTE in MATLAB:
**   [error_status, species_out, observables_out]
**        = ${mex_filebase}( timepoints, species_init, parameters )
**
**   timepoints      : column vector of time points returned by integrator.
**   parameters      : row vector of $n_parameters parameters.
**   species_init    : row vector of $n_species initial species populations.
**
**   error_status    : 0 if the integrator exits without error, non-zero otherwise.
**   species_out     : species population trajectories
**                        (columns correspond to states, rows correspond to time).
**   observables_out : observable trajectories
**                        (columns correspond to observables, rows correspond to time).
*/

/* Library headers */
#include "mex.h"
#include "matrix.h"
#include <stdlib.h>
#include <math.h>
#include <cvode/cvode.h>             /* prototypes for CVODE  */
#include <nvector/nvector_serial.h>  /* serial N_Vector       */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <cvode/cvode_spgmr.h>       /* prototype for CVSpgmr */

/* Problem Dimensions */
#define __N_PARAMETERS__   $n_parameters
#define __N_EXPRESSIONS__  $n_expressions
#define __N_OBSERVABLES__  $n_observables
#define __N_RATELAWS__     $n_reactions
#define __N_SPECIES__      $n_species

/* core function declarations */
void  mexFunction ( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] );
int   check_flag  ( void *flagvalue, char *funcname, int opt );
void  calc_expressions ( N_Vector expressions, double * parameters );
void  calc_observables ( N_Vector observables, N_Vector species, N_Vector expressions );
void  calc_ratelaws    ( N_Vector ratelaws,  N_Vector species, N_Vector expressions, N_Vector observables );
int   calc_species_deriv ( realtype time, N_Vector species, N_Vector Dspecies, void * f_data );

/* user-defined function declarations */
$user_fcn_declarations

/* user-defined function definitions  */
$user_fcn_definitions

/* Calculate expressions */
void
calc_expressions ( N_Vector expressions, double * parameters )
{
$calc_expressions_string   
}

/* Calculate observables */
void
calc_observables ( N_Vector observables, N_Vector species, N_Vector expressions )
{
$calc_observables_string
}

/* Calculate ratelaws */
void
calc_ratelaws ( N_Vector ratelaws, N_Vector species, N_Vector expressions, N_Vector observables )
{  
$calc_ratelaws_string
}


/* Calculate species derivates */
int
calc_species_deriv ( realtype time, N_Vector species, N_Vector Dspecies, void * f_data )
{
    int         return_val;
    N_Vector *  temp_data;
    
    N_Vector    expressions;
    N_Vector    observables;
    N_Vector    ratelaws;

    /* cast temp_data */
    temp_data = (N_Vector*)f_data;
     
    /* sget ratelaws Vector */
    expressions = temp_data[0];
    observables = temp_data[1];
    ratelaws    = temp_data[2];
       
    /* calculate observables */
    calc_observables( observables, species, expressions );
    
    /* calculate ratelaws */
    calc_ratelaws( ratelaws, species, expressions, observables );
                        
    /* calculate derivates */
$calc_derivs_string

    return(0);
}


/*
**   ========
**   main MEX
**   ========
*/
void mexFunction( int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[] )
{
    /* variables */
    double *  return_status;
    double *  species_out;
    double *  observables_out;    
    double *  parameters;
    double *  species_init;
    double *  timepoints; 
    size_t    n_timepoints;
    size_t    i;
    size_t    j;

    /* intermediate data vectors */
    N_Vector  expressions;
    N_Vector  observables;
    N_Vector  ratelaws;

    /* array to hold pointers to data vectors */
    N_Vector  temp_data[3];
    
    /* CVODE specific variables */
    realtype  reltol;
    realtype  abstol;
    realtype  time;
    N_Vector  species;
    void *    cvode_mem;
    int       flag;

    /* check number of input/output arguments */
    if (nlhs != 3)
    {  mexErrMsgTxt("syntax: [err_flag, species_out, obsv_out] = network_mex( timepoints, species_init, params )");  }
    if (nrhs != 3)
    {  mexErrMsgTxt("syntax: [err_flag, species_out, obsv_out] = network_mex( timepoints, species_init, params )");  }


    /* make sure timepoints has correct dimensions */
    if ( (mxGetM(prhs[0]) < 2)  ||  (mxGetN(prhs[0]) != 1) )
    {  mexErrMsgTxt("TIMEPOINTS must be a column vector with 2 or more elements.");  }

    /* make sure species_init has correct dimensions */
    if ( (mxGetM(prhs[1]) != 1)  ||  (mxGetN(prhs[1]) != __N_SPECIES__) )
    {  mexErrMsgTxt("SPECIES_INIT must be a row vector with $n_species elements.");  } 

    /* make sure params has correct dimensions */
    if ( (mxGetM(prhs[2]) != 1)  ||  (mxGetN(prhs[2]) != __N_PARAMETERS__) )
    {  mexErrMsgTxt("PARAMS must be a column vector with $n_parameters elements.");  }

    /* get pointers to input arrays */
    timepoints   = mxGetPr(prhs[0]);
    species_init = mxGetPr(prhs[1]);
    parameters   = mxGetPr(prhs[2]);

    /* get number of timepoints */
    n_timepoints = mxGetM(prhs[0]);

    /* Create an mxArray for output trajectories */
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL );
    plhs[1] = mxCreateDoubleMatrix(n_timepoints, __N_SPECIES__, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(n_timepoints, __N_OBSERVABLES__, mxREAL);

    /* get pointers to output arrays */
    return_status   = mxGetPr(plhs[0]);
    species_out     = mxGetPr(plhs[1]);
    observables_out = mxGetPr(plhs[2]);    
   
    /* initialize intermediate data vectors */
    expressions  = NULL;
    expressions = N_VNew_Serial(__N_EXPRESSIONS__);
    if (check_flag((void *)expressions, "N_VNew_Serial", 0))
    {
        return_status[0] = 1;
        return;
    }

    observables = NULL;
    observables = N_VNew_Serial(__N_OBSERVABLES__);
    if (check_flag((void *)observables, "N_VNew_Serial", 0))
    {
        N_VDestroy_Serial(expressions);
        return_status[0] = 1;
        return;
    }

    ratelaws    = NULL; 
    ratelaws = N_VNew_Serial(__N_RATELAWS__);
    if (check_flag((void *)ratelaws, "N_VNew_Serial", 0))
    {   
        N_VDestroy_Serial(expressions);
        N_VDestroy_Serial(observables);        
        return_status[0] = 1;
        return;
    }
    
    /* set up pointers to intermediate data vectors */
    temp_data[0] = expressions;
    temp_data[1] = observables;
    temp_data[2] = ratelaws;

    /* calculate expressions (expressions are constant, so only do this once!) */
    calc_expressions( expressions, parameters );

        
    /* SOLVE model equations! */
    species   = NULL;
    cvode_mem = NULL;

    /* Set the scalar relative tolerance */
    reltol = $cvode_reltol;
    abstol = $cvode_abstol;

    /* Create serial vector for Species */
    species = N_VNew_Serial(__N_SPECIES__);
    if (check_flag((void *)species, "N_VNew_Serial", 0))
    {  
        N_VDestroy_Serial(expressions);
        N_VDestroy_Serial(observables);
        N_VDestroy_Serial(ratelaws);
        return_status[0] = 1;
        return;
    }
    for ( i = 0; i < __N_SPECIES__; i++ )
    {   NV_Ith_S(species,i) = species_init[i];   }
    
    /* write initial species populations into species_out */
    for ( i = 0; i < __N_SPECIES__; i++ )
    {   species_out[i*n_timepoints] = species_init[i];   }
    
    /* write initial observables populations into species_out */ 
    calc_observables( observables, species, expressions );  
    for ( i = 0; i < __N_OBSERVABLES__; i++ )
    {   observables_out[i*n_timepoints] = NV_Ith_S(observables,i);   }

    /*   Call CVodeCreate to create the solver memory:    
     *   CV_ADAMS or CV_BDF is the linear multistep method
     *   CV_FUNCTIONAL or CV_NEWTON is the nonlinear solver iteration
     *   A pointer to the integrator problem memory is returned and stored in cvode_mem.
     */
    cvode_mem = CVodeCreate($cvode_linear_multistep, $cvode_nonlinear_solver);
    if (check_flag((void *)cvode_mem, "CVodeCreate", 0))
    $cvode_cleanup_memory



    /*   Call CVodeInit to initialize the integrator memory:     
     *   cvode_mem is the pointer to the integrator memory returned by CVodeCreate
     *   rhs_func  is the user's right hand side function in y'=f(t,y)
     *   T0        is the initial time
     *   y         is the initial dependent variable vector
     */
    flag = CVodeInit(cvode_mem, calc_species_deriv, timepoints[0], species);
    if (check_flag(&flag, "CVodeInit", 1))
    $cvode_cleanup_memory
   
    /* Set scalar relative and absolute tolerances */
    flag = CVodeSStolerances(cvode_mem, reltol, abstol);
    if (check_flag(&flag, "CVodeSStolerances", 1))
    $cvode_cleanup_memory   
   
    /* pass params to rhs_func */
    flag = CVodeSetUserData(cvode_mem, &temp_data);
    if (check_flag(&flag, "CVodeSetFdata", 1))
    $cvode_cleanup_memory
    
    /* select linear solver */
    $cvode_linear_solver
    $cvode_cleanup_memory
    
    flag = CVodeSetMaxNumSteps(cvode_mem, $cvode_max_num_steps);
    if (check_flag(&flag, "CVodeSetMaxNumSteps", 1))
    $cvode_cleanup_memory

    flag = CVodeSetMaxErrTestFails(cvode_mem, $cvode_max_err_test_fails);
    if (check_flag(&flag, "CVodeSetMaxErrTestFails", 1))
    $cvode_cleanup_memory

    flag = CVodeSetMaxConvFails(cvode_mem, $cvode_max_conv_fails);
    if (check_flag(&flag, "CVodeSetMaxConvFails", 1))
    $cvode_cleanup_memory

    flag = CVodeSetMaxStep(cvode_mem, $cvode_max_step);
    if (check_flag(&flag, "CVodeSetMaxStep", 1))
    $cvode_cleanup_memory

    /* integrate to each timepoint */
    for ( i=1;  i < n_timepoints;  i++ )
    {
        flag = CVode(cvode_mem, timepoints[i], species, &time, CV_NORMAL);
        if (check_flag(&flag, "CVode", 1))
        {
            N_VDestroy_Serial(expressions);
            N_VDestroy_Serial(observables);           
            N_VDestroy_Serial(ratelaws);
            N_VDestroy_Serial(species);
            CVodeFree(&cvode_mem);
            return_status[0] = 1; 
            return;
        }

        /* copy species output from nvector to matlab array */
        for ( j = 0; j < __N_SPECIES__; j++ )
        {   species_out[j*n_timepoints + i] = NV_Ith_S(species,j);   }
        
        /* copy observables output from nvector to matlab array */
        calc_observables( observables, species, expressions );         
        for ( j = 0; j < __N_OBSERVABLES__; j++ )
        {   observables_out[j*n_timepoints + i] = NV_Ith_S(observables,j);   }      
    }
 
    /* Free vectors */
    N_VDestroy_Serial(expressions);
    N_VDestroy_Serial(observables);  
    N_VDestroy_Serial(ratelaws);        
    N_VDestroy_Serial(species);

    /* Free integrator memory */
    CVodeFree(&cvode_mem);

    return;
}


/*  Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns a flag so check if
 *            flag >= 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer 
 */
int check_flag(void *flagvalue, char *funcname, int opt)
{
    int *errflag;

    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    if (opt == 0 && flagvalue == NULL)
    {
        mexPrintf( "\\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\\n", funcname );    
        return(1);
    }

    /* Check if flag < 0 */
    else if (opt == 1)
    {
        errflag = (int *) flagvalue;
        if (*errflag < 0)
        {
            mexPrintf( "\\nSUNDIALS_ERROR: %s() failed with flag = %d\\n", funcname, *errflag );
            return(1);
        }
    }

    /* Check if function returned NULL pointer - no memory allocated */
    else if (opt == 2 && flagvalue == NULL)
    {
        mexPrintf( "\\nMEMORY_ERROR: %s() failed - returned NULL pointer\\n", funcname );
        return(1);
    }

    return(0);
}
EOF
	close(Mexfile);



    # open Mexfile and begin printing...
	open( Mscript, ">$mscript_path" ) or die "Couldn't open $mscript_path: $!\n";
    print Mscript <<"EOF";
function [err, timepoints, species_out, observables_out ] = ${mscript_filebase}( timepoints, species_init, parameters, suppress_plot )
%${mscript_filebase_caps} Integrate reaction network and plot observables.
%   Integrates the reaction network corresponding to the BioNetGen model
%   '${model_name}' and then (optionally) plots the observable trajectories,
%   or species trajectories if no observables are defined. Trajectories are
%   generated using either default or user-defined parameters and initial
%   species values. Integration is performed by the CVode library interfaced
%   to MATLAB via the MEX interface. Before running this script, the model
%   source in file ${mex_filename} must be compiled (see that file for details).
%   ${mscript_filebase_caps} returns an error value, a vector of timepoints,
%   species trajectories, and observable trajectories.
%   
%   [err, timepoints, species_out, observables_out]
%        = ${mscript_filebase}( timepoints, species_init, parameters, suppress_plot )
%
%   INPUTS:
%   -------
%   timepoints      : column vector of time points returned by integrator.
%   species_init    : row vector of $n_species initial species populations.
%   parameters      : row vector of $n_parameters model parameters.
%   suppress_plot   : 0 if a plot is desired (default), 1 if plot is suppressed.
%
%   Note: to specify default value for an input argument, pass the empty array.
%
%   OUTPUTS:
%   --------
%   err             : 0 if the integrator exits without error, non-zero otherwise.
%   timepoints      : a row vector of timepoints returned by the integrator.
%   species_out     : array of species population trajectories
%                        (columns correspond to species, rows correspond to time).
%   observables_out : array of observable trajectories
%                        (columns correspond to observables, rows correspond to time).
%
%   QUESTIONS about the BNG Mfile generator?  Email justinshogg\@gmail.com



%% Process input arguments

% define any missing arguments
if ( nargin < 1 )
    timepoints = [];
end

if ( nargin < 2 )
    species_init = [];
end

if ( nargin < 3 )
    parameters = [];
end

if ( nargin < 4 )
    suppress_plot = 0;
end


% initialize outputs (to avoid error msgs if script terminates early
err = 0;
species_out     = [];
observables_out = [];


% setup default parameters, if necessary
if ( isempty(parameters) )
   parameters = [ $mscript_param_values ];
end
% check that parameters has proper dimensions
if (  size(parameters,1) ~= 1  |  size(parameters,2) ~= $n_parameters  )
    fprintf( 1, 'Error: size of parameter argument is invalid! Correct size = [1 $n_parameters].\\n' );
    err = 1;
    return;
end

% setup default initial values, if necessary
if ( isempty(species_init) )
   species_init = initialize_species( parameters );
end
% check that species_init has proper dimensions
if (  size(species_init,1) ~= 1  |  size(species_init,2) ~= $n_species  )
    fprintf( 1, 'Error: size of species_init argument is invalid! Correct size = [1 $n_species].\\n' );
    err = 1;
    return;
end

% setup default timepoints, if necessary
if ( isempty(timepoints) )
   timepoints = linspace($t_start,$t_end,$n_steps+1)';
end
% check that timepoints has proper dimensions
if (  size(timepoints,1) < 2  |  size(timepoints,2) ~= 1  )
    fprintf( 1, 'Error: size of timepoints argument is invalid! Correct size = [t 1], t>1.\\n' );
    err = 1;
    return;
end

% setup default suppress_plot, if necessary
if ( isempty(suppress_plot) )
   suppress_plot = 0;
end
% check that suppress_plot has proper dimensions
if ( size(suppress_plot,1) ~= 1  |  size(suppress_plot,2) ~= 1 )
    fprintf( 1, 'Error: suppress_plots argument should be a scalar!\\n' );
    err = 1;
    return;
end

% define parameter labels (this is for the user's reference!)
param_labels = { $mscript_param_names };



%% Integrate Network Model
try 
    % run simulation
    [err, species_out, observables_out] = ${mex_filebase}( timepoints, species_init, parameters );
catch
    fprintf( 1, 'Error: some problem integrating ODE network! (CVODE exitflag %d)\\n', err );
    err = 1;
    return;
end



%% Plot Output, if desired

if ( ~suppress_plot )
    
    % define plot labels
$mscript_plot_labels
    % construct figure
$mscript_make_plot
end



%~~~~~~~~~~~~~~~~~~~~~%
% END of main script! %
%~~~~~~~~~~~~~~~~~~~~~%



% initialize species function
function [species_init] = initialize_species( params )

    species_init = zeros(1,$n_species);
$mscript_species_init
end


end
EOF
	close Mscript;
	print "Wrote Mexfile $mex_path and M-file script $mscript_path.\n";
	return ();
}



###
###
###



sub writeMfile_QueryNames
{
	my $model = shift;
	my $plist = $model->ParamList;
	my $slist = $model->SpeciesList;
	my $err;
	
	my $mscript_param_names;
	my $mscript_param_values;
	($mscript_param_names, $mscript_param_values, $err) = $plist->getMatlabConstantNames();
    if ($err) { return ($err) };
    
    my $mscript_observable_names;
	($mscript_observable_names, $err) = $plist->getMatlabObservableNames();
    if ($err) { return ($err) };
    
    my $mscript_species_names;
	($mscript_species_names, $err) = $slist->getMatlabSpeciesNamesOnly();
    if ($err) { return ($err) };
    
    my $q_mscript = 'QueryNames.m';
    
    open(Q_Mscript,">$q_mscript");
    print Q_Mscript <<"EOF";
function [ param_labels, param_defaults, obs_labels, species_labels] = QueryNames( inputlist )
% % Loads all the parameter labels, parameter defaults, observable labels and species labels in the model
% % If generate_network() was executed, then the nanmes of all species are passed
% % If generate_network() was not executed, then the names of the seed speceis are passed

	param_labels = { $mscript_param_names };
	param_defaults = [ $mscript_param_values ];
	obs_labels = { $mscript_observable_names };
	species_labels = { $mscript_species_names };
end

EOF
	close Q_Mscript;
	print "Wrote M-file script $q_mscript.\n";
	return ();
    
}



###
###
###



sub writeMfile_ParametersObservables
{
	# John Sekar created this subroutine
	my $model = shift;
	my $params = (@_) ? shift : {};
	
	my $err;
	
	#Get ref to parameter list
	my $plist = $model->ParamList;
	
	#Names of M-file
	my $par_mscript = 'ParameterList.m';
	my $obs_mscript = 'ObservableList.m';
	
	#Getting param names and observable names
	my $mscript_param_names;
	my $mscript_param_values;
	($mscript_param_names, $mscript_param_values, $err) = $plist->getMatlabConstantNames();
    if ($err) { return ($err) };
    $mscript_param_names =~ s/[\\]//g;
    
    my $mscript_observable_names;
	($mscript_observable_names, $err) = $plist->getMatlabObservableNames();
    if ($err) { return ($err) };
    $mscript_observable_names =~ s/[\\]//g;
  
    
    #Writing parameter list script
	open( Par_Mscript, ">$par_mscript" ) || die "Couldn't open $par_mscript: $!\n";
    print Par_Mscript <<"EOF";
function [outputlist,defaultvals ] = ParameterList( inputlist )
% Used to manipulate and access parameter names
% If inputlist is empty, the entire list of labels is given as output
% If inputlist is a vector of indices, output is a cell array of parameter
% names corresponding to those indices, returns default error if not found
% If inputlist is a cell array of names, output is a vector of indices
% corresponding to those parameter names, returns zero if not found
	param_labels = { $mscript_param_names };
	param_defaults = [ $mscript_param_values ];
	
    param_num = max(size(param_labels));
	
    if nargin < 1
        outputlist = param_labels;
        defaultvals = param_defaults;
        return;
    end

    defaultvals = zeros(size(inputlist));

    if(isnumeric(inputlist))
        outputlist = cell(size(inputlist));
        
        	
        for i=1:1:max(size(inputlist))
            outputlist{i} = param_labels{inputlist(i)}; 
            defaultvals(i) = param_defaults(inputlist(i));
        end
    end
    
   if(iscellstr(inputlist))
       outputlist = zeros(size(inputlist));
       for i=1:1:max(size(inputlist))
           compare = strcmp(inputlist{i},param_labels);
           if(sum(compare)>0)
               outputlist(i) = find(compare,1);
               if(outputlist(i))
                   defaultvals(i) = param_defaults(outputlist(i));
               end
           end
           
       end
	end
end
EOF
	close Par_Mscript;
	print "Wrote M-file script $par_mscript.\n";
	
	#Writing observable list script
	open( Obs_Mscript, ">$obs_mscript" ) || die "Couldn't open $obs_mscript: $!\n";
    print Obs_Mscript <<"EOF";
function [outputlist ] = ObservableList( inputlist )
% Used to manipulate and access observable names
% If inputlist is empty, the entire list of labels is given as output
% If inputlist is a vector of indices, output is a cell array of observable
% names corresponding to those indices, returns default error if not found
% If inputlist is a cell array of names, output is a vector of indices
% corresponding to those observable names, returns zero if not found
	obs_labels = { $mscript_observable_names };
    obs_num = max(size(obs_labels));

    if nargin < 1
        outputlist = obs_labels;
        return;
    end
    
    if(isnumeric(inputlist))
        outputlist = cell(size(inputlist));
        for i=1:1:max(size(inputlist))
            outputlist{i} = obs_labels{inputlist(i)};
        end
    end
    
   if(iscellstr(inputlist))
       outputlist = zeros(size(inputlist));
       for i=1:1:max(size(inputlist))
           compare = strcmp(inputlist{i},obs_labels);
           if(sum(compare)>0)
               outputlist(i) = find(compare,1);
           else
               outputlist(i) = 0;
           end
       end 
end
EOF
	close Obs_Mscript;
	print "Wrote M-file script $obs_mscript.\n";
	return ();
}



###
###
###



sub writeLatex
{
	my $model = shift @_;
	my $params = @_ ? shift @_ : {};

	return '' if $BNGModel::NO_EXEC;

	unless ( $model->RxnList )
    {   return "writeLatex(): No reactions in current model--nothing to do.";   }

    # parameter list
	my $plist = $model->ParamList;

    # model name
	my $model_name = $model->Name;

	# Strip prefixed path
	my $prefix = ( defined $params->{prefix} ) ? $params->{prefix} : $model->getOutputPrefix();
	my $suffix = ( defined $params->{suffix} ) ? $params->{suffix} : undef;
	unless ( $suffix eq '' ) { $prefix .= "_${suffix}"; }

    # latex filename
	my $file = "${prefix}.tex";

    # open file
    my $Lfile;
	open( $Lfile, ">$file" ) or die "Couldn't open $file: $!\n";

	my $version = BNGversion();
	print$Lfile "% Latex formatted differential equations for model $prefix created by BioNetGen $version\n";

	# Document Header
	print $Lfile <<'EOF';
\documentclass{article}
\begin{document}
EOF

	# Dimensions
	my $Nspecies   = scalar @{$model->SpeciesList->Array};
	my $Nreactions = scalar @{$model->RxnList->Array};
	print $Lfile "\\section{Model Summary}\n";
	printf $Lfile "The model has %d species and %d reactions.\n", $Nspecies, $Nreactions;
	print $Lfile "\n";


	# Stoichiometry matrix
	my %S      = ();
	my @fluxes = ();
	my $irxn   = 1;
	foreach my $rxn ( @{$model->RxnList->Array} )
    {
		# Each reactant contributes a -1
		foreach my $r ( @{$rxn->Reactants} )
        {
			--$S{ $r->Index }{$irxn};
		}

		# Each product contributes a +1
		foreach my $p ( @{$rxn->Products} )
        {
			++$S{ $p->Index }{$irxn};
		}
		my ($flux, $err) = $rxn->RateLaw->toLatexString( $rxn->Reactants, $rxn->StatFactor,
			                                               $model->ParamList );
		if ($err) { return $err; }
		push @fluxes, $flux;
		++$irxn;
	}

	print $Lfile "\\section{Differential Equations}\n";
	print $Lfile "\\begin{eqnarray*}\n";
	foreach my $ispec ( sort {$a <=> $b} keys %S )
    {
		printf $Lfile "\\dot{x_{%d}}&=& ", $ispec;
		my $nrxn = 1;
		foreach my $irxn ( sort { $a <=> $b } keys %{ $S{$ispec} } )
        {
            my $mod;
			my $s = $S{$ispec}{$irxn};
			if ( $s == 1 ) {
				$mod = "+";
			}
			elsif ( $s == -1 ) {
				$mod = "-";
			}
			elsif ( $s > 0 ) {
				$mod = "+$s";
			}
			else {
				$mod = "+($s)";
			}
			if ( ($nrxn % 5) == 0 ) { print $Lfile "\\\\ &&"; }
			if ($s)
            {
				printf $Lfile " %s %s", $mod, $fluxes[ $irxn - 1 ];
				++$nrxn;
			}
		}

		if ( $nrxn == 1 ) {	print $Lfile "0"; }
		print $Lfile "\n\\\\\n";
	}
	print $Lfile "\\end{eqnarray*}\n";
	print $Lfile "\n";

	# Document Footer
	print $Lfile <<'EOF';
\end{document}
EOF
	close $Lfile;
	print "Wrote Latex equations to  $file.\n";
	return;
}

sub writeMfile_all
{

# Author: John Sekar
# Writes 3 files
# initSpecies_$name.m which holds the species initialization relations
# $name.m which is the model simulator (same as writeMfile)
# bngModel_$name.m which is a class which inherits bngModel from the MatLab toolbox	
	my $model = shift @_;
	my $params = @_ ? shift @_ : {};

    # a place to hold errors
    my $err;

    # nothing to do if NO_EXEC is true
	return '' if $BNGModel::NO_EXEC;

    # nothing to do if there are no reactions
	unless ( $model->RxnList )
	{
	    return ( "writeMfile() has nothing to do: no reactions in current model.\n"
	            ."  Did you remember to call generate_network() before attempting to\n"
	            ."  write network output?");
	}

    # get reference to parameter list
	my $plist = $model->ParamList;

	
	# get model name
	my $model_name = $model->Name;

    # Build output file name
	# ..use prefix if defined, otherwise use model name
	my $prefix = ( defined $params->{prefix} ) ? $params->{prefix} : $model->getOutputPrefix();
	# ..add suffix, if any
	my $suffix = ( defined $params->{suffix} ) ? $params->{suffix} : undef;
	if ( $suffix )
	{   $prefix .= "_${suffix}";   }

    # split prefix into volume, path and filebase
    my ($vol, $path, $filebase) = File::Spec->splitpath($prefix);
    
	# define m-script file name
    my $mscript_filebase = "${filebase}";
    my $mscript_filename = "${mscript_filebase}.m";
	my $mscript_path     = File::Spec->catpath($vol,$path,$mscript_filename);
    my $mscript_filebase_caps = uc $mscript_filebase;

    # configure options (see Matlab documentation of functions ODESET and ODE15S)
    my $odeset_abstol = 1e-4;
    if ( exists $params->{'atol'} )
    {   $odeset_abstol = $params->{'atol'};  }
    
    my $odeset_reltol = 1e-8;
    if ( exists $params->{'rtol'} )
    {   $odeset_reltol = $params->{'rtol'};  } 

    my $odeset_stats = 'off';
    if ( exists $params->{'stats'} )
    {   $odeset_stats = $params->{'stats'};  } 

    my $odeset_bdf = 'off';
    if ( exists $params->{'bdf'} )
    {   $odeset_bdf = $params->{'bdf'};  }

    my $odeset_maxorder = 5;
    if ( exists $params->{'maxOrder'} )
    {   $odeset_maxorder = $params->{'maxOrder'};  } 

    # time options for mscript
    my $t_start = 0;
    if ( exists $params->{'t_start'} )
    {   $t_start = $params->{'t_start'};  }  

    my $t_end = 10;
    if ( exists $params->{'t_end'} )
    {   $t_end = $params->{'t_end'};  } 

    my $n_steps = 20;
    if ( exists $params->{'n_steps'} )
    {   $n_steps = $params->{'n_steps'};  } 

    # configure time step dependent options
    my $odeset_maxstep = undef;
    if ( exists $params->{'max_step'} )
    {   $odeset_maxstep = $params->{'max_step'};  }     
    
    # construct ODESET function call
    my $mscript_call_odeset;
    if ( defined $odeset_maxstep )
    {
        $mscript_call_odeset = "opts = odeset( 'RelTol',   $odeset_reltol,   ...\n"
                              ."               'AbsTol',   $odeset_abstol,   ...\n"
                              ."               'Stats',    '$odeset_stats',  ...\n"
                              ."               'BDF',      '$odeset_bdf',    ...\n"
                              ."               'MaxOrder', $odeset_maxorder, ...\n"
                              ."               'MaxStep',  $odeset_maxstep    );\n";
    }
    else
    {
        $mscript_call_odeset = "opts = odeset( 'RelTol',   $odeset_reltol,   ...\n"
                              ."               'AbsTol',   $odeset_abstol,   ...\n"
                              ."               'Stats',    '$odeset_stats',  ...\n"
                              ."               'BDF',      '$odeset_bdf',    ...\n"
                              ."               'MaxOrder', $odeset_maxorder   );\n";    
    }

    # Index parameters associated with Constants, ConstantExpressions and Observables
    ($err) = $plist->indexParams();
    if ($err) { return $err };

    # and retrieve a string of expression definitions
    my $n_parameters = $plist->countType( 'Constant' );
    my $n_expressions = $plist->countType( 'ConstantExpression' ) + $n_parameters;
    (my $calc_expressions_string, $err) = $plist->getMatlabExpressionDefs();    
    if ($err) { return $err };

    # get list of parameter names and defintions for matlab
	my $mscript_param_names;
	my $mscript_param_values;
	($mscript_param_names, $mscript_param_values, $err) = $plist->getMatlabConstantNames();
    if ($err) { return $err };

    # get number of species
    my $n_species = scalar @{$model->SpeciesList->Array};
     
	# retrieve a string of observable definitions
    my $n_observables = scalar @{$model->Observables};
    my $calc_observables_string;
    ($calc_observables_string, $err) = $plist->getMatlabObservableDefs();
    if ($err) { return $err };
    
    # get list of observable names for matlab
	my $mscript_observable_names;
	($mscript_observable_names, $err) = $plist->getMatlabObservableNames();
    if ($err) { return $err };
    
    # Construct user-defined functions
    my $user_fcn_declarations = '';
    my $user_fcn_definitions  = '';
	foreach my $param ( @{ $model->ParamList->Array } )
	{
		if ( $param->Type eq 'Function' )
		{
		    # get reference to the actual Function
		    my $fcn = $param->Ref;
		    
		    # don't write function if it depends on a local observable evaluation (this is useless
		    #   since CVode can't do local evaluations)
		    next if ( $fcn->checkLocalDependency($plist) );
		    		    
		    # get function definition			    
		    my $fcn_defn = $fcn->toMatlabString( $plist, {fcn_mode=>'define', indent=>''} );

		    # add definition to the user_fcn_definitions string
		    $user_fcn_definitions .= $fcn_defn . "\n";
        }
	}
	
    # index reactions
    ($err) = $model->RxnList->updateIndex( $plist );
    if ($err) { return $err };

	# retrieve a string of reaction rate definitions
	my $n_reactions = scalar @{$model->RxnList->Array};
    my $calc_ratelaws_string;
    ($calc_ratelaws_string, $err) = $model->RxnList->getMatlabRateDefs( $plist );
    if ($err) { return $err };
    

    # get stoichiometry matrix (sparse encoding in a hashmap)
	my $stoich_hash = {};
	($err) = $model->RxnList->calcStoichMatrix( $stoich_hash );

	# retrieve a string of species deriv definitions
    my $calc_derivs_string;
    ($calc_derivs_string, $err) = $model->SpeciesList->toMatlabString( $model->RxnList, $stoich_hash, $plist );
    if ($err) { return $err };   	


    # get list of species names and initial value expressions for matlab
	my $mscript_species_names;
	my $mscript_species_init;
	($mscript_species_names, $mscript_species_init, $err) = $model->SpeciesList->getMatlabSpeciesNames( $model );
    if ($err) { return $err }; 


    ## Set up MATLAB Plot
    # fontsizes
    my $title_fontsize = 14;
    my $axislabel_fontsize = 12;
    my $legend_fontsize = 10;

    # generate code snippets for plotting observables or species
    my $mscript_plot_labels;
    my $mscript_make_plot;

    # get ylabel (either Number of Concentration)
    my $ylabel;
    if ( $model->SubstanceUnits eq 'Number' )
    {   $ylabel = 'number';   }
    elsif ( $model->SubstanceUnits eq 'Concentration' )
    {   $ylabel = 'concentration';   }
    else
    {   return "writeMfile(): I could not identify model substance units!";   }

    
    if ( @{$model->Observables} )
    {   # plot observables
        $mscript_plot_labels = "    observable_labels = { $mscript_observable_names };\n";
        
        $mscript_make_plot = "    plot(timepoints,observables_out);\n"
                            ."    title('${mscript_filebase} observables','fontSize',${title_fontsize},'Interpreter','none');\n"
                            ."    axis([${t_start} timepoints(end) 0 inf]);\n"
                            ."    legend(observable_labels,'fontSize',${legend_fontsize},'Interpreter','none');\n"
                            ."    xlabel('time','fontSize',${axislabel_fontsize},'Interpreter','none');\n"
                            ."    ylabel('${ylabel}','fontSize',${axislabel_fontsize},'Interpreter','none');\n";
    
    }
    else
    {   # plot species
        $mscript_plot_labels = "    species_labels = { $mscript_species_names };\n";
    
        $mscript_make_plot = "    plot(timepoints,species_out);\n"
                            ."    title('${mscript_filebase} species','fontSize',${title_fontsize},'Interpreter','none');\n"
                            ."    axis([${t_start} timepoints(end) 0 inf]);\n"
                            ."    legend(species_labels,'fontSize',${legend_fontsize},'Interpreter','none');\n"
                            ."    xlabel('time','fontSize',${axislabel_fontsize},'Interpreter','none');\n"
                            ."    ylabel('${ylabel}','fontSize',${axislabel_fontsize},'Interpreter','none');\n";
    }
    


    # open Mfile and begin printing...
	open( Mscript, ">$mscript_path" ) || die "Couldn't open $mscript_path: $!\n";
    print Mscript <<"EOF";
function [err, timepoints, species_out, observables_out ] = ${mscript_filebase}( timepoints, species_init, parameters, suppress_plot )
%${mscript_filebase_caps} Integrate reaction network and plot observables.
%   Integrates the reaction network corresponding to the BioNetGen model
%   '${model_name}' and then (optionally) plots the observable trajectories,
%   or species trajectories if no observables are defined. Trajectories are
%   generated using either default or user-defined parameters and initial
%   species values. Integration is performed by the MATLAB stiff solver
%   'ode15s'. ${mscript_filebase_caps} returns an error value, a vector of timepoints,
%   species trajectories, and observable trajectories.
%   
%   [err, timepoints, species_out, observables_out]
%        = $mscript_filebase( timepoints, species_init, parameters, suppress_plot )
%
%   INPUTS:
%   -------
%   species_init    : row vector of $n_species initial species populations.
%   timepoints      : column vector of time points returned by integrator.
%   parameters      : row vector of $n_parameters model parameters.
%   suppress_plot   : 0 if a plot is desired (default), 1 if plot is suppressed.
%
%   Note: to specify default value for an input argument, pass the empty array.
%
%   OUTPUTS:
%   --------
%   err             : 0 if the integrator exits without error, non-zero otherwise.
%   timepoints      : a row vector of timepoints returned by the integrator.
%   species_out     : array of species population trajectories
%                        (columns correspond to species, rows correspond to time).
%   observables_out : array of observable trajectories
%                        (columns correspond to observables, rows correspond to time).
%
%   QUESTIONS about the BNG Mfile generator?  Email justinshogg\@gmail.com



%% Process input arguments

% define any missing arguments
if ( nargin < 1 )
    timepoints = [];
end

if ( nargin < 2 )
    species_init = [];
end

if ( nargin < 3 )
    parameters = [];
end

if ( nargin < 4 )
    suppress_plot = 0;
end


% initialize outputs (to avoid error msgs if script terminates early
err = 0;
species_out     = [];
observables_out = [];


% setup default parameters, if necessary
if ( isempty(parameters) )
   parameters = [ $mscript_param_values ];
end
% check that parameters has proper dimensions
if (  size(parameters,1) ~= 1  |  size(parameters,2) ~= $n_parameters  )
    fprintf( 1, 'Error: size of parameter argument is invalid! Correct size = [1 $n_parameters].\\n' );
    err = 1;
    return;
end

% setup default initial values, if necessary
if ( isempty(species_init) )
   species_init = initialize_species( parameters );
end
% check that species_init has proper dimensions
if (  size(species_init,1) ~= 1  |  size(species_init,2) ~= $n_species  )
    fprintf( 1, 'Error: size of species_init argument is invalid! Correct size = [1 $n_species].\\n' );
    err = 1;
    return;
end

% setup default timepoints, if necessary
if ( isempty(timepoints) )
   timepoints = linspace($t_start,$t_end,$n_steps+1)';
end
% check that timepoints has proper dimensions
if (  size(timepoints,1) < 2  |  size(timepoints,2) ~= 1  )
    fprintf( 1, 'Error: size of timepoints argument is invalid! Correct size = [t 1], t>1.\\n' );
    err = 1;
    return;
end

% setup default suppress_plot, if necessary
if ( isempty(suppress_plot) )
   suppress_plot = 0;
end
% check that suppress_plot has proper dimensions
if ( size(suppress_plot,1) ~= 1  |  size(suppress_plot,2) ~= 1 )
    fprintf( 1, 'Error: suppress_plots argument should be a scalar!\\n' );
    err = 1;
    return;
end

% define parameter labels (this is for the user's reference!)
param_labels = { $mscript_param_names };



%% Integrate Network Model
 
% calculate expressions
[expressions] = calc_expressions( parameters );

% set ODE integrator options
$mscript_call_odeset

% define derivative function
rhs_fcn = @(t,y)( calc_species_deriv( t, y, expressions ) );

% simulate model system (stiff integrator)
try 
    [timepoints, species_out] = ode15s( rhs_fcn, timepoints, species_init', opts );
catch
    err = 1;
    fprintf( 1, 'Error: some problem encounteredwhile integrating ODE network!\\n' );
    return;
end

% calculate observables
observables_out = zeros( length(timepoints), $n_observables );
for t = 1 : length(timepoints)
    observables_out(t,:) = calc_observables( species_out(t,:), expressions );
end


%% Plot Output, if desired

if ( ~suppress_plot )
    
    % define plot labels
$mscript_plot_labels
    % construct figure
$mscript_make_plot
end


%~~~~~~~~~~~~~~~~~~~~~%
% END of main script! %
%~~~~~~~~~~~~~~~~~~~~~%


% initialize species function
function [species_init] = initialize_species( params )

    species_init = zeros(1,$n_species);
$mscript_species_init
end


% user-defined functions
$user_fcn_definitions


% Calculate expressions
function [ expressions ] = calc_expressions ( parameters )

    expressions = zeros(1,$n_expressions);
$calc_expressions_string   
end



% Calculate observables
function [ observables ] = calc_observables ( species, expressions )

    observables = zeros(1,$n_observables);
$calc_observables_string
end


% Calculate ratelaws
function [ ratelaws ] = calc_ratelaws ( species, expressions, observables )

    ratelaws = zeros(1,$n_observables);
$calc_ratelaws_string
end

% Calculate species derivates
function [ Dspecies ] = calc_species_deriv ( time, species, expressions )
    
    % initialize derivative vector
    Dspecies = zeros($n_species,1);
    
    % update observables
    [ observables ] = calc_observables( species, expressions );
    
    % update ratelaws
    [ ratelaws ] = calc_ratelaws( species, expressions, observables );
                        
    % calculate derivatives
$calc_derivs_string
end


end
EOF

	close(Mscript);
	print "Wrote M-file script $mscript_path.\n";
	
# define m-script file name
     $mscript_filebase = "initSpecies_${filebase}";
     $mscript_filename = "${mscript_filebase}.m";
	 $mscript_path     = File::Spec->catpath($vol,$path,$mscript_filename);
     $mscript_filebase_caps = uc $mscript_filebase;
	 
	 open( Mscript, ">$mscript_path" ) || die "Couldn't open $mscript_path: $!\n";
	print Mscript <<"EOF";
	function [species_init] = initialize_species( params )

    species_init = zeros(1,$n_species);
$mscript_species_init
end
	
EOF

	close(Mscript);
	print "Wrote M-file script $mscript_path.\n";

	$mscript_filebase = "bngModel_${filebase}";
     $mscript_filename = "${mscript_filebase}.m";
	 $mscript_path     = File::Spec->catpath($vol,$path,$mscript_filename);
     $mscript_filebase_caps = uc $mscript_filebase;
	
	 open( Mscript, ">$mscript_path" ) || die "Couldn't open $mscript_path: $!\n";
	print Mscript <<"EOF";
classdef ${mscript_filebase} < bngModel

	properties 
	end
	
	methods
	
	function obj = ${mscript_filebase}
	obj = obj\@bngModel;
	
	species_labels = { $mscript_species_names };
	param_labels =  { $mscript_param_names };
	param_defaults =  [ $mscript_param_values ] ;
	obs_labels = { $mscript_observable_names };
	simulators = { \'$filebase\' };
	
	obj.Param = bngParam(param_labels,param_defaults);
    obj.Obs = bngObs(obs_labels);
    obj.Species = bngSpecies(species_labels,\'initSpecies_${filebase}\');
    obj.simulators = simulators;
	
	end
	
	end


end
	
EOF

	close(Mscript);
	print "Wrote M-file script $mscript_path.\n";

return();
}

###
###
###

1;
