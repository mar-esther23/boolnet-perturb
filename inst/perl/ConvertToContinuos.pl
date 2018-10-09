#!/usr/bin/perl;

#################################################################################################################################
# Developed by: Juan Antonio Arias Del Angel 
# Last version: 1.2 (May 17th 2017)
# Contact: jariasdelangel@gmail.com
#
#
# Description: This program is created with the purpouse to automatize the continuous mapping of boolean network over which further analysis might be performed using R. 
#
# Input: The program takes 3 arguments as input. 
#
# (1) The first argument is a plain text file containing a boolean network encoded according to R-library BoolNet specifications. 
# (2) A string specificating the type of fuzzy logic which will be employed to generate a continuous model. 'Zadeh' specifies Zadeh's logic based in max and min functions and 'Probabilistic' is based in sums and multiplications. 
# (3) A string specificating the type of equation employed to describe the rate of change. "Villarreal" and "SQUAD" specifies for the equations showed by Villarreal et al. 2012 and Ortíz-Gutierrez et al., 2015, respecitvely. 
# (4) A string specificating the separator character between targets and factos as used in the BoolNet format. 
################################################################################################################################

use strict;

######################################################
# In this variables it is stored the information regarding to the logic and equations employed to generate the continuous model.

my $logic = $ARGV[1];
my $equation = $ARGV[2];
my $sep = $ARGV[3];
my $mutants = $ARGV[4];
my $input = $ARGV[5];

######################################################

# Read the BoolNet file and stores the lines.
# Particularly in line 44 the whitespaces are eliminated. 

if($equation eq 'SQUAD' || $equation eq 'Villarreal'){

    open(net, $ARGV[0]);

    my @BooleanNetwork;
    my $BoolFunc;

    while(<net>){

        if(!/targets/ && !/\#/){
            $BoolFunc = $_;
            chomp $BoolFunc;

            $BoolFunc .= '\\';

            $BoolFunc =~ s/ //gi;

            push(@BooleanNetwork, $BoolFunc);
        }


    } close(net);

    ######################################################
    my @BoolFunc;
    my $squad = "dX = ((-exp(0.5*h)+exp(-h*(w_X)))/((1-exp(0.5*h))*(1+exp(-h*(w_X-0.5)))))-(alphaX*X)";
    my $villarreal = "dX = 1/(1+(exp(-2*h*(w_X-b)))) - (alphaX*X)";
    my $eq = ();
    my @nodes = ();
    $mutants = '-' . $mutants . '-';
    $mutants =~ s/,/-/g;
    ######################################################
    # Print the beginning of the R function to simulate the continuous model as employed in R-library deSolve.
print("BoolODE = list()\n");
    print("BoolODE\$func = function(t, state, parameters) {\n");
    print("\twith(as.list(c(state, parameters)),{\n");



    ############################################################################################################
    ############################################################################################################
    print("\t\t\# Input Nodes\n");

    # According to the parameter chosen. The inputs nodes are specified according to fuzzy or probabilistic logic.


    if($logic eq 'Zadeh'){


        for($BoolFunc = 0; $BoolFunc < scalar(@BooleanNetwork); $BoolFunc++){
            @BoolFunc = split($sep, $BooleanNetwork[$BoolFunc]);
            my $characters = length($BoolFunc[1]);

            my $FuzzyExpression = &getFuzzyContent($BoolFunc[1], 0);
            $FuzzyExpression =~ s/!/1-/gi;

            print("\t\tw\_$BoolFunc[0] = $FuzzyExpression\n");

        }
    } 

    if($logic eq 'Probabilistic'){

        for($BoolFunc = 0; $BoolFunc < scalar(@BooleanNetwork); $BoolFunc++){
            @BoolFunc = split($sep, $BooleanNetwork[$BoolFunc]);
            my $characters = length($BoolFunc[1]);

            my $FuzzyExpression = &getProbabilisticContent($BoolFunc[1], 0);
            $FuzzyExpression =~ s/!/1-/gi;

            print("\t\tw\_$BoolFunc[0] = $FuzzyExpression\n");

        }
    }

    ############################################################################################################
    ############################################################################################################
    print "\n\n";

    print("\t\t\# Rates of Change\n");
    if($equation eq 'SQUAD'){
        for($BoolFunc = 0; $BoolFunc < scalar(@BooleanNetwork); $BoolFunc++){
            $eq = $squad;
            @BoolFunc = split($sep, $BooleanNetwork[$BoolFunc]);
            $eq =~ s/X/$BoolFunc[0]/g;
            if((($BoolFunc[0] . '\\') eq $BoolFunc[1]) && $input == 1 ){
                $eq =~ s/= (.+)/= 0/g;
            }
            if(index($mutants, ('-' . $BoolFunc[0] . '-')) != -1){
                $eq =~ s/= (.+)/= 0/g;
            }
            print "\t\t$eq\n";

            push(@nodes, 'd' . $BoolFunc[0]);
        }
    }

    if($equation eq 'Villarreal'){
        for($BoolFunc = 0; $BoolFunc < scalar(@BooleanNetwork); $BoolFunc++){
            $eq = $villarreal;
            @BoolFunc = split($sep, $BooleanNetwork[$BoolFunc]);
            $eq =~ s/X/$BoolFunc[0]/g;
            if((($BoolFunc[0] . '\\') eq $BoolFunc[1]) && $input == 1 ){
                $eq =~ s/= (.+)/= 0/g;
            }
            if($BoolFunc[0] =~ /$mutants/){
                $eq =~ s/= (.+)/= 0/g;
            }
            if(index($mutants, ('-' . $BoolFunc[0] . '-')) != -1){
                $eq =~ s/= (.+)/= 0/g;
            }
            print "\t\t$eq\n";

            push(@nodes, 'd' . $BoolFunc[0]);
        }
    }

    ############################################################################################################
    ############################################################################################################
    # Print the last of the R function to simulate the continuous model as employed in R-library deSolve

    my $nodes = join(',', @nodes);

    print("\n\n");
    print("\t\tlist(c($nodes))\n");
    print("\t})\n}\n");

    # Print a vector containing the default parameters.
    # alphaXi = 1
    # h = 20
    # b = 0.5

    print("BoolODE\$parameters = c(h = 20, b = 0.5");

    for($BoolFunc = 0; $BoolFunc < scalar(@BooleanNetwork); $BoolFunc++){
        @BoolFunc = split($sep, $BooleanNetwork[$BoolFunc]);
        print ", alpha$BoolFunc[0] = 1";
    }


    print("\)\n\n");

    # Print a vector containing the default initial state.
    # Xi = 1

    print("BoolODE\$state = c(");

    for($BoolFunc = 0; $BoolFunc < scalar(@BooleanNetwork); $BoolFunc++){
        @BoolFunc = split($sep, $BooleanNetwork[$BoolFunc]);
        if($BoolFunc == 0){
            print "$BoolFunc[0] = 1"
        } else {
            print ", $BoolFunc[0] = 1";
        }
    }

    print("\)\n");

    print("BoolODE\$times = seq(0, 20,0.1)\n");

    #################################################################################################################################
    #######################                                                                                   #######################
    #######################          Subroutines employed to transform the boolean into fuzzy network         #######################
    #######################                                                                                   #######################
    #################################################################################################################################

    sub getFuzzyContent{

        my $string = $_[0];
        my $pos = $_[1];
        my @variables = ();
        my @characters = ();
        my $position = 0;
        my $fuzzy = ();
        my $nextfuzzy = ();
        my $i = ();
        my $c = ();


        @characters = split('', $BoolFunc[1]);
        for($c = $pos; $c < scalar(@characters); $c++){


    ######################################################
    ######################################################
    # Actions to be taken when a closing bracket is found. 

            if($characters[$c] eq ')'){
                if(scalar(@variables) == 1){
                    $fuzzy = $variables[0];
                }
                elsif(scalar(@variables) == 3){
                    if($variables[1] eq '&'){
                        $fuzzy = "min(" . $variables[0] . "," . $variables[2] . ")";
                    } 

                    if($variables[1] eq '|'){
                        $fuzzy = "max(" . $variables[0] . "," . $variables[2] . ")";
                    }
                }
                elsif(scalar(@variables) >= 5){
                    if($variables[1] eq '&'){
                        $fuzzy = "min(" . $variables[0] . "," . $variables[2] . ")";
                    } 

                    if($variables[1] eq '|'){
                        $fuzzy = "max(" . $variables[0] . "," . $variables[2] . ")";
                    }


                    for($i = 3; $i < scalar(@variables); $i += 2){
                        if($variables[$i] eq '&'){
                            $nextfuzzy = "min(" . $fuzzy . "," . $variables[$i + 1] . ")";
                        } 

                        if($variables[$i] eq '|'){
                            $nextfuzzy = "max(" . $fuzzy . "," . $variables[$i + 1] . ")";
                        }

                        $fuzzy = $nextfuzzy;
                    }



                }
                my $result = $fuzzy . "\t" . $c;
                return($result);
            } 

    ######################################################
    ######################################################
    # Actions to be taken when an openning bracket is found. 
    # This is the recursive part of the function.


            elsif($characters[$c] eq '('){
                my $data = &getFuzzyContent($string, $c + 1);
                my @data = split("\t", $data);
                $variables[$position] .= $data[0];
                $c = $data[1]; 
            } 

    ######################################################
    ######################################################
    # Actions to be taken when binary logic operators are found. 

            elsif($characters[$c] eq '&' || $characters[$c] eq '|'){
                $position++;
                $variables[$position] .= $characters[$c];
                $position++;
            } 

    ######################################################
    ######################################################
    # Actions to be taken when the end of the line is reaached. 

            elsif ($characters[$c] eq '\\'){
               # $variables[$position] .= $characters[$c];

                if(scalar(@variables) == 1){
                    $fuzzy = $variables[0];
                }
                elsif(scalar(@variables) == 3){
                    if($variables[1] eq '&'){
                        $fuzzy = "min(" . $variables[0] . "," . $variables[2] . ")";
                    } 

                    if($variables[1] eq '|'){
                        $fuzzy = "max(" . $variables[0] . "," . $variables[2] . ")";
                    }
                }
                elsif(scalar(@variables) >= 5){
                    if($variables[1] eq '&'){
                        $fuzzy = "min(" . $variables[0] . "," . $variables[2] . ")";
                    } 

                    if($variables[1] eq '|'){
                        $fuzzy = "max(" . $variables[0] . "," . $variables[2] . ")";
                    }

                    for($i = 3; $i < scalar(@variables); $i += 2){
                        if($variables[$i] eq '&'){
                            $nextfuzzy = "min(" . $fuzzy . "," . $variables[$i + 1] . ")";
                        } 

                        if($variables[$i] eq '|'){
                            $nextfuzzy = "max(" . $fuzzy . "," . $variables[$i + 1] . ")";
                        }

                        $fuzzy = $nextfuzzy;
                    }



                }
                return($fuzzy);

            }

    ######################################################
    ######################################################
    # Actions to be taken when other characters are found.

            else {
                $variables[$position] .= $characters[$c];
            }
        }



    }

    ############################################################################################################
    ############################################################################################################
    ############################################################################################################
    ############################################################################################################

    sub getProbabilisticContent{

        my $string = $_[0];
        my $pos = $_[1];
        my @variables = ();
        my @characters = ();
        my $position = 0;
        my $fuzzy = ();
        my $nextfuzzy = ();
        my $i = ();
        my $c = ();
        my $term = ();


        @characters = split('', $BoolFunc[1]);
        for($c = $pos; $c < scalar(@characters); $c++){


    ######################################################
    ######################################################
    # Actions to be taken when a closing bracket is found.
    # An important difference with the other method is that parenthesis are set between variables. 

            if($characters[$c] eq ')'){

                for($term = 0; $term < scalar(@variables); $term++){
                    $variables[$term] =~ s/!/1-/gi;
                }

                if(scalar(@variables) == 1){
                    $fuzzy = "(" . $variables[0] . ")";
                }
                elsif(scalar(@variables) == 3){
                    if($variables[1] eq '&'){
                        $fuzzy = "(" . $variables[0] . ")" . "*" . "(" . $variables[2] . ")";
                    } 

                    if($variables[1] eq '|'){
                        $fuzzy = "(" . $variables[0] . ")" . "+" . "(" . $variables[2] . ")" . "-"  . "(" . $variables[0] . ")" . "*" . "(" . $variables[2] . ")";
                    }
                }
                elsif(scalar(@variables) >= 5){
                    if($variables[1] eq '&'){
                        $fuzzy = "(" . $variables[0] . ")" . "*" . "(" . $variables[2] . ")";
                    } 

                    if($variables[1] eq '|'){
                        $fuzzy = "(" . $variables[0] . ")" . "+" . "(" . $variables[2] . ")" . "-"  . "(" . $variables[0] . ")" . "*" . "(" . $variables[2] . ")";
                    }


                    for($i = 3; $i < scalar(@variables); $i += 2){
                        if($variables[$i] eq '&'){
                            $nextfuzzy = "(" . $fuzzy . ")" . "*" . "(" .  $variables[$i + 1] . ")";
                        } 

                        if($variables[$i] eq '|'){
                            $nextfuzzy = "(" . $fuzzy . ")" . "+" . "(" . $variables[$i + 1] . ")" . "-" . "(" . $fuzzy . ")" . "*" . "(" . $variables[$i + 1] . ")";
                        }

                        $fuzzy = $nextfuzzy;
                    }



                }
                my $result = "(" . $fuzzy . ")" . "\t" . $c;
                return($result);
            } 

    ######################################################
    ######################################################
    # Actions to be taken when an openning bracket is found.
    # This is the recursive part of the function.

            elsif($characters[$c] eq '('){
                my $data = &getProbabilisticContent($string, $c + 1);
                my @data = split("\t", $data);
                $variables[$position] .= $data[0];
                $c = $data[1]; 
            } 

    ######################################################
    ######################################################
    # Actions to be taken when binary logic operators are found.

            elsif($characters[$c] eq '&' || $characters[$c] eq '|'){
                $position++;
                $variables[$position] .= $characters[$c];
                $position++;
            } 

    ######################################################
    ######################################################
    # Actions to be taken when the end of the line is reaached. 
    # An important difference with the other method is that parenthesis are set between variables. 


            elsif ($characters[$c] eq '\\'){
               # $variables[$position] .= $characters[$c];

                for($term = 0; $term < scalar(@variables); $term++){
                    $variables[$term] =~ s/!/1-/gi;
                }

                if(scalar(@variables) == 1){
                    $fuzzy = "(" . $variables[0] . ")";
                }
                elsif(scalar(@variables) == 3){
                    if($variables[1] eq '&'){
                        $fuzzy = "(" . $variables[0] . ")" . "*" . "(" . $variables[2] . ")";
                    } 

                    if($variables[1] eq '|'){
                        $fuzzy = "(" . $variables[0] . ")" . "+" . "(" . $variables[2] . ")" . "-"  . "(" . $variables[0] . ")" . "*" . "(" . $variables[2] . ")";
                    }
                }
                elsif(scalar(@variables) >= 5){
                    if($variables[1] eq '&'){
                        $fuzzy = "(" . $variables[0] . ")" . "*" . "(" . $variables[2] . ")";
                    } 

                    if($variables[1] eq '|'){
                        $fuzzy = "(" . $variables[0] . ")" . "+" . "(" . $variables[2] . ")" . "-"  . "(" . $variables[0] . ")" . "*" . "(" . $variables[2] . ")";
                    }


                    for($i = 3; $i < scalar(@variables); $i += 2){
                        if($variables[$i] eq '&'){
                            $nextfuzzy = "(" . $fuzzy . ")" . "*" . "(" .  $variables[$i + 1] . ")";
                        } 

                        if($variables[$i] eq '|'){
                            $nextfuzzy = "(" . $fuzzy . ")" . "+" . "(" . $variables[$i + 1] . ")" . "-" . "(" . $fuzzy . ")" . "*" . "(" . $variables[$i + 1] . ")";
                        }

                        $fuzzy = $nextfuzzy;
                    }



                }
                return($fuzzy);

            }
    ######################################################
    ######################################################
    ######################################################
    ######################################################
    # Actions to be taken when other characters are found.

            else {
                $variables[$position] .= $characters[$c];
            }
        }



    }
}

#################################################################################################################################
#################################################################################################################################

if($equation eq 'HillCube'){
    my $node = ();
    my @regs = ();
    my @true = ();
    my @nodes = ();
    my $nodes = ();

    my $Hill_global = "";

    my @Hill_column = ();

    my @Hill_column_variable = ();
    my $Hill_column_variable = ();
    
    ############################################################################################################
    ############################################################################################################
    
    print("BoolODE = list()\n");
    print("BoolODE[[1]] = function(t, state, parameters) {\n");
    print("\twith(as.list(c(state, parameters)),{\n");
    print("\t\t\# Input Nodes\n");
    
    ############################################################################################################
    ############################################################################################################
    open(HillCube, $ARGV[0]);

    while(<HillCube>){
        if(/(.+)\t(.+)\t(.+)/){
            $node = $1;
            @regs = split(",", $2);
            @true = split(",", $3);
            
            push(@nodes, 'd' . $node);
            print "\t$node\.t = 0\n";
            
            @Hill_column = ();
            for(my $true_row = 0; $true_row < scalar(@true); $true_row++){
                my @condition = &dec2bin($true[$true_row], scalar(@regs));
                @Hill_column_variable = ();

                for(my $column = 0; $column < scalar(@condition); $column++){
                    if($condition[$column] == 0){
                        push(@Hill_column_variable, "(1 - $regs[$column]**n/(K**n + $regs[$column]**n))");
                    }
                    else{
                        push(@Hill_column_variable, "($regs[$column]**n/(K**n + $regs[$column]**n))");
                    }
                }

                $Hill_column_variable = join('*', @Hill_column_variable);
                push(@Hill_column, $Hill_column_variable);
                print "\t$node\.t = $node\.t + $Hill_column_variable\n\n";
                
            }
            print "\td$node = ($node\.t - alpha*$node)\n\n";
            #$Hill_global = join(' + ', @Hill_column);
            #print "d$node = (1/tau_$node)*($Hill_global - alpha\_$node*$node)\n";

        }
        
       
    } close(HillCube);
    
    $nodes = join(',', @nodes);
        
    print("\n\n");
    print("\t\tlist(c($nodes))\n");
    print("\t})\n}\n");
        

}
#################################################################################################################################
#################################################################################################################################


sub dec2bin{
    my $dec = $_[0];
    my $len = $_[1];
    my $bit = 0;
    my @bin = ();
    
    for(my $pos = 0; $pos < $len; $pos++){
        $bit = $dec & 1;
        $dec = $dec >> 1;
        push(@bin, $bit);
        
    }
    
    return(reverse @bin);
}