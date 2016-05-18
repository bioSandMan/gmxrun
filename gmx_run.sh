#!/bin/bash

# Script authored by FM Ytreberg <fmytreberg@gmail.com>
# Last modified Jan 2015
#
# For usage information run with '-h' flag.

# Requires:
#   (i)   input.pdb
#   (ii)  em.mdp, pr.mdp, and md.mdp (can be copied from ~/res/md/gmx/mdp/)
#   (iii) gromacs
#   (iv)  (optional) martini ff (http://md.chem.rug.nl/cgmartini/)
#   (v)   (optional) DSSP needed for martini
#   (vi)  (optional) sirah ff (http://www.sirahff.com/)
#   (vii) (optional) .gro and .itp files for methane or ethane solvation

### May need to change variables below ###
mtheth="/home/marty/res/md/gmx/solvents"
gmx="/usr/local/gromacs-5.0/bin/gmx-5.0"
mdrun="/usr/local/gromacs-5.0/bin/gmx-5.0 mdrun"
martini="/home/marty/apps/martini"
dssp="/home/marty/apps/dssp-2.0.4/dssp"
sirah="/home/marty/apps/sirah.ff"
if [ ! -z `hostname | grep compute` ] || [ ! -z `hostname | grep coan` ]; then
    mtheth="/home/marty/res/md/gmx/solvents"
    gmx="/share/apps/gromacs-5.0/bin/gmx-5.0"
    mdrun="/opt/openmpi/bin/mpirun -machinefile nodefile.txt\
            /share/apps/gromacs-5.0/bin/mdrun-5.0_mpi"
    gpumdrun="/share/apps/gromacs-5.0/bin/mdrun-5.0_gpu"
    martini="/share/apps/martini"
    dssp="/share/apps/dssp-2.0.4/dssp"
    sirah="/share/apps/sirah.ff"
elif [ ! -z `hostname | grep fortytwo` ]; then
    mtheth="/mnt/home/mytrebe/res/md/gmx/solvents"
    gmx="gmx_mpi" # PBS script will need to load module
    mdrun="mpirun -machinefile nodefile.txt mdrun_mpi"
fi
#########################################

#
# PARSE COMMAND LINE
#

show_help() {
cat << EOF
usage: ${0##*/} [options]

With -b this script will take an input pdb file and then generate .top and .gro
files. Full run creates a waterbox, fills it with water and 0.15 salt
concentration, performs energy minimization, position restrained and
unrestrained simulations, respectively. Full run requires em.mdp, pr.mdp and
md.mdp files be in the same directory as input file. With -r this script will
take input files .gro, .top and .mdp with same base name and run a simulation.
With -a this script will take input files .gro, .top and .mdp with same base
name and run an alchemical free energy simulation. The .mdp file is required
to contain init-lambda-state = ILS.
    
optional arguments:
  -h          display this help and exit
  -c          clean files and exit. WARNING: will remove all files in directory
              that differ from start of script run
  -b          (default=false, requires -i and -f) build and run equilibration
              simulations from a single PDB or PQR file
  -r          (default=false, requires -i) run a simulation
  -a          (default=false, requires -i) run alchemical free energy simulation
  -i INPUT    with -b this is input PDB or PQR file. with -r or -a this is
              basename for input files
  -f FFIELD   with -b this is forcefield choice:
              'amber'=AMBER99SB*-ILDNP+TIP3P,  'gromos'=GROMOS54a7+SPC,
              'charmm'=CHARMM22*-CMAP+TIP3P,   'gbsa'=AMBER99SB*-ILDNP+GBSA,
              'martini'=MARTINI+W,             'pmartini'=MARTINI+PW,
              'sirah'=SIRAH+W,
              'meth'=AMBER99SB*-ILDNP+methane, 'eth'=AMBER99SB*-ILDNP+ethane
  -t GRO/PDB  with -b this will generate .top and .gro/.pdb for protein with
              no solvent and no sims
  -o OUTBASE  (default=out) base name for output files
  -w WATERPAD (default=1.0) distance between edge of box and solute in nm.
              choose >=0.8
  -n NCORE    (default=auto) number of cores to use. do not use this flag for
              MPI or GPU runs on a cluster.
  -p POFFSET  (default=0) specify gmx mdrun -pinoffset flag. this is
              important if running more than one job on a single machine.
              should be number of running threads gmx mdrun -nt. (if one -nt 10
              job is already running use -p 10, if two -nt 8 jobs are already
              running use -p 16)
  -g          (default=false) use GPU version of mdrun
EOF
}                

output="out"
watpad="1.0"
np=0
pinoffset=0
flagsum=0
buildflag=false
runflag=false
alchemflag=false
topflag=false
cleanflag=false
gpuflag=false
while getopts ":hcbrai:f:t:o:w:n:p:g" opt; do
    case $opt in
        h)
            show_help
            exit
            ;;
        c)
            cleanflag=true
            ((flagsum=flagsum+1))
            ;;
        b)
            buildflag=true
            ((flagsum=flagsum+1))
            ;;
        r)
            runflag=true
            ((flagsum=flagsum+1))
            ;;
        a)
            alchemflag=true
            ((flagsum=flagsum+1))
            ;;
        i)
            in=$OPTARG
            ;;
        f)
            ff=$OPTARG
            ;;
        t)
            topflag=$OPTARG
            ;;
        o)
            output=$OPTARG
            ;;
        w)
            watpad=$OPTARG
            ;;
        n)
            np=$OPTARG
            ;;
        p)
            pinoffset=$OPTARG
            ;;
        g)
            gpuflag=true
            ;;
        \?)
            echo "error: invalid option: -$OPTARG"
            exit
            ;;
    esac
done

cleanfile="gmx_run.clean"
if [ "$OPTIND" -eq 1 ]; then  # Show help if no arguments are specified.
    show_help
    exit
elif [ "$cleanflag" == true ]; then
    rm -f $(<"$cleanfile")
    rm -f $cleanfile nohup.out
    echo "Script files cleaned"
    exit
fi

# Check flag dependencies.
if [ "$flagsum" -eq 0 ]; then
    echo "error: must specify one of these flags: -c, -r, -b, -a"
    exit
elif [ "$flagsum" -gt 1 ]; then
    echo "error: can only specify one of these flags: -c, -r, -b, -a"
    exit
fi
if [ "$buildflag" == true ]; then
    if [ -z "$ff" ] || [ -z "$in" ]; then
        echo "error: must specify both -i and -f for -b option"
        exit
    fi
fi

# Look for necesary files.
if [ "$buildflag" == true ]; then
    if [ ! -e "$in" ]; then
        echo "error: $in file not found"
        exit
    fi
    if [ "$topflag" == false ]; then
        if [ "$ff" == "gbsa" ]; then
            if [ ! -e "em.mdp" ] || [ ! -e "md.mdp" ]; then
	            echo "error: -b option requires em.mdp and md.mdp files"
	            exit
            fi
        else
            if [ ! -e "em.mdp" ] || [ ! -e "pr.mdp" ] || [ ! -e "md.mdp" ]; then
	            echo "error: -b option requires em.mdp pr.mdp and md.mdp files"
	            exit
            fi
        fi
    fi
else
    if [ ! -e "$in.mdp" ] || [ ! -e "$in.gro" ] || [ ! -e "$in.top" ]; then
        echo "error: -r option requires $in.mdp, $in.gro and $in.top"\
                "files"
        exit
    fi
fi

#
# SET FILE NAMES
#

# Generate list of files that are in directory before run.
ls_before=`ls -t`

# Create a log file, or quit if one already exists.
log="gmx_run.log"
if [ -e "$log" ]; then
    echo "error: $log file found: will not erase"
    exit
fi
echo "** Running $0 $@ **" > $log
echo >> $log

# Set names for output files and mdrun.
if [ "$gpuflag" == true ]; then
    mdrun="$gpumdrun -ntomp $np -pin on" 
else
    mdrun="$mdrun -nt $np -ntomp 1 -pin on -pinoffset $pinoffset -pinstride 1"
fi
top="$output.top"
gro="$output.gro"
ndx="$output.ndx"
pdb="$output.pdb"


#
# MAIN PROGRAM
#

if [ "$buildflag" == true ]; then
    echo "Running with -b option"

    # Generate dodecahdral water box
    if [ "$ff" == "amber" ]; then
        echo "  using AMBER99SB*-ILDNP and TIP3P"
        if [ "$topflag" == "PDB" ]; then
            $gmx pdb2gmx -vsite hydrogens -ff amber99sb-star-ildnp -water none\
                    -ignh -f $in -o $pdb -p $top -i $output.posre >> $log 2>&1 
            echo "Topology and structure files generated for protein"
            exit
        elif [ "$topflag" == "GRO" ]; then
            $gmx pdb2gmx -vsite hydrogens -ff amber99sb-star-ildnp -water none\
                    -ignh -f $in -o $gro -p $top -i $output.posre >> $log 2>&1 
            echo "Topology and structure files generated for protein"
            exit
        elif [ "$topflag" != false ]; then
            echo "error: invalid entry for -t flag"
            exit
        fi
        $gmx pdb2gmx -vsite hydrogens -ff amber99sb-star-ildnp -water tip3p\
                -ignh -f $in -o $gro -p $top -i $output.posre >> $log 2>&1 
        $gmx editconf -f $gro -o $gro -bt dodecahedron -d $watpad >> $log 2>&1
        $gmx solvate -cs -cp $gro -o $gro -p $top >> $log 2>&1
    elif [ "$ff" == "charmm" ]; then
        echo "  using CHARMM22*-CMAP and TIP3P"
        if [ "$topflag" == "PDB" ]; then
            $gmx pdb2gmx -vsite hydrogens -ff charmm22star -water none -ignh\
                    -f $in -o $pdb -p $top -i $output.posre >> $log 2>&1 
            echo "Topology and structure files generated for protein"
            exit
        elif [ "$topflag" == "GRO" ]; then
            $gmx pdb2gmx -vsite hydrogens -ff charmm22star -water none -ignh\
                    -f $in -o $gro -p $top -i $output.posre >> $log 2>&1 
            echo "Topology and structure files generated for protein"
            exit
        elif [ "$topflag" != false ]; then
            echo "error: invalid entry for -t flag"
            exit
        fi
        $gmx pdb2gmx -vsite hydrogens -ff charmm22star -water tip3p -ignh\
                -f $in -o $gro -p $top -i $output.posre >> $log 2>&1 
        $gmx editconf -f $gro -o $gro -bt dodecahedron -d $watpad >> $log 2>&1
        $gmx solvate -cs -cp $gro -o $gro -p $top >> $log 2>&1
    elif [ "$ff" == "sirah" ]; then
        echo "  using SIRAH and W"
        if [ "$topflag" != false ]; then
            echo "error: topology generation not implemented for forcefield"
            exit
        fi
        ln -s $sirah .
        ln -s $sirah/residuetypes.dat .
        ln -s $sirah/specbond.dat .
        $sirah/tools/CGCONV/cgconv.pl -i $in -o ${gro}.pdb >> $log 2>&1
        echo 1 | $gmx pdb2gmx -ff sirah -f ${gro}.pdb -o $gro -p $top\
                -i $output.posre >> $log 2>&1 
        $gmx editconf -f $gro -o $gro -bt dodecahedron -d $watpad >> $log 2>&1
        $gmx solvate -cs $sirah/wt416.gro -cp $gro -o $gro -p $top >> $log 2>&1
        nwat=`grep "Generated solvent" $log | awk '{print $7}'`
        sed -i '$a\WT4                '$nwat'' $top 
    elif [ "$ff" == "pmartini" ]; then
        echo "  using MARTINI and PW"
        if [ "$topflag" != false ]; then
            echo "error: topology generation not implemented for forcefield"
            exit
        fi
        $gmx pdb2gmx -ff oplsaa -water none -ignh -f $in -o tmp.pdb\
                -i $output.posre >> $log 2>&1 
        rm -f topol* posre*
        $martini/martinize.py -f tmp.pdb -o $top -x ${gro}.pdb\
                -dssp $dssp -p backbone -ff martini21p >> $log 2>&1
        rm -f tmp.pdb
        $gmx editconf -f ${gro}.pdb -o $gro -bt dodecahedron -d $watpad\
                >> $log 2>&1
        $gmx solvate -cs $martini/waterP.gro -cp $gro -o $gro -p $top\
                >> $log 2>&1
        nwat=`grep "Generated solvent" $log | awk '{print $7}'`
        sed -i '1c\#include "'$martini'/martini_v2.P.itp"' $top 
        sed -i '2i\#include "'$martini'/martini_v2.1P_aminoacids.itp"' $top 
        sed -i '$a\PW             '$nwat'' $top 
    elif [ "$ff" == "martini" ]; then
        echo "  using MARTINI and W"
        if [ "$topflag" != false ]; then
            echo "error: topology generation not implemented for forcefield"
            exit
        fi
        $gmx pdb2gmx -ff oplsaa -water none -ignh -f $in -o tmp.pdb\
                -i $output.posre >> $log 2>&1 
        rm -f topol* posre*
        $martini/martinize.py -f tmp.pdb -o $top -x ${gro}.pdb\
                -dssp $dssp -p backbone -ff martini22 >> $log 2>&1
        rm -f tmp.pdb
        $gmx editconf -f ${gro}.pdb -o $gro -bt dodecahedron -d $watpad\
                >> $log 2>&1
        $gmx solvate -cs $martini/water.gro -cp $gro -o $gro -p $top\
                >> $log 2>&1
        nwat=`grep "Generated solvent" $log | awk '{print $4}'`
        sed -i '1c\#include "'$martini'/martini_v2.2.itp"' $top 
        sed -i '2i\#include "'$martini'/martini_v2.2_aminoacids.itp"' $top 
        sed -i '$a\W             '$nwat'' $top 
    elif [ "$ff" == "gbsa" ]; then
        echo "  using AMBER99SB*-ILDNP and GBSA"
        if [ "$topflag" == "PDB" ]; then
            $gmx pdb2gmx -vsite hydrogens -ff amber99sb-star-ildnp -water none\
                    -ignh -f $in -o $pdb -p $top -i $output.posre >> $log 2>&1 
            echo "Topology and structure files generated for protein"
            exit
        elif [ "$topflag" == "GRO" ]; then
            $gmx pdb2gmx -vsite hydrogens -ff amber99sb-star-ildnp -water none\
                    -ignh -f $in -o $gro -p $top -i $output.posre >> $log 2>&1 
            echo "Topology and structure files generated for protein"
            exit
        elif [ "$topflag" != false ]; then
            echo "error: invalid entry for -t flag"
            exit
        fi
        $gmx pdb2gmx -vsite hydrogens -ff amber99sb-star-ildnp -water none\
                -ignh -f $in -o $gro -p $top -i $output.posre >> $log 2>&1 
    elif [ "$ff" == "gromos" ]; then
        echo "  using GROMOS 54a7 and SPC"
        if [ "$topflag" == "PDB" ]; then
            $gmx pdb2gmx -vsite hydrogens -ff gromos54a7 -water none\
                    -ignh -f $in -o $pdb -p $top -i $output.posre >> $log 2>&1 
            echo "Topology and structure files generated for protein"
            exit
        elif [ "$topflag" == "GRO" ]; then
            $gmx pdb2gmx -vsite hydrogens -ff gromos54a7 -water none\
                    -ignh -f $in -o $gro -p $top -i $output.posre >> $log 2>&1 
            echo "Topology and structure files generated for protein"
            exit
        elif [ "$topflag" != false ]; then
            echo "error: invalid entry for -t flag"
            exit
        fi
        $gmx pdb2gmx -vsite hydrogens -ff gromos54a7 -water spc -ignh\
                -f $in -o $gro -p $top -i $output.posre >> $log 2>&1 
        $gmx editconf -f $gro -o $gro -bt dodecahedron -d $watpad >> $log 2>&1
        $gmx solvate -cs -cp $gro -o $gro -p $top >> $log 2>&1
    elif [ "$ff" == "meth" ]; then
        echo "  using AMBER99SB*-ILDNP and methane"
        if [ "$topflag" == "PDB" ]; then
            $gmx pdb2gmx -ff amber99sb-star-ildnp -water none -ignh -f $in\
                    -o $pdb -p $top -i $output.posre >> $log 2>&1 
            echo "Topology and structure files generated for protein"
            exit
        elif [ "$topflag" == "GRO" ]; then
            $gmx pdb2gmx -ff amber99sb-star-ildnp -water none -ignh -f $in\
                    -o $gro -p $top -i $output.posre >> $log 2>&1 
            echo "Topology and structure files generated for protein"
            exit
        elif [ "$topflag" != false ]; then
            echo "error: invalid entry for -t flag"
            exit
        fi
        # Can't use -vsite option since methane .itp doesn't define them
        $gmx pdb2gmx -ff amber99sb-star-ildnp -water none -ignh\
                -f $in -o $gro -p $top -i $output.posre >> $log 2>&1 
        $gmx editconf -f $gro -o $gro -bt dodecahedron -d $watpad >> $log 2>&1
        $gmx solvate -cs $mtheth/mth.gro -cp $gro -o $gro -p $top >> $log 2>&1
        nsolv=`grep "Generated solvent" $log | awk '{print $7}'`
        sed -i '/amber99sb-star-ildnp.ff/a #include "'$mtheth'/mth.itp"' $top
        sed -i '$a\methane           '$nsolv'' $top 
    elif [ "$ff" == "eth" ]; then
        echo "  using AMBER99SB*-ILDNP and ethane"
        if [ "$topflag" == "PDB" ]; then
            $gmx pdb2gmx -ff amber99sb-star-ildnp -water none -ignh -f $in\
                    -o $pdb -p $top -i $output.posre >> $log 2>&1 
            echo "Topology and structure files generated for protein"
            exit
        elif [ "$topflag" == "GRO" ]; then
            $gmx pdb2gmx -ff amber99sb-star-ildnp -water none -ignh -f $in\
                    -o $gro -p $top -i $output.posre >> $log 2>&1 
            echo "Topology and structure files generated for protein"
            exit
        elif [ "$topflag" != false ]; then
            echo "error: invalid entry for -t flag"
            exit
        fi
        # Can't use -vsite option since ethane .itp doesn't define them
        $gmx pdb2gmx -ff amber99sb-star-ildnp -water none -ignh\
                -f $in -o $gro -p $top -i $output.posre >> $log 2>&1 
        $gmx editconf -f $gro -o $gro -bt dodecahedron -d $watpad >> $log 2>&1
        $gmx solvate -cs $mtheth/eth.gro -cp $gro -o $gro -p $top >> $log 2>&1
        nsolv=`grep "Generated solvent" $log | awk '{print $7}'`
        sed -i '/amber99sb-star-ildnp.ff/a #include "'$mtheth'/eth.itp"' $top
        sed -i '$a\ethane           '$nsolv'' $top 
    else
        echo "error: invalid forcefield choice: $ff"
        exit 1
    fi

    # Neutralize water box
    if [ "$ff" != "gbsa" ] && [ "$ff" != "meth" ] && [ "$ff" != "eth" ]; then
        echo "  adding ions to make 0.15 mol/liter concentration"
        $gmx grompp -f em -po em.mdout -c $gro -p $top -o em >> $log 2>&1
        if [ "$ff" == "martini" ]; then
            echo W | $gmx genion -s em -o $gro -conc 0.15 -neutral\
                    -pname NA -nname CL >> $log 2>&1
            sed -i '3i\#include "'$martini'/martini_v2.0_ions.itp"' $top
            nna=`tail $log | grep 'Will try to add' | awk '{print $5}'`
            ncl=`tail $log | grep 'Will try to add' | awk '{print $9}'`
            ((nwat=nwat-nna-ncl))
            sed -i '$c\W              '$nwat'' $top 
            sed -i '$a\NA             '$nna'' $top 
            sed -i '$a\CL             '$ncl'' $top 
        elif [ "$ff" == "pmartini" ]; then
            echo PW | $gmx genion -s em -o $gro -conc 0.15 -neutral\
                    -pname NA -nname CL >> $log 2>&1
            sed -i '3i\#include "'$martini'/martini_v2.0_ions.itp"' $top
            nna=`tail $log | grep 'Will try to add' | awk '{print $5}'`
            ncl=`tail $log | grep 'Will try to add' | awk '{print $9}'`
            ((nwat=nwat-nna-ncl))
            sed -i '$c\PW              '$nwat'' $top 
            sed -i '$a\NA             '$nna'' $top 
            sed -i '$a\CL             '$ncl'' $top 
        elif [ "$ff" == "sirah" ]; then
            echo WT4 | $gmx genion -s em -o $gro -conc 0.15 -neutral\
                    -pname NaW -nname ClW >> $log 2>&1
            nna=`tail $log | grep 'Will try to add' | awk '{print $5}'`
            ncl=`tail $log | grep 'Will try to add' | awk '{print $9}'`
            ((nwat=nwat-nna-ncl))
            sed -i '$c\WT4             '$nwat'' $top 
            sed -i '$a\NaW             '$nna'' $top 
            sed -i '$a\ClW             '$ncl'' $top 
        else
            echo SOL | $gmx genion -s em -o $gro -p $top -conc 0.15 -neutral\
                    -pname NA -nname CL >> $log 2>&1
        fi
    fi

    # Make index file
    echo q | $gmx make_ndx -f $gro -o $ndx >> $log 2>&1

    # Perform energy minimization
    echo "  performing energy minimization with em.mdp"
    $gmx grompp -f em -po em.mdout -c $gro -p $top -o em >> $log 2>&1
    $mdrun -v -cpo em -s em -o em -x em -c em.last -e em -g em >> $log 2>&1
    
    # Perform position restrained and unrestrained simulations
    if [ "$ff" != "gbsa" ]; then
        echo "  performing position restrained simulation with pr.mdp"
        $gmx grompp -f pr -po pr.mdout -c em.last -p $top -o pr >> $log 2>&1
        $mdrun -v -cpo pr -s pr -o pr -x pr -c pr.last -e pr -g pr >> $log 2>&1
        echo "  performing unrestrained simulation with md.mdp"
        $gmx grompp -f md -po md.mdout -c pr.last -p $top -o md >> $log 2>&1
        $mdrun -v -cpo md -s md -o md -x md -c md.last -e md -g md >> $log 2>&1
    else
        echo "  performing unrestrained simulation with md.mdp"
        $gmx grompp -f md -po md.mdout -c em.last -p $top -o md >> $log 2>&1
        $mdrun -v -cpo md -s md -o md -x md -c md.last -e md -g md >> $log 2>&1
    fi

    cp md.last.gro $gro  # Update the output .gro file.

elif [ "$runflag" == true ]; then
    echo "Running with -r option"
    echo "  performing simulation with $in.mdp"
    if [ -e "$in.ndx" ]; then
        $gmx grompp -n $in -f $in -po $output.mdout -c $in -p $in -o $output\
                >> $log 2>&1
    else
        $gmx grompp -f $in -po $output.mdout -c $in -p $in -o $output\
                >> $log 2>&1
    fi
    $mdrun -cpo $output -s $output -o $output -x $output -c $output.last\
        -e $output -g $output >> $log 2>&1

elif [ "$alchemflag" == true ]; then
    echo "Running with -a option"
    nlambda=`grep coul-lambdas $in.mdp | wc -w`
    ((nlambda=nlambda-2))
    lastgro="$in.gro"
    echo "  performing alchemical sims with $nlambda lambda values"
    for ((i=0; i<$nlambda; i++)); do
        echo "    running lambda column $i"
        ibase="${output}.${i}"
        sed 's/ILS/'$i'/' $in.mdp > ${ibase}.mdp
        $gmx grompp -maxwarn 3 -f $ibase -po $ibase.mdout -c $lastgro\
                -p $in -o $ibase >> $log 2>&1
        $mdrun -s $ibase -o $ibase -x $ibase -c $ibase.last -e $ibase\
                -g $ibase -cpo $ibase -dhdl $ibase.dhdl >> $log 2>&1
        lastgro="$ibase.last.gro"
    done

    # Process output files using pymbar
    dt=`grep 'dt ' $in.mdp | awk '{print $3}'`
    nsteps=`grep 'nsteps ' $in.mdp | awk '{print $3}'`
    begin=`echo "$dt*$nsteps/2" | bc`
    temp=`grep ref_t $in.mdp | awk '{print $3}'`
    python /home/marty/apps/pymbar-examples/alchemical-free-energy/alchemical-gromacs.py -s $begin -t $temp -p $output -v >> $log 2>&1
    echo "MBAR estimate in kJ/mol is (half of trajectory discarded):"
    grep TOTAL pymbar_results.txt | awk '{print $17, $18, $19}'

fi


# Produce list of files generated since beginning of script (for -c option)
ls_after=`ls -t`
comm -23 <(echo "$ls_after" | sort) <(echo "$ls_before" | sort)\
        >> $cleanfile

echo "Script complete. Can remove generated files with -c option."
