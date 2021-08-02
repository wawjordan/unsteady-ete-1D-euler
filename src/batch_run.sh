#!./bin/bash
input="input.nml"
temp="temp.txt"
sdir="../results/summaries/"
start=`date +%s`
shock_str=""
cfl_str=""
limiter_str=""
p0=300.0
T0=600.0
prat=0.4
shock=1
cons=T

flux=3 # 1 2 3
limiter=1
eps_MUSCL=1.0 #0.0 1.0
#kappa_MUSCL=-1.0

eps_roe=0.01
beta_lim=1.0 #1.0 1.5 2.0
cfl=0.1
k2=0.5 #0.5 0.4 0.3 0.25
k4=0.03125 #0.03125 0.02 0.015625
maxk=200000
Sout=200000
Rout=100
disp_out=1000
mkdir -p "$sdir"
for kappa_MUSCL in -1.0 0.0 0.5 1.0 #-1.0 0.0 0.5 1.0
#for k2 in 0.5 0.4 0.3 0.25 
do
 for limiter in 1 2 3 4
 #for k4 in 0.03125 0.02 0.015625
 do
    if [ $flux -eq 1 ]; then
      summary="${sdir}summary_K${k2}_K_${k4}.dat"
    else
      summary="${sdir}summary_K${kappa_MUSCL}_lim_${limiter}.dat"
    fi
    touch $summary
    if [ -f $summary ]; then
      rm -f "$summary"
    fi
  for imax in 128 #16 32 64 128 256 512
  do

if [ $shock -eq 0 ]; then
  shock_str="isentropic"
else
  shock_str="normal shock"
fi
if [ $flux -eq 1 ]; then
  flux_str="central scheme"
elif [ $flux -eq 2 ]; then
  flux_str="van Leer flux"
elif [ $flux -eq 3 ]; then
  flux_str="Roe's flux"
fi

if [ $flux -ne 1 ]; then
if [ $limiter -eq 1 ]; then
  limiter_str="van Leer limiter"
elif [ $limiter -eq 2 ]; then
  limiter_str="van Albada limiter"
elif [ $limiter -eq 3 ]; then
  limiter_str="minmod limiter"
elif [ $limiter -eq 4 ]; then
  limiter_str="beta ($beta_lim) limiter"
fi
fi
echo ""
echo "!==============================================================================!"
echo "| N=$imax | $shock_str | P_rat=$prat "
echo "--------------------------------------------------------------------------------"
echo "| $flux_str | $limiter_str | CFL=$cfl "
echo "!==============================================================================!"

echo "&grid" > $input
echo "  imax = $imax" >> $input
echo "  xmin = -1.0" >> $input
echo "  xmax =  1.0" >> $input
echo "  n_ghost_cells = 2" >> $input
echo "/" >> $input
echo "" >> $input
echo "&geometry" >> $input
echo "  Astar = 0.2" >> $input
echo "/" >> $input
echo "" >> $input
echo "&constants" >> $input
echo "  gamma = 1.4" >> $input
echo "/" >> $input
echo "" >> $input
echo "&initial" >> $input
echo "  p0 = $p0" >> $input
echo "  T0 = $T0" >> $input
echo "  p_ratio = $prat" >> $input
echo "  shock = $shock" >> $input
echo "/" >> $input
echo "" >> $input
echo "&numerical" >> $input
echo "  CFL = $cfl" >> $input
echo "  max_iter = $maxk" >> $input
echo "/" >> $input
echo "" >> $input
echo "&flux" >> $input
echo "  flux_scheme = $flux" >> $input
echo "  limiter_scheme = $limiter" >> $input
echo "  beta_lim = $beta_lim" >> $input
echo "  k2 = $k2" >> $input
echo "  k4 = $k4" >> $input
echo "/" >> $input
echo "" >> $input
echo "&output" >> $input
echo "  soln_save = $Sout" >> $input
echo "  res_save = $Rout" >> $input
echo "  res_out  = $disp_out" >> $input
echo "  cons     = $cons" >> $input
echo "/" >> $input
echo "" >> $input
echo "&reconstruction" >> $input
echo "  epsM = $eps_MUSCL" >> $input
echo "  kappaM = $kappa_MUSCL" >> $input
echo "/" >> $input

./../build/bin/test_program
#echo "N$imax" >> $summary
cat "$temp" >> $summary

  done        
 done
done
rm -f temp.txt
end=`date +%s`
runtime=$((end-start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"
