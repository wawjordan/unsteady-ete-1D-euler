#!./bin/bash
input="input.nml"
summary="../results/summary.dat"
temp="temp.txt"
start=`date +%s`
flux_str=""
limiter_str=""
#touch $summary
#if [ -f $summary ]; then
#  rm -f "$summary"
#fi
imax=512 #16 32 64 128 256 512
p0=300.0
T0=600.0
flux=1
limiter=2
beta_lim=2.0
cfl=1.0 #0.1 0.5 0.9
eps_roe=0.1
eps_MUSCL=1.0
kappa_MUSCL=-1.0
maxk=600
Sout=16
Rout=1
disp_out=100
cons=F
if [ $flux -eq 1 ]; then
  flux_str="van Leer flux"
elif [ $flux -eq 2 ]; then
  flux_str="Roe's flux"
fi

if [ $limiter -eq 1 ]; then
  limiter_str="van Leer limiter"
elif [ $limiter -eq 2 ]; then
  limiter_str="van Albada limiter"
elif [ $limiter -eq 3 ]; then
  limiter_str="minmod limiter"
elif [ $limiter -eq 4 ]; then
  limiter_str="beta ($beta_lim) limiter"
fi
echo ""
echo "!==============================================================================!"
echo "| N = $imax | $flux_str | $limiter_str | CFL=$cfl "
echo "!==============================================================================!"

echo "&grid" > $input
echo "  imax = $imax" >> $input
echo "  xmin = -1.0" >> $input
echo "  xmax =  1.0" >> $input
echo "  n_ghost = 2" >> $input
echo "/" >> $input
echo "" >> $input
echo "&geometry" >> $input
echo "  areaStar = 0.2" >> $input
echo "/" >> $input
echo "" >> $input
echo "&constants" >> $input
echo "  gamma = 1.4" >> $input
echo "/" >> $input
echo "" >> $input
echo "&initial" >> $input
echo "  p0 = $p0" >> $input
echo "  T0 = $T0" >> $input
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
echo "  eps_roe = $eps_roe" >> $input
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
#cat "$temp" >> $summary
rm -f temp.txt
end=`date +%s`
runtime=$((end-start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"
