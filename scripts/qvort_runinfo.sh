#!/bin/bash
echo -E '\title{Qvort run info}' > run_info.tex
echo -E '\date{\today}' >> run_info.tex
echo -E '\documentclass[12pt]{article}' >> run_info.tex
echo -E '\usepackage{alltt}' >> run_info.tex
echo -E '\usepackage{color}' >> run_info.tex
echo -E '\begin{document}' >> run_info.tex
echo -E '\maketitle' >> run_info.tex
echo -E 'number of timesteps: ' >> run_info.tex
awk -F"nsteps" '{print $2}' run.in | sed '/^$/d' >>run_info.tex
echo -E ' ' >> run_info.tex
echo -E 'particle separation ($\delta$):' >> run_info.tex
awk -F"delta" '{print $2}' run.in | sed '/^$/d' >>run_info.tex
echo -E 'cm' >> run_info.tex
echo -E ' ' >> run_info.tex
echo -E 'timestep ($\delta t$):' >> run_info.tex
awk -F"dt" '{print $2}' run.in | sed '/^$/d' >>run_info.tex
echo -E 's' >> run_info.tex
echo -E ' ' >> run_info.tex
if grep -q corea ./run.in ; then 
  echo -E '\centerline{\bf Vortex information}' >> run_info.tex
  echo -E 'core size ($a$):' >> run_info.tex
  awk -F"corea" '{print $2}' run.in | sed '/^$/d' >>run_info.tex
  echo -E 'cm' >> run_info.tex
  echo -E ' ' >> run_info.tex
fi
if grep -q quant_circ ./run.in ; then 
  echo -E 'quantum of circulation ($\Gamma$):' >> run_info.tex
  awk -F"quant_circ" '{print $2}' run.in | sed '/^$/d' >>run_info.tex
  echo -E '${\rm cm}^2/{\rm s}$' >> run_info.tex
fi
if grep -q box_size ./run.in ; then
  echo -E 'box size ($D$):' >> run_info.tex
  awk -F"box_size" '{print $2}' run.in | sed '/^$/d' >>run_info.tex
  echo -E 'cm' >> run_info.tex
  echo -E ' ' >> run_info.tex
fi
echo -E '\centerline{\bf Initial conditions}' >> run_info.tex
if grep -q hayder_wave ./run.in ; then
  echo -E '5 waves with prespecified amplitudes, multiplied by prefactor ($A$) below:' >> run_info.tex
  echo -E ' ' >> run_info.tex
fi
if grep -q wave_amp ./run.in ; then
  echo -E 'wave amplitude ($A$):' >> run_info.tex
  awk -F"wave_amp" '{print $2}' run.in | sed '/^$/d' >>run_info.tex
  echo -E '$\delta$ cm' >> run_info.tex
  echo -E ' ' >> run_info.tex
fi
if grep -q wave_start ./run.in ; then
  echo -E 'normalised wave number ($k$):' >> run_info.tex
  awk -F"wave_start" '{print $2}' run.in | sed '/^$/d' >>run_info.tex
  echo -E ' ' >> run_info.tex
fi
echo -E ' ' >> run_info.tex
echo -E '\centerline{\bf Velocity}' >> run_info.tex
echo -E 'velocity field:' >> run_info.tex
awk -F"velocity" '{print $2}' run.in | sed '/^$/d' >>run_info.tex
echo -E ', critical opening angle, $\theta_{\rm max}$:' >> run_info.tex
awk -F"tree_theta" '{print $2}' run.in | sed '/^$/d' >>run_info.tex
echo -E ' ' >> run_info.tex
echo -E 'boundary conditions:' >> run_info.tex
awk -F"boundary" '{print $2}' run.in | sed '/^$/d' >>run_info.tex
echo -E ' ' >> run_info.tex
echo -E '\end{document}' >> run_info.tex
latex run_info.tex
dvipdf run_info.dvi 
rm run_info.tex run_info.dvi run_info.log run_info.aux
