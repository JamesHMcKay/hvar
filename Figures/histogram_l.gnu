unset key
set xlabel 'l'
set ylabel 'Frequency'
set term postscript monochrome enhanced "Times-New-Roman" 16
set output "../Figures/Figures/hist_l.eps"
set xrange[0:360]
binwidth=0.1
set boxwidth binwidth
width=pi;
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot '../Figures/data/Figure_8_hist_l.txt' using (bin($1,binwidth)):(360) smooth freq with boxes
#set terminal epslatex color
#set output "hist_l.tex"
#replot

