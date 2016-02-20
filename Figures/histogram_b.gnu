unset key
set xlabel 'b'
set ylabel 'Frequency'
set term postscript monochrome enhanced "Times-New-Roman" 16
set output "../Figures/Figures/hist_b.eps"
set xrange[-40:40]
binwidth=0.1
set boxwidth binwidth
width=pi;
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot '../Figures/data/Figure_8_hist_b.txt' using (bin($1,binwidth)):(180) smooth freq with boxes


