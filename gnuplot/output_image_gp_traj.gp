set terminal png size 1500,1000 enhanced font 'Verdana,10'
set output "data.png"
set datafile separator ","
set style line 1 lc rgb "red"
set style line 2 lc rgb "black"
set xrange [xstart:xend]
set yrange [ystart:yend]
set key off
set view map
set parametric
set xlabel "x (nm)"
set ylabel "y (nm)"
set cblabel "U (K)"
set label sprintf("t = %4.4f us",time) at (xstart+(xend-xstart)/10),(ystart+(yend-ystart)/10) front
splot "data.txt" u 1:2:3 w pm3d,"dots.txt" u 1:2:3:($3>0.5?1:2) w points pointtype 7 pointsize 1 lc variable
