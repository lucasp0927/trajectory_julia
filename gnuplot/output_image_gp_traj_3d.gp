set terminal png size 1500,1000 enhanced font 'Verdana,10'

set output "data.png"
set multiplot layout 2,2 columnsfirst
set datafile separator ","
set style line 1 lc rgb "red"
set style line 2 lc rgb "black"
set xrange [xstart:xend]
set yrange [ystart:yend]
set zrange [zstart:zend]
#set key off
#set view map
#set parametric
set xlabel "x (nm)"
set ylabel "y (nm)"
set zlabel "z (nm)"

set palette model CMY rgbformulae 7,5,15
#set dgrid3d 30,30
#set hidden3d
set ticslevel 0.0
set label sprintf("t = %4.4f us",time) at (xstart+(xend-xstart)/10),(ystart+(yend-ystart)/10) front
#splot "dots.txt" u 1:2:3:($3>0.5?1:2) w points pointtype 7 pointsize 1 lc variable

set view equal xyz
splot "dots.txt" u 1:2:3 w points pointtype 7 pointsize 1
set view 0,0
splot "dots.txt" u 1:2:3 w points pointtype 7 pointsize 1
set view 90,0
splot "dots.txt" u 1:2:3 w points pointtype 7 pointsize 1
set view 90,90
splot "dots.txt" u 1:2:3 w points pointtype 7 pointsize 1
