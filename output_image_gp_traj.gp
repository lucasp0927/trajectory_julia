set terminal pngcairo size 1300,1000 enhanced font 'Verdana,10'
set output "data.png"
set datafile separator ","
#unset colorbox
set xrange [xstart:xend]
set yrange [ystart:yend]
set key off
set view map
set parametric
set label sprintf("t = %4.4f us",time) at (xstart+(xend-xstart)/10),(ystart+(yend-ystart)/10) front
splot "data.txt" u 1:2:3 w pm3d,"dots.txt" u 1:2:3 w points pointtype 7 pointsize 1
#splot "data.txt" matrix u 1:2:3,"dots.txt" u 1:2:3 w points pointtype 7 pointsize 1
