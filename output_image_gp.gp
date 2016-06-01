set terminal pngcairo size 1300,1000 enhanced font 'Verdana,10'
set output "data.png"
set datafile separator ","
set xrange [xstart:xend]
set yrange [ystart:yend]
set key off
set view map
splot "data.txt" u 1:2:3 w pm3d
