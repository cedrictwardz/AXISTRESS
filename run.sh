#make axitra
#gfortran cut.f90 -o cut
#clean everything
make clean
rm funcAmp.txt cut
rm axi.x axi.y axi.z axi.e axi.n axi.head axi.sis axi.res
rm st*bin
#cp axi.data1 axi.data
#./axitra
#cp axi.data2 axi.data
#./axitra

#./cut
