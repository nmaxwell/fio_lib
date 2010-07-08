
#rm *.o
#rm *.a




# osx:
# i386, x86_64

g++  -Wall -c -fno-common -arch i386 -O3 exact_fio.cpp -o exact_fio.o
gcc -dynamiclib  -arch i386 -current_version 1.0  exact_fio.o  -o libexact_fio.dylib 
mv libexact_fio.dylib /usr/local/lib/libexact_fio.dylib

gcc -x c -std=c99 -Wall -c -fno-common -arch i386 -O3 bfio_interpolation.cc -o bfio_interpolation.o
gcc -dynamiclib  -arch i386 -current_version 1.0  bfio_interpolation.o  -o libbfio_interpolation.dylib 
mv libbfio_interpolation.dylib /usr/local/lib/libbfio_interpolation.dylib

gcc -x c -std=c99 -Wall -c -fno-common -arch i386 -O3 bfio_prototype.cc -o bfio_prototype.o
gcc -dynamiclib  -arch i386 -current_version 1.0  bfio_prototype.o  -o libbfio_prototype.dylib  -lbfio_interpolation
mv libbfio_prototype.dylib /usr/local/lib/libbfio_prototype.dylib







rm *.o
##rm *.out
