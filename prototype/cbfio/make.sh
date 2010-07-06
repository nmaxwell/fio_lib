
#rm *.o
#rm *.a




# osx:
#gcc -x c -std=c99 -Wall -c -fno-common -arch i386 -O3 exact_fio.cc -o exact_fio.o
#x86_64, i386
#gcc -dynamiclib  -arch i386 -current_version 1.0  exact_fio.o  -o libexact_fio.dylib 
#mv libexact_fio.dylib /usr/local/lib/libexact_fio.dylib




# linux:
#g++ -x c -std=c99 -O3 -fPIC -Wall -c exact_fio.cc -o exact_fio.o
#g++ -shared -Wl,-soname,libexact_fio.so.1 -o libexact_fio.so.1.0  exact_fio.o

#cp libexact_fio.so.1.0 /usr/local/lib/
#cp libexact_fio.so.1.0 /usr/lib/
#ln -sf /usr/local/lib/libexact_fio.so.1.0 /usr/local/lib/libexact_fio.so
#ln -sf /usr/local/lib/libexact_fio.so.1.0 /usr/local/lib/libexact_fio.so.1
#ln -sf /usr/lib/libexact_fio.so.1.0 /usr/lib/libexact_fio.so
#ln -sf /usr/lib/libexact_fio.so.1.0 /usr/lib/libexact_fio.so.1





gcc -x c -std=c99 -O3 -fPIC -Wall -c bfio_interpolation.cc -o bfio_interpolation.o
gcc -shared -Wl,-soname,libbfio_interpolation.so.1 -o libbfio_interpolation.so.1.0  bfio_interpolation.o

ln -sf libbfio_interpolation.so.1.0 libbfio_interpolation.so.1
ln -sf libbfio_interpolation.so.1.0 libbfio_interpolation.so

cp libbfio_interpolation.so.1.0 /usr/local/lib/
cp libbfio_interpolation.so.1.0 /usr/lib/
ln -sf /usr/local/lib/libbfio_interpolation.so.1.0 /usr/local/lib/libbfio_interpolation.so
ln -sf /usr/local/lib/libbfio_interpolation.so.1.0 /usr/local/lib/libbfio_interpolation.so.1
ln -sf /usr/lib/libbfio_interpolation.so.1.0 /usr/lib/libbfio_interpolation.so
ln -sf /usr/lib/libbfio_interpolation.so.1.0 /usr/lib/libbfio_interpolation.so.1






gcc -x c -std=c99 -O3 -fPIC -Wall -c bfio_prototype.cc -o bfio_prototype.o
gcc -shared -Wl,-soname,libbfio_prototype.so.1 -o libbfio_prototype.so.1.0  bfio_prototype.o

ln -sf libbfio_prototype.so.1.0 libexact_fio.so.1
ln -sf libbfio_prototype.so.1.0 libexact_fio.so

cp libbfio_prototype.so.1.0 /usr/local/lib/
cp libbfio_prototype.so.1.0 /usr/lib/
ln -sf /usr/local/lib/libbfio_prototype.so.1.0 /usr/local/lib/libbfio_prototype.so
ln -sf /usr/local/lib/libbfio_prototype.so.1.0 /usr/local/lib/libbfio_prototype.so.1
ln -sf /usr/lib/libbfio_prototype.so.1.0 /usr/lib/libbfio_prototype.so
ln -sf /usr/lib/libbfio_prototype.so.1.0 /usr/lib/libbfio_prototype.so.1







rm *.o
rm *.so
rm *.so.1
rm *.so.1.0
##rm *.out
