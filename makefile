all:SimGoldbeter
SimGoldbeter:	
	gfortran -o 3Comp1DimCellMov main-TwoKe.f \
	anfang-TwoKe.f \
	rs-TwoTypesKe.f \
	out.f ODE-Merson.f ic-Source.f StartingTime.f \
	-L/usr/bmp/pgplot-5.2/ -lpgplot \
	-L/usr/bmp/slatec-4.1/lib -lslatec \
	-L/usr/bmp/lapack-3.4.0 \
	-lX11 \
	-lpng
