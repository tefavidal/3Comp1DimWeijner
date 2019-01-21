all:SimGoldbeter
SimGoldbeter:	
	gfortran -o 3Comp1DimCellMov main-LowKe.f \
	anfang-LowKe.f \
	rs-LowKeWithSource.f \
	CellMovement.f \
	Lap-NoFlux.f \
	out.f ODE-Merson.f ic-Source.f StartingTime.f \
