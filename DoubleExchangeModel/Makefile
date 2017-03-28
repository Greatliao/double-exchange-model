all: double_exchange_mode double_exchange_mode.f90
	#ifort -w -r8 double_exchange_mode.f90 -o double_exchange_mode -L/opt/intel/composer_xe_2013/mkl/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
	ifort  double_exchange_mode.f90

clean:
	rm dem
