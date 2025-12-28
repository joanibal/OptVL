# Master makefile for OptVL. The actual makefile you want is:
# src/build/Makefile

default:
# Check if the config.mk file is in the config dir.
	@if [ ! -f "config/config.mk" ]; then \
	echo "WARNING: this method of compiling is deprecated use meson insted!"; \
	echo "Before compiling, copy an existing config file from the "; \
	echo "config/defaults/ directory to the config/ directory and  "; \
	echo "rename to config.mk. For example:"; \
	echo " ";\
	echo "  cp config/defaults/config.LINUX_GFORTRAN.mk config/config.mk"; \
	echo " ";\
	echo "The modify this config file as required. With the config file"; \
	echo " specified, rerun 'make' and the build will start"; \
	else make optvl_build;\
	fi;

clean:
	rm -fr src/build/*.mod
	rm -fr src/build/*.o
	rm -fr src/build/*.a
	rm -fr src/build/*.so
	rm -fr optvl/*.so
	rm -f *~ config.mk;

optvl_build:
	ln -sf config/config.mk config.mk;
	(cd src/build/ && make)
