# Makefile file for PyNEST
#
# author: Moritz Helias, november 2005
# moritz@bccn.uni-freiburg.de
# ported for prerelease 1.9.0, Moritz Helias, june 2006
#
# Extended and partly rewritten at the CodeJam 2007
# Jochen Martin Eppler, Marc-Oliver Gewaltig, Moritz Helias
#
# Rewritten to make use of Python's distutils package, 01/31/08
# Sven Schrader, Jochen Martin Eppler
#

defs= @SLI_THREAD_DEFS@ $defs
MAKEFLAGS= @MAKE_FLAGS@

EXTRA_DIST= setup.py.in\
	do_tests.py.in \
	nest examples \
	LICENSE.PyNEST\
	README.PyNEST \
	cynestkernel.h\
	cynestkernel.cpp \
	checkpython.py \
	pyobjectdatum.h \
	cython_neuron.h \
	cython_neuron.cpp \
	datumtopythonconverter.h \
	datumtopythonconverter.cpp\
	kernel.pyx\
	classes.pxd


# 'all' does sometimes miss changes, therefore we 'clean' first
all:	clean
	@PYTHON_EXEC@ setup.py build

install-exec-hook:
	@PYTHON_EXEC@ setup.py install --prefix=$(DESTDIR)@PYNEST_PREFIX@
	if test `echo @PYTHON_VERSION@ | cut -d. -f2` -ge 6; then \
	  mkdir -p $(DESTDIR)@PYNEST_PREFIX@/lib/python@PYTHON_VERSION@/dist-packages; \
	  cd $(DESTDIR)@PYNEST_PREFIX@/lib/python@PYTHON_VERSION@/dist-packages; \
	     ln -f -s ../site-packages/nest; \
	  cd $(DESTDIR)@PYNEST_PREFIX@/lib/python@PYTHON_VERSION@/dist-packages; \
	     ln -f -s ../site-packages/PyNEST*.egg-info; \
	fi

install-data-hook: install-exec
	@INSTALL_PROGRAM@ -m 644 @PKGSRCDIR@/cynest/README.PyNEST $(DESTDIR)@PKGDOCDIR@/
	@INSTALL_PROGRAM@ -m 644 @PKGSRCDIR@/cynest/LICENSE.PyNEST $(DESTDIR)@PKGDOCDIR@/
	mkdir -p $(DESTDIR)@PKGDATADIR@/extras
	mkdir -p $(DESTDIR)@PYNEST_PREFIX@/share/doc/nest/cythondoc

# Neuron.pyx version
	if test `echo @PYTHON_VERSION@ | cut -d. -f1` -ge 3; then \
		cp $(VPATH)/nest/Neuron3.x.pyx $(DESTDIR)@PYNEST_PREFIX@/include/Neuron.pyx; \
	else \
		cp $(VPATH)/nest/Neuron2.x.pyx $(DESTDIR)@PYNEST_PREFIX@/include/Neuron.pyx; \
	fi
	
# graphics tools compilation and ressources copying
	$(VPATH)/modules/graphics/graphics_simulator/compile.sh "$(VPATH)/modules/graphics/graphics_simulator/"
	mv $(VPATH)/modules/graphics/graphics_simulator/g_simulator $(DESTDIR)@PYNEST_PREFIX@/lib/python@PYTHON_VERSION@/site-packages/cynest/g_simulator
	mkdir -p $(DESTDIR)@PYNEST_PREFIX@/lib/python@PYTHON_VERSION@/site-packages/cynest/ressources
	cp -r $(VPATH)/modules/graphics/graphics_simulator/ressources $(DESTDIR)@PYNEST_PREFIX@/lib/python@PYTHON_VERSION@/site-packages/cynest


	cp $(VPATH)/../doc/How_to_Install_CyNEST.pdf "$(DESTDIR)@PYNEST_PREFIX@/share/doc/nest/CyNEST Installation.pdf"
	cp $(VPATH)/../doc/cythondoc/user_doc/user_doc.pdf "$(DESTDIR)@PYNEST_PREFIX@/share/doc/nest/cythondoc/Custom Neurons - User Documentation.pdf"
	cp $(VPATH)/../doc/cythondoc/developer_doc/developer_doc.pdf "$(DESTDIR)@PYNEST_PREFIX@/share/doc/nest/cythondoc/Custom Neurons - Developer Documentation.pdf"
	@INSTALL_PROGRAM@ do_tests.py $(DESTDIR)@PKGDATADIR@/extras/
	for dir in @SLI_SUBDIRS@; do \
	  $(MAKE) -k -C @PKGBUILDDIR@/$$dir install-pynest-component 2> /dev/null || true ; done

clean:
	rm -rf build

datadir=@PKGDOCDIR@/examples/cynest

data_DATA=\
	examples/balancedneuron.py \
	examples/brunel-alpha-nest.py \
	examples/brunel-alpha-numpy.py \
	examples/brunel-delta-nest.py \
	examples/BrodyHopfield.py \
	examples/brette-gerstner-fig-2c.py \
	examples/brette-gerstner-fig-3d.py \
	examples/CampbellSiegert.py \
	examples/mc_neuron.py \
	examples/multimeter.py \
	examples/one-neuron.py \
	examples/one-neuron-with-noise.py \
	examples/one-neuron-with-sine-wave.py \
	examples/pulsepacket.py \
	examples/testiaf.py \
	examples/tsodyks_facilitating.py \
	examples/tsodyks_depressing.py \
	examples/vinit_example.py
