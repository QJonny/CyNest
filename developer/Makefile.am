# Automake file for the Developer Module
#
# Hans Ekkehard Plesser, June 2006 based on
# Jochen Martin Eppler, February 2004 based on
# Marc-Oliver Gewaltig, July 2000
#
# To avoid double work, developermodule is built only 
# for compile-time linking, not for dynamic loading 
# 
# Logically, libdevelopermodule is a convenience library,
# but when using the noinst_ prefix, libdevelopermodule
# is built as .a only, and module registration does not
# work in the final executable. Therefore, the library
# is built as a normal library.
# 
# Headers from this directory are not to be installed upon
# make install. They are therefore included in _SOURCES.
#
# IMPORTANT --- IMPORTANT --- IMPORTANT --- IMPORTANT
#
# After changes to developer/Makefile.am (including changes
# made by svn update), you MUST run ./bootstrap.sh in the
# toplevel source directory. Otherwise, new sources in the
# developer directory will not be compiled. See also #211.
#
# IMPORTANT --- IMPORTANT --- IMPORTANT --- IMPORTANT

defs= @SLI_THREAD_DEFS@ $defs
MAKEFLAGS= @MAKE_FLAGS@

noinst_LTLIBRARIES= libdevelopermodule.la

libdevelopermodule_la_CXXFLAGS= @AM_CXXFLAGS@
libdevelopermodule_la_SOURCES= \
			developermodule.cpp developermodule.h \
			a2eif_cond_exp.cpp a2eif_cond_exp.h \
			a2eif_cond_exp_HW.cpp a2eif_cond_exp_HW.h \
			spike_queue.cpp spike_queue.h\
			selective_connector.h\
			markram_connection.cpp markram_connection.h \
			schemmel_s1_connection.cpp schemmel_s1_connection.h \
			maturing_connection.cpp maturing_connection.h\
			maturing_connector.cpp maturing_connector.h\
			maturing_connection_fr.cpp maturing_connection_fr.h\
			maturing_connector_fr.cpp maturing_connector_fr.h\
			stdp_triplet_connection.h stdp_triplet_connection.cpp\
			annealing_connection.h annealing_connection.cpp\
			iaf_chxk_2008.cpp iaf_chxk_2008.h\
			iaf_chs_2007.cpp iaf_chs_2007.h\
			environment.cpp environment.h\
			environment_cliff.cpp environment_cliff.h\
                        theta_neuron.cpp theta_neuron.h\
                        theta_neuron_ps.cpp theta_neuron_ps.h\
			static_connection_mult0.h\
                        iaf_psc_delta_canon_cvv.cpp iaf_psc_delta_canon_cvv.h\
			iaf_psc_delta_nodelay.cpp iaf_psc_delta_nodelay.h\
                        stdp_connection_rsnn_spikepairing_hom.cpp stdp_connection_rsnn_spikepairing_hom.h\
			pif_psc_delta_canon_cvv.cpp pif_psc_delta_canon_cvv.h\
	 		inh_gamma_generator.cpp inh_gamma_generator.h\
			ac_gamma_generator.cpp ac_gamma_generator.h\
			ht_neuron_fs.cpp ht_neuron_fs.h\
			iaf_psc_alpha_mod.cpp iaf_psc_alpha_mod.h\
			iaf_psc_alpha_norec.cpp iaf_psc_alpha_norec.h\
		        iaf_psc_alpha_dynth.cpp iaf_psc_alpha_dynth.h\
			spike_dilutor.cpp spike_dilutor.h\
			lossy_connection.cpp lossy_connection.h\
                        individual_spike_generator.cpp individual_spike_generator.h\
		        tsodyks2_stdp_doublet_connection.h tsodyks2_stdp_doublet_connection.cpp



# Please leave this separating blank line. It is needed for the build tools on the PETA fx machine to work. MH 2010-11-4
#			dopamine_connection.cpp dopamine_connection.h\
#			policy_connection.cpp policy_connection.h\
#			dopamine_th_connection.cpp dopamine_th_connection.h\
#			policy_th_connection.cpp policy_th_connection.h\
#			iaf_psc_exp_canon.cpp iaf_psc_exp_canon.h
#	 		pp_psc_delta.cpp pp_psc_delta.h\
#			iaf_psc_delta_canon_stepcurr.cpp iaf_psc_delta_canon_stepcurr.h\
#			iaf_psc_alpha_currents.cpp iaf_psc_alpha_currents.h\
#			iaf_psc_delta_currents.cpp iaf_psc_delta_currents.h\
#			iaf_psc_alpha_canon_nr.cpp iaf_psc_alpha_canon_nr.h\
# 			ac_poisson_generator.cpp ac_poisson_generator.h\
#	 		ginzburg_neuron.cpp ginzburg_neuron.h\
#	 		iaf_cond_delta.cpp iaf_cond_delta.h



AM_CPPFLAGS= -I$(top_srcdir)/libnestutil\
	  -I$(top_srcdir)/librandom\
	  -I$(top_srcdir)/sli\
	  -I$(top_srcdir)/nestkernel\
	  -I$(top_srcdir)/models\
	  @INCLTDL@ \
          @GSL_CFLAGS@ @MUSIC_INCLUDE@ @MPI_INCLUDE@

nobase_pkgdata_DATA= sli/developer-init.sli

# Documentation extraction on the developer module was never properly
# called. Since the developer module is going away at some point, it
# doesn't make sense to fix it. Leaving this in place in case if this
# decision is revisited (see #394).

HELPDIRS="@PKGSRCDIR@/$(subdir):@PKGSRCDIR@/$(subdir)/sli:@PKGSRCDIR@/$(subdir)/unittests"

install-slidoc:
	mkdir -p $(DESTDIR)@PKGDOCDIR@/help
	$(DESTDIR)$(exec_prefix)/bin/sli --userargs=$(HELPDIRS) $(top_srcdir)/lib/sli/install-help.sli 2>&1 >> @INSTALL_HELP_LOG@

install-data-hook:
	mkdir -p $(DESTDIR)@PKGDOCDIR@/unittests
	cp @PKGSRCDIR@/$(subdir)/unittests/*.sli $(DESTDIR)@PKGDOCDIR@/unittests

EXTRA_DIST= unittests examples sli
