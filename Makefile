# The "make" targets are:
# 	wheel: build a Python wheel in "dist" directory.
# 	app-install: build wheel (if needed) and install in ChimeraX.
# 	test: run ChimeraX
# 	debug: run ChimeraX with debugging flag set
# 	clean: remove files used in building wheel
# 	distclean: remove files used in building wheel and license file

# These parameters may be changed as needed.

# ChimeraX bundle names must start with "ChimeraX_"
# to avoid clashes with package names in pypi.python.org.
# When uploaded to the ChimeraX toolshed, the bundle
# will be displayed without the ChimeraX- prefix.
# BUNDLE_NAME = ChimeraX-privateer
# BUNDLE_VERSION = 0.1
# ChimeraX bundles should only include packages
# that install as chimerax.package_name.
# General Python packages should be uploaded to
# pypi.python.org rather than the ChimeraX toolshed.
# PKG_NAME = chimerax.privateer

# Define where ChimeraX is installed.

ifeq ($(OS),Windows_NT)     # is Windows_NT on XP, 2000, 7, Vista, 10...
detected_OS := Windows
else
detected_OS := $(shell uname -s)  # same as "uname -s" but seems to return "Linux  " rather than just "Linux"
endif


ifeq ($(detected_OS),Windows)
# Windows
CHIMERAX_APP = "/c/Program Files/ChimeraX"
endif

ifeq ($(detected_OS),Darwin)
# Mac
CHIMERAX_APP = /Applications/ChimeraX.app
endif


ifeq ($(detected_OS),Linux  ) # For some reason the OS name has trailing white space
CHIMERAX_APP = chimerax
endif


# ==================================================================
# Theoretically, no changes are needed below this line

# Platform-dependent settings.  Should not need fixing.
# For Windows, we assume Cygmakewin is being used.
ifeq ($(detected_OS),Windows)
CHIMERAX_EXE = $(CHIMERAX_APP)/bin/ChimeraX.exe
endif
ifeq ($(detected_OS),Darwin)
CHIMERAX_EXE = $(CHIMERAX_APP)/Contents/bin/ChimeraX
export MACOSX_DEPLOYMENT_TARGET=10.13
endif
ifeq ($(detected_OS),Linux  ) # For some reason the OS name has trailing white space
CHIMERAX_EXE = $(CHIMERAX_APP)
endif

BUNDLE_BASE_NAME = $(subst ChimeraX-,,$(BUNDLE_NAME))
SOURCE = src
SRCS = $(SOURCE)/*.py #$(SOURCE)/*.cpp


:DEFAULT_GOAL := wheel

#
# Actual make dependencies!
#

wheel $(WHEEL): bundle_info.xml $(SRCS)
	$(CHIMERAX_EXE) --nogui --safemode --cmd "devel build . ; exit"

install app-install:	$(WHEEL)
	$(CHIMERAX_EXE) --nogui --safemode --cmd "devel install . ; exit"

uninstall app-uninstall:	$(WHEEL)
	$(CHIMERAX_EXE) --nogui --safemode --cmd "toolshed uninstall $(BUNDLE_BASE_NAME) ; exit"

docs:
	$(CHIMERAX_EXE) -m sphinx docs/source src/docs/user

test:
	$(CHIMERAX_EXE)

debug:
	$(CHIMERAX_EXE) --debug

clean:
	$(CHIMERAX_EXE) --nogui --safemode --cmd "devel clean . ; exit"

.PHONY: docs
