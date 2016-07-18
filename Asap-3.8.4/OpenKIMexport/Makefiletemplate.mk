#
# CDDL HEADER START
#
# The contents of this file are subject to the terms of the Common Development
# and Distribution License Version 1.0 (the "License").
#
# You can obtain a copy of the license at
# http://www.opensource.org/licenses/CDDL-1.0.  See the License for the
# specific language governing permissions and limitations under the License.
#
# When distributing Covered Code, include this CDDL HEADER in each file and
# include the License file in a prominent location with the name LICENSE.CDDL.
# If applicable, add the following below this CDDL HEADER, with the fields
# enclosed by brackets "[]" replaced with your own identifying information:
#
# Portions Copyright (c) [yyyy] [name of copyright owner]. All rights reserved.
#
# CDDL HEADER END
#

#
# Copyright (c) 2012, Regents of the University of Minnesota.  All rights reserved.
#
# Contributors:
#    Ellad B. Tadmor
#    Ryan S. Elliott
#    Jakob Schiotz


# load all basic KIM make configuration
include ../Makefile.KIM_Config

# Set all model driver specific information in the include file, so that this
# makefile can be common for all Asap-derived drivers (OK, only one yet, but more
# to come!)
# The include file sets MODEL_DRIVER_NAME and MODEL_SOURCES
include Modeldefinition.mk

LOCALOBJ = $(patsubst %.cpp,%.o,$(MODEL_SOURCES))

include Depend.mk

# APPEND to compiler option flag lists
#FFLAGS   +=
#CFLAGS  +=
#CXXFLAGS +=
#LDFLAGS  +=

# load remaining KIM make configuration
include $(KIM_DIR)/$(builddir)/Makefile.ModelDriver
