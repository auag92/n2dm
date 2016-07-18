# If not provided, use default distribution tag
%{!?disttag: %define disttag %{dist}}

# Set niflheim to 0 if we are not on niflheim
%{!?niflheim: %define niflheim 0}

# The following options are supported:
#   --with prefix=prefix # location of root
#   --with compiler=compiler
#   --with openmpi=openmpi
#   --with openmpidir=openmpidir
#   --with[out] default_version - makes this version the default
#   --with[out] parallel
# with the following defaults:
%define prefix /
%define compiler gfortran
%define openmpi openmpi
%define openmpidir %{_prefix}
%define default_version 0
%define parallel 0

%{?_with_prefix:%define prefix %(set -- %{_with_prefix}; echo $1 | grep -v with | sed 's/=//')}
%{?_with_compiler:%define compiler %(set -- %{_with_compiler}; echo $1 | grep -v with | sed 's/=//')}
%{?_with_openmpi:%define openmpi %(set -- %{_with_openmpi}; echo $1 | grep -v with | sed 's/=//')}
%{?_with_openmpidir:%define openmpidir %(set -- %{_with_openmpidir}; echo $1 | grep -v with | sed 's/=//')}

# --with/--without processing
# first, error if conflicting options are used
%{?_with_default_version: %{?_without_default_version: %{error: both _with_default_version and _without_default_version}}}
%{?_with_parallel: %{?_without_parallel: %{error: both _with_parallel and _without_parallel}}}

# did we find any --with options?
%{?_with_default_version: %define default_version 1}
%{?_with_parallel: %define parallel 1}

# did we find any --without options?
%{?_without_default_version: %define default_version 0}
%{?_without_parallel: %define parallel 0}

%define serial_version serial_version
%if ! %{parallel}
%define openmpidir %{serial_version}
%endif

%ifarch ppc64 sparc64 x86_64 ia64
%define bitness 64
%define is64bit 1
%else
%define bitness 32
%define is64bit 0
%endif

%define real_name Asap
%define ppkg_name %{real_name}
%define ppkg_release 1

Summary: Asap is a tool for doing atomic-scale computer simulations
Name: %{ppkg_name}
Version: 2.20.0
Release: %{ppkg_release}.%{disttag}
License: GPL
Group: Applications/Engineering
#URL: http://www.camp.dtu.dk/Software
URL: http://dcwww.fysik.dtu.dk/software/
Source: %{real_name}-%{version}.tar.gz
BuildRoot: %{_tmppath}/%{name}-%{version}-%{release}-root
Provides: %{ppkg_name}
Requires: campos-ase
Obsoletes: Asap-parallel
AutoReq: 0
####################################
# deal with compiler
%if %{compiler} == none
%define compiler "no compiler used"
%define compiler_string none
%else
%if %{compiler} == gfortran
Requires: libgfortran >= 4.1.1
%define compiler_string gnu
%else
%if %{compiler} == pgf90
%{error: %{compiler} compiler unsupported}
%quit
%define compiler_string pgi
%if %{niflheim}
Requires: pgi-compat = 6.2.4
%endif
%else
%if %{compiler} == pathf90
%define compiler_string pathscale
%if %{niflheim}
Requires: pathscale-compat = 2.5
%endif
%else
%if %{compiler} == ifort
%define compiler_string intel
%else
%{error: %{compiler} compiler unsupported}
%quit
%endif
%endif
%endif
%endif
%endif
# end of: deal with compiler
####################################
# deal with openmpi
%if %{parallel}
%if %{openmpi} == lam
%{error: %{openmpi} openmpi unsupported}
%quit
# Lam before 7.1.1-5 is missing:
# -shared library support
# -fPIC compilation flag
BuildRequires: lam-devel >= 2:7.1.1-5
Requires: lam >= 2:7.1.1-5
%else
%if %{openmpi} == openmpi
BuildRequires: python-scientific
BuildRequires: %{openmpidir}/bin/mpipython openmpi >= 1.2.3-1.%{compiler}.%{disttag}
Requires: python-scientific
Requires: %{openmpidir}/bin/mpipython openmpi >= 1.2.3-1.%{compiler}.%{disttag}
%else
%{error: %{openmpi} openmpi unsupported}
%quit
%endif
%endif
%endif
# end of: deal with openmpi
####################################

%description
Atomic SimulAtion Program - As Soon As Possible.
Asap is a tool for doing atomic-scale computer simulations (mainly
molecular dynamics) using classical potentials (mainly Effective
Medium Theory).
Built with:
%if %{parallel}
%{openmpidir}/bin/mpipython,
%endif
%{compiler_string}.

Author:  Jakob Schiotz et al. 

%prep
%setup -n %{real_name}-%{version}

%build

%define __find_requires %{nil}

#CFLAGS="%{optflags}"; export CFLAGS
#CXXFLAGS="%{optflags}"; export CXXFLAGS
####################################
# deal with openmpi
%if %{parallel}
%if %{openmpi} == openmpi
. %{openmpidir}/bin/mpivars-*.sh
which mpicc
%else
%{error: %{openmpi} openmpi unsupported}
%quit
%endif
%endif
# end of: deal with openmpi
####################################
# deal with compiler
%if %{compiler} == none
compiler_string=none
%else
%if %{compiler} == gfortran
compiler_string=gnu
%else
%if %{compiler} == pgf90
%{error: %{compiler} compiler unsupported}
%quit
compiler_string=pgi
%else
%if %{compiler} == pathf90
compiler_string=pathscale
%else
%if %{compiler} == ifort
compiler_string=intel
%else
%{error: %{compiler} compiler unsupported}
%quit
%endif
%endif
%endif
%endif
%endif
# end of: deal with compiler
####################################
# detect the compiler
compiler_detected=`python detect.py compiler`
compiler_matches=`python -c "print '${compiler_detected}' == '${compiler_string}'"`
if [ "${compiler_matches}" == "True" ]
then
make depend
%if %{parallel}
make all
%else
make serial
%endif
else
%{error: compiler detected does not match the one requested}
%quit
fi

%install
%define __find_requires %{nil}
%{__rm} -rf %{buildroot}
DESTDIR="%{buildroot}/%{prefix}"; export DESTDIR
%if %{parallel}
make install-parallel
%else
make install
%endif

%if %{parallel}
%define bindir %{prefix}/usr/bin
mv %{buildroot}/%{prefix}/usr/local/bin %{buildroot}/%{prefix}/usr/local/bin.bak
mkdir -p %{buildroot}/%{bindir}
mv %{buildroot}/%{prefix}/usr/local/bin.bak/asap-niflheim %{buildroot}/%{bindir}/%{ppkg_name}-%{version}.%{ppkg_release}
rm -rf %{buildroot}/%{prefix}/usr/local/bin.bak
# make it the default version
%if %{default_version}
ln -s %{bindir}/%{ppkg_name}-%{version}.%{ppkg_release} %{buildroot}/%{bindir}/asap-niflheim
%endif
%endif


%clean
%{__rm} -rf %{buildroot}


%files
%defattr(-, root, root, 0755)
# %doc doc/ LICENSE PKG-INFO README.txt
# %{_bindir}/*
%{_libdir}/python*/site-packages/%{real_name}
%{_libdir}/python*/site-packages/asapserial2.py
%{_libdir}/python*/site-packages/_asapserial2.so
%if %{parallel}
%{_libdir}/python*/site-packages/asapparallel2.py
%{_libdir}/python*/site-packages/_asapparallel2.so
%if %{default_version}
%{bindir}/asap-niflheim
%endif
%{bindir}/%{ppkg_name}-%{version}.%{ppkg_release}
%endif
