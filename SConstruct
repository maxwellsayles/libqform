from glob import glob
import os

ccflags=["-Wall", "-Werror", "-DNDEBUG", "-O3"]

uname = os.uname()
if uname[0] == 'Darwin' and uname[4] == 'i386':
	ccflags.append('-mdynamic-no-pic')

source_files = glob('*.c')
source_files += glob('dbreps/*.c')

StaticLibrary(target='qform',
              source=source_files,
              CPPPATH=['..'],
              CCFLAGS=ccflags)

SConscript(dirs = ['tests', 'timing'])


