#ifndef PLATFORM_H
#define PLATFORM_H

#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#	define WIN32_LEAN_AND_MEAN             // Exclude rarely-used stuff from Windows headers
#endif
#	include <windows.h>

#undef min
#undef max

#endif

#endif // PLATFORM_H