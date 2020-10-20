/*
BSD 3-Clause License

Copyright (c) 2019, Mihai Francu
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef UTILS_H
#define UTILS_H

#include <stdlib.h>
#include <stdio.h>

#ifdef _WIN32
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif
#include<iostream>

inline std::string GetCurrentWorkingDir(void) {
	char buff[FILENAME_MAX];
	GetCurrentDir(buff, FILENAME_MAX);
	std::string current_working_dir(buff);
	return current_working_dir;
}

//#define DISABLE_PRINTF
#ifndef DISABLE_PRINTF
void Printf(const char* text, ...);
#else
inline void Printf(const char* text, ...) { }
#endif

#ifndef DISABLE_PRINTF
#include <stdarg.h>
#include <stdio.h>
#include "Platform.h"
#if defined(LINUX)
#	define Printf printf
#else
extern FILE* out;
inline void Printf(const char* text, ...)
{
	va_list arg_list;
	const int strSize = 8192;
	static char str[strSize];

	va_start(arg_list, text);
#ifdef WIN32
	vsprintf_s(str, strSize, text, arg_list);
#else
	vsprintf(str, text, arg_list);
#endif
	va_end(arg_list);

#if defined(_WIN32)
	OutputDebugStringA(str);
#endif
	std::cout << str;
	if (out)
		fprintf(out, "%s", str);
#ifdef ANDROID_NDK
	__android_log_print(ANDROID_LOG_INFO, "GameEngine", "%s", str);
#endif
}
#endif
#endif // !LINUX

#if defined(_DEBUG) && !defined(ANDROID_NDK)
//#	include <assert.h>
#	define ASSERT(_condition) if (!(_condition)) __debugbreak(); //DebugBreak(); //assert(_condition)
#else
#	define ASSERT(_condition)
#endif

#define ALIGN16 __declspec(align(16))

#ifdef ANDROID_NDK
#	include <android/log.h>
#	define HWND int
#	define sprintf_s(a, b, c, ...) sprintf(a, c, ## __VA_ARGS__)

struct POINT
{
	int x, y;
};

#	define VK_UP 201
#	define VK_DOWN 202
#	define VK_LEFT 203
#	define VK_RIGHT 204
#	define VK_LBUTTON 205
#endif // ANDROID_NDK

#ifdef ANDROID_NDK
#	define INLINE inline
#else
#	define INLINE __forceinline
#endif

inline unsigned int GetRandomColor(int offset = 60)
{
	int modulo = 256 - offset;
	int red = offset + rand() % modulo;
	int green = offset + rand() % modulo;
	int blue = offset + rand() % modulo;
	return 0xff000000 | (red & 0xff) << 16 | (green & 0xff) << 8 | (blue & 0xff);
}

// TODO: move to math?
inline float GetRandomReal01()
{
	const int r = rand();
	return (float)r / (float)RAND_MAX;
}

inline float GetRandomReal11()
{
	return -1.0f + 2.0f * GetRandomReal01();
}

#endif // UTILS_H
