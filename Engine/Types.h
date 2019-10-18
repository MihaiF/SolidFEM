#ifndef TYPES_H
#define TYPES_H

// NB: these are platform/processor specific
typedef unsigned char uint8;
typedef signed char int8;
typedef signed short int16;
typedef unsigned short uint16;
typedef signed int int32;
typedef unsigned int uint32;
typedef signed long long int64;
typedef unsigned long long uint64;
typedef char char8;
typedef wchar_t char16;

typedef unsigned int uint;

// 2^(n-1)-1
#define int16_MAX		0x7fff
// -2^(n-1) or 1 - 2^(n-1)
#define int16_MIN		0xffff

#endif // TYPES_H