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

#ifndef STRIDED_VECTOR_H
#define STRIDED_VECTOR_H

// Custom references
#include <Engine/Types.h>


template<class T>
class StridedVector {
public:

	StridedVector(T* buffer, uint elements, uint stride)
	{
		// TODO remove 'elements' as it is not really used for anything?
		// -> The buffer might even change in size after the stridedvector has been initialized

		// Assert that buffer is not pointing to null
		mBuffer = buffer;
		mElements = elements;
		mStride = stride;
	}

	~StridedVector()
	{
		// Do nothing
	}

	T* GetAt(uint element)
	{
		// Assert that element < mElements

		char* ptr = (char*)mBuffer + element * mStride;
		T* pVal = (T*)ptr;
		return pVal;
	}

	T* mBuffer;
	uint mElements;
	uint mStride;

private:

};

#endif // STRIDED_VECTOR_H
