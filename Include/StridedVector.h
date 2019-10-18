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
