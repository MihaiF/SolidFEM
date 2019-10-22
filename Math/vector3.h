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

#ifndef VECTOR3_H
#define VECTOR3_H

#include <math.h>

namespace Math {

template<typename Real>
struct Vector3T
{
    union 
    {
        struct { Real x, y, z; };
        Real v[3];
    };

	Vector3T()
	{
		SetZero();
	}

	Vector3T(Real a, Real b, Real c)
	{
		v[0] = a;
		v[1] = b;
		v[2] = c;
	}

	explicit Vector3T(Real a)
	{
		v[0] = v[1] = v[2] = a;
	}

	void SetZero()
	{
		v[0] = v[1] = v[2] = 0.f;
	}

	void Set(Real a, Real b, Real c)
	{
		v[0] = a;
		v[1] = b;
		v[2] = c;
	}

	Real operator [](int i) const { return v[i]; }
	Real& operator [](int i) { return v[i]; }

	Real operator ()(int i) const { return v[i]; }
	Real& operator ()(int i) { return v[i]; }

	Real X() const { return v[0]; }
	Real Y() const { return v[1]; }
	Real Z() const { return v[2]; }
	Real& X() { return v[0]; }
	Real& Y() { return v[1]; }
	Real& Z() { return v[2]; }

	void Add(const Vector3T& other)
	{
		for (int i = 0; i < 3; i++)
			v[i] += other[i];
	}

	void Subtract(const Vector3T& other)
	{
		for (int i = 0; i < 3; i++)
			v[i] -= other[i];
	}

	Vector3T& operator +=(const Vector3T& other)
	{
		Add(other);
		return *this;
	}

	Vector3T& operator -=(const Vector3T& other)
	{
		Subtract(other);
		return *this;
	}

	Vector3T& operator *=(Real s)
	{
		Scale(s);
		return *this;
	}

	Real Dot(const Vector3T& v) const
	{
		return X() * v.X() + Y() * v.Y() + Z() * v.Z();
	}

	void Scale(Real s)
	{
		for (int i = 0; i < 3; i++)
			v[i] *= s;
	}

	void Flip()
	{
		for (int i = 0; i < 3; i++)
			v[i] = -v[i];
	}

	Real LengthSquared() const
	{
		return Dot(*this);
	}

	Real Length() const
	{
		return sqrt(LengthSquared());
	}

	void Normalize()
	{
		Real len = Length();
		//ASSERT(len != 0);
		Scale(1.f / len);
	}

	Vector3T Cross(const Vector3T& v)
	{
		return Vector3T(Y() * v.Z() - v.Y() * Z(), -X() * v.Z() + Z() * v.X(), X() * v.Y() - v.X() * Y());
	}

	Vector3T Perpendicular()
	{
		//Vector3T up(1, 0, 0);
		// TODO: check if dot is close to 1
		Vector3T abs = this->GetAbs();
		int max = abs.GetMaxComponent();
		int max1 = (max + 1) % 3;
		Vector3T up(0);
		up[max1] = 1;
		Vector3T r = this->Cross(up);
		r.Normalize();
		return r;
	}

	int GetMaxComponent()
	{
		Real max = v[0];
		int idx = 0;
		for (int i = 1; i < 3; i++)
		{
			if (v[i] > max)
			{
				max = v[i];
				idx = i;
			}
		}
		return idx;
	}

	Vector3T GetAbs()
	{
		return Vector3T(abs(v[0]), abs(v[1]), abs(v[2]));
	}

	Vector3T Multiply(const Vector3T& u)
	{
		return Vector3T(v[0] * u.v[0], v[1] * u.v[1], v[2] * u.v[2]);
	}
};

template<typename Real>
inline Vector3T<Real> operator +(const Vector3T<Real>& a, const Vector3T<Real>& b)
{
	Vector3T<Real> ret;
	for (int i = 0; i < 3; i++)
		ret[i] = a[i] + b[i];
	return ret;
}

template<typename Real>
inline Vector3T<Real> operator -(const Vector3T<Real>& a, const Vector3T<Real>& b)
{
	Vector3T<Real> ret;
	for (int i = 0; i < 3; i++)
		ret[i] = a[i] - b[i];
	return ret;
}

template<typename Real>
inline Vector3T<Real> operator *(float s, const Vector3T<Real>& v)
{
	return Vector3T<Real>(s * v.X(), s * v.Y(), s * v.Z());
}

template<typename Real>
inline Vector3T<Real> operator *(const Vector3T<Real>& v, float s)
{
	return Vector3T<Real>(s * v.X(), s * v.Y(), s * v.Z());
}

template<typename Real>
inline Vector3T<Real> operator *(double s, const Vector3T<Real>& v)
{
	return Vector3T<Real>(s * v.X(), s * v.Y(), s * v.Z());
}

template<typename Real>
inline Vector3T<Real> operator *(const Vector3T<Real>& v, double s)
{
	return Vector3T<Real>(s * v.X(), s * v.Y(), s * v.Z());
}

template<typename Real>
inline Real operator *(const Vector3T<Real>& u, const Vector3T<Real>& v)
{
	return u.Dot(v);
}

template<typename Real>
inline Vector3T<Real> operator -(const Vector3T<Real>& v)
{
	return Vector3T<Real>(-v.X(), -v.Y(), -v.Z());
}

template<typename Real>
inline Vector3T<Real> cross(const Vector3T<Real>& v1, const Vector3T<Real>& v2)
{
	Vector3T<Real> ret = v1;
	return ret.Cross(v2);
}

template<typename Real>
inline Real dot(const Vector3T<Real>& v1, const Vector3T<Real>& v2)
{
    Vector3T<Real> ret = v1;
    return ret.Dot(v2);
}

template<typename Real>
inline Real triple(const Vector3T<Real>& v1, const Vector3T<Real>& v2, const Vector3T<Real>& v3)
{
	return v3.Dot(cross(v1, v2));
}

template<typename Real>
inline Vector3T<Real> vmin(const Vector3T<Real>& v1, const Vector3T<Real>& v2)
{
	return Vector3T<Real>(min(v1.X(), v2.X()), min(v1.Y(), v2.Y()), min(v1.Z(), v2.Z()));
}

template<typename Real>
inline Vector3T<Real> vmax(const Vector3T<Real>& v1, const Vector3T<Real>& v2)
{
	return Vector3T<Real>(max(v1.X(), v2.X()), max(v1.Y(), v2.Y()), max(v1.Z(), v2.Z()));
}

typedef Vector3T<float> Vector3;

inline Vector3 DoubleToFloat(const Vector3T<double>& v)
{
	Vector3 vf;
	vf.x = (float)v.x;
	vf.y = (float)v.y;
	vf.z = (float)v.z;
	return vf;
}

inline Vector3 DoubleToFloat(const Vector3T<float>& v)
{
	return v;
}
} // namespace

#endif // VECTOR3_H

