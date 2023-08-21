#include <vector>
#include <fstream>
#include <strstream>
#include <math.h>

float FastInverseSquareRoot(float x)
{
    /*float f = x;
    long i = *(long *)&x;  // bits from float memory box to long memory box

    // log2(1 + x) ~ bit representation of x
    // i = 1 / 2^23 * (Mantissa_x + 2^23 * Exponent) + error - 127
    // you can shift bits in long
    i = 0x5f3759df - (i >> 1);

    f = *(float *)i;
    f = f * (1.5f - (0.5f * x * f * f));  // Newton's iteration first step
    return f;*/
    union
    {
        float f;
        uint32_t i;
    } conv = {.f = x};
    conv.i = 0x5f3759df - (conv.i >> 1);
    conv.f *= 1.5F - (x * 0.5F * conv.f * conv.f);
    return conv.f;
}

struct vec3d
{
    float x = 0.0f, y = 0.0f, z = 0.0f, w = 1.0f;

    const vec3d operator+(const vec3d rhs)
    {
        vec3d out;
        out.x = x + rhs.x;
        out.y = y + rhs.y;
        out.z = z + rhs.z;
        return out;
    }

    const vec3d operator-(const vec3d rhs)
    {
        vec3d out;
        out.x = x - rhs.x;
        out.y = y - rhs.y;
        out.z = z - rhs.z;
        return out;
    }

    void operator+=(const vec3d rhs)
    {
        this->x += rhs.x;
        this->y += rhs.y;
        this->z += rhs.z;
    }

    void operator-=(const vec3d rhs)
    {
        this->x -= rhs.x;
        this->y -= rhs.y;
        this->z -= rhs.z;
    }

    const vec3d operator*(const float rhs)
    {
        vec3d out;
        out.x = x * rhs;
        out.y = y * rhs;
        out.z = z * rhs;
        return out;
    }

    const vec3d operator/(const float rhs)
    {
        vec3d out;
        out.x = x / rhs;
        out.y = y / rhs;
        out.z = z / rhs;
        return out;
    }

    void operator/=(const float rhs)
    {
        this->x /= rhs;
        this->y /= rhs;
        this->z /= rhs;
    }

    void operator*=(const float rhs)
    {
        this->x *= rhs;
        this->y *= rhs;
        this->z *= rhs;
    }

    float GetLengthSqared()
    {
        return x * x + y * y + z * z;
    }

    float GetLength()
    {
        return sqrtf(x * x + y * y + z * z);
    }

    void Normalize()
    {
        // float reverseL = 1.0f / GetLength();
        *(this) *= FastInverseSquareRoot(GetLengthSqared());
        // this->x *= reverseL;
        // this->y *= reverseL;
        // this->z *= reverseL;
    }

    float DotProduct(vec3d rhs)
    {
        return x * rhs.x + y * rhs.y + z * rhs.z;
    }

    vec3d CrossProduct(vec3d rhs)
    {
        vec3d output;
        output.x = y * rhs.z - z * rhs.y;
        output.y = z * rhs.x - x * rhs.z;
        output.z = x * rhs.y - y * rhs.x;
        return output;
    }
};

struct mat4x4
{
    float m[4][4] = {0.0f};

    vec3d operator*(const vec3d &rhs)
    {
        vec3d output;
        output.x = m[0][0] * rhs.x + m[0][1] * rhs.y + m[0][2] * rhs.z + m[0][3] * rhs.w;
        output.y = m[1][0] * rhs.x + m[1][1] * rhs.y + m[1][2] * rhs.z + m[1][3] * rhs.w;
        output.z = m[2][0] * rhs.x + m[2][1] * rhs.y + m[2][2] * rhs.z + m[2][3] * rhs.w;
        output.w = m[3][0] * rhs.x + m[3][1] * rhs.y + m[3][2] * rhs.z + m[3][3] * rhs.w;
        return output;
    }

    mat4x4 operator*(const mat4x4 &rhs)
    {
        mat4x4 output;
        for (u_int i = 0; i < 4; i++)
            for (u_int j = 0; j < 4; j++)
                for (u_int t = 0; t < 4; t++)
                    output.m[i][j] += m[i][t] * rhs.m[t][j];

        return output;
    }
};

mat4x4 Matrix_MakeIdentity()
{
    mat4x4 output;
    output.m[0][0] = 1.0f;
    output.m[1][1] = 1.0f;
    output.m[2][2] = 1.0f;
    output.m[3][3] = 1.0f;
    return output;
}

mat4x4 Matrix_MakeTranslation(vec3d move)
{
    mat4x4 output;
    output.m[0][0] = 1.0f;
    output.m[1][1] = 1.0f;
    output.m[2][2] = 1.0f;
    output.m[0][3] = move.x;
    output.m[1][3] = move.y;
    output.m[2][3] = move.z;
    output.m[3][3] = 1.0f;
    return output;
}

mat4x4 Matrix_MakeRotationAroundX(float angle)
{
    mat4x4 output;
    output.m[0][0] = 1.0f;
    output.m[1][1] = cosf(angle);
    output.m[2][1] = sinf(angle);
    output.m[2][2] = cosf(angle);
    output.m[1][2] = -sinf(angle);
    output.m[3][3] = 1.0f;
    return output;
}

mat4x4 Matrix_MakeRotationAroundY(float angle)
{
    mat4x4 output;
    output.m[1][1] = 1.0f;
    output.m[0][0] = cosf(angle);
    output.m[0][2] = -sinf(angle);
    output.m[2][2] = cosf(angle);
    output.m[2][0] = sinf(angle);
    output.m[3][3] = 1.0f;
    return output;
}

mat4x4 Matrix_MakeRotationAroundZ(float angle)
{
    mat4x4 output;
    output.m[2][2] = 1.0f;
    output.m[0][0] = cosf(angle);
    output.m[1][0] = sinf(angle);
    output.m[1][1] = cosf(angle);
    output.m[0][1] = -sinf(angle);
    output.m[3][3] = 1.0f;
    return output;
}

mat4x4 Matrix_MakeProjection(float aspectRatio, float FoV, float zNear, float zFar)
{
    mat4x4 output;
    float f = 1.0f / tanf(FoV * 0.5f * M_PI / 180.f);
    float q = zFar / (zFar - zNear);

    output.m[0][0] = aspectRatio * f;
    output.m[1][1] = f;
    output.m[2][2] = q;
    output.m[2][3] = -zNear * q;
    output.m[3][2] = -1.0f;
    output.m[3][3] = 0.0f;
    return output;
}

mat4x4 Matrix_MakePointAt(vec3d pos, vec3d target, vec3d up)
{
    /*// New forward direction
    vec3d newForward = target - pos;
    newForward.Normalize();

    // New upward direction
    vec3d a = newForward * up.DotProduct(newForward);
    vec3d newUp = up - a;
    newUp.Normalize();

    vec3d newRight = newUp.CrossProduct(newForward);

    mat4x4 output;
	output.m[0][0] = newRight.x;	output.m[1][0] = newRight.y;	output.m[2][0] = newRight.z;	output.m[3][0] = 0.0f;
	output.m[0][1] = newUp.x;		output.m[1][1] = newUp.y;		output.m[2][1] = newUp.z;		output.m[3][1] = 0.0f;
	output.m[0][2] = newForward.x;	output.m[1][2] = newForward.y;	output.m[2][2] = newForward.z;	output.m[3][2] = 0.0f;
	output.m[0][3] = pos.x;			output.m[1][3] = pos.y;			output.m[2][3] = pos.z;			output.m[3][3] = 1.0f;
    return output;*/

    // Calculate new forward direction
	vec3d newForward = target - pos;
	newForward.Normalize();

	// Calculate new Up direction
	vec3d a = newForward * up.DotProduct(newForward);
	vec3d newUp = up - a;
	newUp.Normalize();

	// New Right direction is easy, its just cross product
	vec3d newRight = newUp.CrossProduct(newForward);

	// Construct Dimensioning and Translation Matrix	
	mat4x4 matrix;
	matrix.m[0][0] = newRight.x;	matrix.m[1][0] = newRight.y;	matrix.m[2][0] = newRight.z;	matrix.m[3][0] = 0.0f;
	matrix.m[0][1] = newUp.x;		matrix.m[1][1] = newUp.y;		matrix.m[2][1] = newUp.z;		matrix.m[3][1] = 0.0f;
	matrix.m[0][2] = newForward.x;	matrix.m[1][2] = newForward.y;	matrix.m[2][2] = newForward.z;	matrix.m[3][2] = 0.0f;
	matrix.m[0][3] = pos.x;			matrix.m[1][3] = pos.y;			matrix.m[2][3] = pos.z;			matrix.m[3][3] = 1.0f;
	return matrix;
}

mat4x4 Matrix_MakeInverseForRotTrans(mat4x4 m)
{
    /*mat4x4 output;
    output.m[0][0] = input.m[0][0]; output.m[1][0] = input.m[0][1]; output.m[2][0] = input.m[0][2]; output.m[3][0] = 0.0f;
    output.m[0][1] = input.m[1][0]; output.m[1][1] = input.m[1][1]; output.m[2][1] = input.m[1][2]; output.m[3][1] = 0.0f; 
    output.m[0][2] = input.m[2][0]; output.m[1][2] = input.m[2][1]; output.m[2][2] = input.m[2][2]; output.m[3][2] = 0.0f;

    //m03 = -dotProduct(TranslateVector, NewRight)
    //m13 = -dotProduct(TranslateVector, newUp)
    //m23 = -dotProduct(TranslateVector, newForward)
    output.m[0][3] = -(input.m[0][3] * input.m[0][0] + input.m[1][3] * input.m[0][1] + input.m[2][3] * input.m[0][2]);
    output.m[1][3] = -(input.m[0][3] * input.m[1][0] + input.m[1][3] * input.m[1][1] + input.m[2][3] * input.m[1][2]);
    output.m[2][3] = -(input.m[0][3] * input.m[2][0] + input.m[1][3] * input.m[2][1] + input.m[2][3] * input.m[2][2]);

    output.m[3][3] = 1.0f;
    return output;*/

    mat4x4 matrix;
	matrix.m[0][0] = m.m[0][0]; matrix.m[1][0] = m.m[0][1]; matrix.m[2][0] = m.m[0][2]; matrix.m[3][0] = 0.0f;
	matrix.m[0][1] = m.m[1][0]; matrix.m[1][1] = m.m[1][1]; matrix.m[2][1] = m.m[1][2]; matrix.m[3][1] = 0.0f;
	matrix.m[0][2] = m.m[2][0]; matrix.m[1][2] = m.m[2][1]; matrix.m[2][2] = m.m[2][2]; matrix.m[3][2] = 0.0f;
	
    matrix.m[0][3] = -(m.m[0][3] * matrix.m[0][0] + m.m[1][3] * matrix.m[0][1] + m.m[2][3] * matrix.m[0][2]);
	matrix.m[1][3] = -(m.m[0][3] * matrix.m[1][0] + m.m[1][3] * matrix.m[1][1] + m.m[2][3] * matrix.m[1][2]);
	matrix.m[2][3] = -(m.m[0][3] * matrix.m[2][0] + m.m[1][3] * matrix.m[2][1] + m.m[2][3] * matrix.m[2][2]);
	
    matrix.m[3][3] = 1.0f;
	return matrix;
}