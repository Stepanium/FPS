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

struct triangle
{
    vec3d p[3];

    olc::Pixel color;

    vec3d getNormal()
    {
        vec3d line1 = p[1] - p[0];
        vec3d line2 = p[2] - p[0];
        vec3d normal;
        normal = line1.CrossProduct(line2);
        normal *= FastInverseSquareRoot(normal.GetLengthSqared());
        return normal;
    }
};

struct mesh
{
    std::vector<triangle> tris;

    bool LoadFromFile(std::string filename)
    {
        std::ifstream f(filename);
        if (!f.is_open())
            return false;

        // Local cache of vertecies
        std::vector<vec3d> verts;

        while (!f.eof())
        {
            char line[128];
            f.getline(line, 128);

            std::strstream s;
            s << line;

            char junk;
            if (line[0] == 'v')
            {
                vec3d v;
                s >> junk >> v.x >> v.y >> v.z;
                verts.push_back(v);
            }
            else if (line[0] == 'f')
            {
                triangle t;
                int f[3];
                s >> junk >> f[0] >> f[1] >> f[2];
                tris.push_back({verts[f[0] - 1], verts[f[1] - 1], verts[f[2] - 1]});
            }
        }
        return true;
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
    mat4x4 matrix;
	matrix.m[0][0] = m.m[0][0]; matrix.m[1][0] = m.m[0][1]; matrix.m[2][0] = m.m[0][2]; matrix.m[3][0] = 0.0f;
	matrix.m[0][1] = m.m[1][0]; matrix.m[1][1] = m.m[1][1]; matrix.m[2][1] = m.m[1][2]; matrix.m[3][1] = 0.0f;
	matrix.m[0][2] = m.m[2][0]; matrix.m[1][2] = m.m[2][1]; matrix.m[2][2] = m.m[2][2]; matrix.m[3][2] = 0.0f;
	
    //m03 = -dotProduct(TranslateVector, NewRight)
    //m13 = -dotProduct(TranslateVector, newUp)
    //m23 = -dotProduct(TranslateVector, newForward)
    matrix.m[0][3] = -(m.m[0][3] * matrix.m[0][0] + m.m[1][3] * matrix.m[0][1] + m.m[2][3] * matrix.m[0][2]);
	matrix.m[1][3] = -(m.m[0][3] * matrix.m[1][0] + m.m[1][3] * matrix.m[1][1] + m.m[2][3] * matrix.m[1][2]);
	matrix.m[2][3] = -(m.m[0][3] * matrix.m[2][0] + m.m[1][3] * matrix.m[2][1] + m.m[2][3] * matrix.m[2][2]);
	
    matrix.m[3][3] = 1.0f;
	return matrix;
}

	vec3d Vector_IntersectPlane(vec3d &plane_p, vec3d &plane_n, vec3d &lineStart, vec3d &lineEnd)
	{
		plane_n.Normalize();
		float plane_d = -plane_n.DotProduct(plane_p);
		float ad = lineStart.DotProduct(plane_n);
		float bd = lineEnd.DotProduct(plane_n);
		float t = (-plane_d - ad) / (bd - ad);
		vec3d lineStartToEnd = lineEnd - lineStart;
		vec3d lineToIntersect = lineStartToEnd * t;
		return lineStart + lineToIntersect;
	}

/*int Triangle_ClipAgainstPlane(vec3d plane_n, vec3d plane_p, triangle &in_tri, triangle &out_tri1, triangle &out_tri2)
{
    //make sure plain normal is normalized
    plane_n.Normalize();

    // return signed distance from point to plane
    // if dist is positive p is "inside" of the plane
    auto dist = [&](vec3d &p)
    {
        return plane_n.DotProduct(p) - plane_n.DotProduct(plane_p);
    };

    // temp arrays of inside and outside points
    vec3d* inside_points[3]; int nInsidePointCount = 0;
    vec3d* outside_points[3]; int nOutsidePointCount = 0;

    float d0 = dist(in_tri.p[0]);
    float d1 = dist(in_tri.p[1]);
    float d2 = dist(in_tri.p[2]);

    if (d0 >= 0.0f) { inside_points[nInsidePointCount++] = &in_tri.p[0]; }
    else { outside_points[nOutsidePointCount++] = &in_tri.p[0]; }

    if (d1 >= 0.0f) { inside_points[nInsidePointCount++] = &in_tri.p[1]; }
    else { outside_points[nOutsidePointCount++] = &in_tri.p[1]; }

    if (d2 >= 0.0f) { inside_points[nInsidePointCount++] = &in_tri.p[2]; }
    else { outside_points[nOutsidePointCount++] = &in_tri.p[2]; }

    if (nInsidePointCount = 0)
    {
        // all point of the triangle are outside of the plane
        // don't render the triangle
        return 0; // no returned triangles are valid
    }

    if (nInsidePointCount = 3)
    {
        // the whole of triangle is inside
        out_tri1 = in_tri;
        return 1; // one returned triangle is valid
    }

    if (nInsidePointCount = 1)
    {
        // one new triangle inside
        out_tri1.color = in_tri.color;

        out_tri1.p[0] = *inside_points[0];
        out_tri1.p[1] = Vector_IntersectsPlane(plane_n, plane_p, *outside_points[0], *inside_points[0]);
        out_tri1.p[2] = Vector_IntersectsPlane(plane_n, plane_p, *outside_points[1], *inside_points[0]);

        return 1; // one returned triangle is valid
    }

    if (nInsidePointCount = 2)
    {
        // two new triangles inside
        out_tri1.color = in_tri.color;

        out_tri1.p[0] = *inside_points[0];
        out_tri1.p[1] = *inside_points[1];
        out_tri1.p[2] = Vector_IntersectsPlane(plane_n, plane_p, *outside_points[0], *inside_points[0]);

        out_tri2.color = in_tri.color;

        out_tri2.p[0] = *inside_points[1];
        out_tri2.p[1] = out_tri1.p[2];
        out_tri2.p[2] = Vector_IntersectsPlane(plane_n, plane_p, *outside_points[0], *inside_points[1]);

        return 2; // two of the returned triangles are valid
    }

}*/

int Triangle_ClipAgainstPlane(vec3d plane_p, vec3d plane_n, triangle &in_tri, triangle &out_tri1, triangle &out_tri2)
	{
		// Make sure plane normal is indeed normal
		plane_n.Normalize();

		// Return signed shortest distance from point to plane, plane normal must be normalised
		auto dist = [&](vec3d &p)
		{
			// vec3d n = p / p.GetLength();
			return (plane_n.DotProduct(p) - plane_n.DotProduct(plane_p));
		};

		// Create two temporary storage arrays to classify points either side of plane
		// If distance sign is positive, point lies on "inside" of plane
		vec3d* inside_points[3];  int nInsidePointCount = 0;
		vec3d* outside_points[3]; int nOutsidePointCount = 0;

		// Get signed distance of each point in triangle to plane
		float d0 = dist(in_tri.p[0]);
		float d1 = dist(in_tri.p[1]);
		float d2 = dist(in_tri.p[2]);

		if (d0 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[0]; }
		else { outside_points[nOutsidePointCount++] = &in_tri.p[0]; }
		if (d1 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[1]; }
		else { outside_points[nOutsidePointCount++] = &in_tri.p[1]; }
		if (d2 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[2]; }
		else { outside_points[nOutsidePointCount++] = &in_tri.p[2]; }

		// Now classify triangle points, and break the input triangle into 
		// smaller output triangles if required. There are four possible
		// outcomes...

		if (nInsidePointCount == 0)
		{
			// All points lie on the outside of plane, so clip whole triangle
			// It ceases to exist

			return 0; // No returned triangles are valid
		}

		if (nInsidePointCount == 3)
		{
			// All points lie on the inside of plane, so do nothing
			// and allow the triangle to simply pass through
			out_tri1 = in_tri;

			return 1; // Just the one returned original triangle is valid
		}

		if (nInsidePointCount == 1 && nOutsidePointCount == 2)
		{
			// Triangle should be clipped. As two points lie outside
			// the plane, the triangle simply becomes a smaller triangle

			// Copy appearance info to new triangle
			out_tri1.color = olc::BLUE; // in_tri.color;

			// The inside point is valid, so keep that...
			out_tri1.p[0] = *inside_points[0];

			// but the two new points are at the locations where the 
			// original sides of the triangle (lines) intersect with the plane
			out_tri1.p[1] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0]);
			out_tri1.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[1]);

			return 1; // Return the newly formed single triangle
		}

		if (nInsidePointCount == 2 && nOutsidePointCount == 1)
		{
			// Triangle should be clipped. As two points lie inside the plane,
			// the clipped triangle becomes a "quad". Fortunately, we can
			// represent a quad with two new triangles

			// Copy appearance info to new triangles
			out_tri1.color = olc::GREEN; // in_tri.color;

			out_tri2.color = olc::RED; // in_tri.color;

			// The first triangle consists of the two inside points and a new
			// point determined by the location where one side of the triangle
			// intersects with the plane
			out_tri1.p[0] = *inside_points[0];
			out_tri1.p[1] = *inside_points[1];
			out_tri1.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0]);

			// The second triangle is composed of one of he inside points, a
			// new point determined by the intersection of the other side of the 
			// triangle and the plane, and the newly created point above
			out_tri2.p[0] = *inside_points[1];
			out_tri2.p[1] = out_tri1.p[2];
			out_tri2.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[1], *outside_points[0]);

			return 2; // Return two newly formed triangles which form a quad
		}
	}
