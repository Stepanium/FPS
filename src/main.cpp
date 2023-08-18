#define OLC_PGE_APPLICATION

#include "olcPixelGameEngine.h"
#include <vector>
#include <fstream>
#include <strstream>

// g++ -o main main.cpp -lX11 -lGL -lpthread -lpng -lstdc++fs -std=c++17

struct vec3d
{
    float x = 0.0f, y = 0.0f, z = 0.0f;

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

    float GetLength()
    {
        return sqrtf(x * x + y * y + z * z);
    }

    void Normalize()
    {
        float reverseL = 1.0f / GetLength();
        this->x *= reverseL;
        this->y *= reverseL;
        this->z *= reverseL;
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
    vec3d points[3];

    vec3d getNormal()
    {
        vec3d line1 = points[1] - points[0];
        vec3d line2 = points[2] - points[0];
        vec3d normal;
        normal = line1.CrossProduct(line2);
        normal.Normalize();
        return normal;
    }
};

struct mesh
{
    std::vector<triangle> tris;

    bool LoadFromFile(std::string filename)
    {
    }
};

struct mat4x4
{
    float m[4][4] = {0.0f};
    float w = 1.0f;

    vec3d operator*(const vec3d rhs)
    {
        vec3d output;
        output.x = m[0][0] * rhs.x + m[0][1] * rhs.y + m[0][2] * rhs.z + m[0][3];
        output.y = m[1][0] * rhs.x + m[1][1] * rhs.y + m[1][2] * rhs.z + m[1][3];
        output.z = m[2][0] * rhs.x + m[2][1] * rhs.y + m[2][2] * rhs.z + m[2][3];
        w = m[3][0] * rhs.x + m[3][1] * rhs.y + m[3][2] * rhs.z + m[3][3];
        return output;
    }
};

void MakeRotationMatrixAroundX(mat4x4 &input, float angle)
{
    input.m[0][0] = 1.0f;
    input.m[1][1] = cosf(angle);
    input.m[2][1] = sinf(angle);
    input.m[2][2] = cosf(angle);
    input.m[1][2] = -sinf(angle);
    input.m[3][3] = 1.0f;
}

void MakeRotationMatrixAroundY(mat4x4 &input, float angle)
{
    input.m[1][1] = 1.0f;
    input.m[0][0] = cosf(angle);
    input.m[0][2] = sinf(angle);
    input.m[2][2] = cosf(angle);
    input.m[2][0] = -sinf(angle);
    input.m[3][3] = 1.0f;
}

void MakeRotationMatrixAroundZ(mat4x4 &input, float angle)
{
    input.m[2][2] = 1.0f;
    input.m[0][0] = cosf(angle);
    input.m[1][0] = sinf(angle);
    input.m[1][1] = cosf(angle);
    input.m[0][1] = -sinf(angle);
    input.m[3][3] = 1.0f;
}

void MakeProjectionMatrix(mat4x4 &input, float aspectRatio, float f, float zNear, float zFar)
{
    float q = zFar / (zFar - zNear);
    input.m[0][0] = aspectRatio * f;
    input.m[1][1] = f;
    input.m[2][2] = q;
    input.m[2][3] = -zNear * q;
    input.m[3][2] = 1.0f;
    input.m[3][3] = 0.0f;
}

vec3d ProjectToScreen(vec3d input, mat4x4 matProj)
{
    vec3d output = matProj * input;

    if (matProj.w != 0.0f)
    {
        float wInverse = 1.0f / matProj.w;
        output *= wInverse;
    }
    return output;
}

class FPS : public olc::PixelGameEngine
{
private:
    mesh meshCube;
    float aspectRatio;
    float FoV;
    float f;
    float zNear;
    float zFar;
    float q;

    vec3d vCamera;
    vec3d vLight;

    float fTheta;

public:
    FPS()
    {
        sAppName = "Simple and fun FPS";
    }

private:
    bool OnUserCreate() override
    {
        meshCube.tris = {

            // SOUTH
            {0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 1.0f, 1.0f, 0.0f},
            {0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f},

            // EAST
            {1.0f, 0.0f, 0.0f, 1.0f, 1.0f, 0.0f, 1.0f, 1.0f, 1.0f},
            {1.0f, 0.0f, 0.0f, 1.0f, 1.0f, 1.0f, 1.0f, 0.0f, 1.0f},

            // NORTH
            {1.0f, 0.0f, 1.0f, 1.0f, 1.0f, 1.0f, 0.0f, 1.0f, 1.0f},
            {1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 1.0f, 0.0f, 0.0f, 1.0f},

            // WEST
            {0.0f, 0.0f, 1.0f, 0.0f, 1.0f, 1.0f, 0.0f, 1.0f, 0.0f},
            {0.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f},

            // TOP
            {0.0f, 1.0f, 0.0f, 0.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f},
            {0.0f, 1.0f, 0.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 0.0f},

            // BOTTOM
            {1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f},
            {1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f},

        };

        aspectRatio = (float)ScreenHeight() / (float)ScreenWidth();

        FoV = 90.0f;
        f = 1.0f / tanf(FoV * 0.5f * M_PI / 180.f);

        zNear = 0.1f;
        zFar = 1000.0f;
        q = zFar / (zFar - zNear);

        fTheta = 0.0f;
        vLight = {-1.0f, -1.0f, -1.0f};
        vLight.Normalize();

        return true;
    }

    bool OnUserUpdate(float fElapsedTime) override
    {
        fTheta += fElapsedTime;

        Clear(olc::BLACK);

        // Drawing triangles
        for (auto tri : meshCube.tris)
        {
            triangle triRotatedX, triRotatedXY, triRotatedXYZ, triTranslated, triProjected;

            mat4x4 matRotX, matRotY, matRotZ, matProj;
            MakeRotationMatrixAroundX(matRotX, fTheta);
            MakeRotationMatrixAroundY(matRotY, fTheta * 0.1f);
            MakeRotationMatrixAroundZ(matRotZ, fTheta * 0.3f);
            MakeProjectionMatrix(matProj, aspectRatio, f, q, zNear);

            triRotatedX.points[0] = matRotX * tri.points[0];
            triRotatedX.points[1] = matRotX * tri.points[1];
            triRotatedX.points[2] = matRotX * tri.points[2];

            triRotatedXY.points[0] = matRotY * triRotatedX.points[0];
            triRotatedXY.points[1] = matRotY * triRotatedX.points[1];
            triRotatedXY.points[2] = matRotY * triRotatedX.points[2];

            triRotatedXYZ.points[0] = matRotZ * triRotatedXY.points[0];
            triRotatedXYZ.points[1] = matRotZ * triRotatedXY.points[1];
            triRotatedXYZ.points[2] = matRotZ * triRotatedXY.points[2];

            triTranslated = triRotatedXYZ;

            triTranslated.points[0].z += 3.0f;
            triTranslated.points[1].z += 3.0f;
            triTranslated.points[2].z += 3.0f;

            vec3d normal = triTranslated.getNormal();
            vec3d ray = triTranslated.points[0] - vCamera;
            ray.Normalize();
            // If visible from camera's point of view
            if (normal.DotProduct(ray) < 0.0f)
            {
                // Project from 3D to 2D space
                triProjected.points[0] = ProjectToScreen(triTranslated.points[0], matProj);
                triProjected.points[1] = ProjectToScreen(triTranslated.points[1], matProj);
                triProjected.points[2] = ProjectToScreen(triTranslated.points[2], matProj);

                // Translate to the center of the screen
                triProjected.points[0].x += 1.0f;
                triProjected.points[0].y += 1.0f;
                triProjected.points[1].x += 1.0f;
                triProjected.points[1].y += 1.0f;
                triProjected.points[2].x += 1.0f;
                triProjected.points[2].y += 1.0f;

                // Scale to screen width and height
                triProjected.points[0].x *= 0.5f * (float)ScreenWidth();
                triProjected.points[0].y *= 0.5f * (float)ScreenHeight();
                triProjected.points[1].x *= 0.5f * (float)ScreenWidth();
                triProjected.points[1].y *= 0.5f * (float)ScreenHeight();
                triProjected.points[2].x *= 0.5f * (float)ScreenWidth();
                triProjected.points[2].y *= 0.5f * (float)ScreenHeight();

                float brightness = normal.DotProduct(vLight);
                if (brightness < 0.0f)
                    brightness = 0.0f;
                //brightness *= 0.9f;
                //brightness += 0.1;

                FillTriangle(
                    triProjected.points[0].x, triProjected.points[0].y,
                    triProjected.points[1].x, triProjected.points[1].y,
                    triProjected.points[2].x, triProjected.points[2].y,
                    olc::Pixel(255 * brightness, 200 * brightness, 200 * brightness));

                DrawTriangle(
                    triProjected.points[0].x, triProjected.points[0].y,
                    triProjected.points[1].x, triProjected.points[1].y,
                    triProjected.points[2].x, triProjected.points[2].y,
                    olc::BLACK);
            }
        }

        return true;
    }
};

int main()
{
    FPS game;
    if (game.Construct(256, 256, 4, 4))
        game.Start();

    return 0;
}
