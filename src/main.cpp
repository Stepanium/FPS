#define OLC_PGE_APPLICATION

#include "olcPixelGameEngine.h"

#include <algorithm>
#include "MyMath.h"
// g++ -o main main.cpp -lX11 -lGL -lpthread -lpng -lstdc++fs -std=c++17

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

class FPS : public olc::PixelGameEngine
{
private:
    mesh meshCube;
    float aspectRatio;
    float FoV;
    float zNear;
    float zFar;

    vec3d vCamera;
    vec3d vLookDir;
    vec3d vUp;
    vec3d vRight;
    float fYaw;
    float fTheta;

    vec3d light_direction;

    mat4x4 matProj;

public:
    FPS()
    {
        sAppName = "Simple and fun FPS";
    }

private:
    bool OnUserCreate() override
    {
        meshCube.LoadFromFile("axis.obj");

        aspectRatio = (float)ScreenHeight() / (float)ScreenWidth();
        FoV = 90.0f;
        zNear = 0.1f;
        zFar = 1000.0f;
        // Making projection matrix
        matProj = Matrix_MakeProjection(aspectRatio, FoV, zNear, zFar);

        vCamera = {0.0f, 0.0f, 0.0f};
        fTheta = 0.0f;
        fYaw = 0.0f;
        vLookDir = {0.0f, 0.0f, 1.0f};
        vUp = {0.0f, 1.0f, 0.0f};
        vRight = {1.0f, 0.0f, 0.0f};

        return true;
    }

    bool OnUserUpdate(float fElapsedTime) override
    {

        // INPUT//
        // TRANSLATION//
        if (GetKey(olc::Key::W).bHeld)
            vCamera += vLookDir * (8.0f * fElapsedTime);
        if (GetKey(olc::Key::S).bHeld)
            vCamera -= vLookDir * (8.0f * fElapsedTime);
        if (GetKey(olc::Key::A).bHeld)
            vCamera += vRight * (8.0f * fElapsedTime);
        if (GetKey(olc::Key::D).bHeld)
            vCamera -= vRight * (8.0f * fElapsedTime);
        //////

        // ROTATION//
        if (GetKey(olc::Key::LEFT).bHeld)
            fTheta -= 1.0f * fElapsedTime;
        if (GetKey(olc::Key::RIGHT).bHeld)
            fTheta += 1.0f * fElapsedTime;
        if (GetKey(olc::Key::UP).bHeld)
            fYaw -= 1.0f * fElapsedTime;
        if (GetKey(olc::Key::DOWN).bHeld)
            fYaw += 1.0f * fElapsedTime;
        //////

        // Object transltion and rotation
        mat4x4 matRotX, matRotY, matRotXY, matTrans, matWorld;
        matRotX = Matrix_MakeRotationAroundX(0.0f);
        matRotY = Matrix_MakeRotationAroundY(0.0f);
        matRotXY = matRotY * matRotX;
        matTrans = Matrix_MakeTranslation({0.0f, 0.0f, 15.0f});
        matWorld = matTrans * matRotXY;
        //////

        // CAMERA//
        //  Camera translation and rotation
        matRotX = Matrix_MakeRotationAroundX(fYaw);
        matRotY = Matrix_MakeRotationAroundY(fTheta);
        matRotXY = matRotY * matRotX;

        // New forward direction
        vLookDir = {0.0f, 0.0f, 1.0f};
        vLookDir = matRotXY * vLookDir;

        // New upward direction
        vUp = {0.0f, 1.0f, 0.0f};
        vec3d a = vLookDir * vUp.DotProduct(vLookDir);
        vUp -= a;
        vUp.Normalize();

        // New rightside direction
        vRight = vUp.CrossProduct(vLookDir);
        //////

        // MATRICIES FOR TRIS//
        // making inverse matrix for camera movement
        mat4x4 matCamera = Matrix_MakePointAt(vCamera, vCamera + vLookDir, vUp);
        mat4x4 matView = Matrix_MakeInverseForRotTrans(matCamera);
        //////

        // LIGHTING DIR//
        light_direction = {0.0f, 0.0f, -1.0f};
        light_direction.Normalize();
        //////

        std::vector<triangle> vTrianglesToRaster;

        Clear(olc::BLACK);
        // Calculating triangles
        for (auto tri : meshCube.tris)
        {
            triangle triTransformed, triProjected, triViewed;

            // World Matrix Transform
            triTransformed.p[0] = matWorld * tri.p[0];
            triTransformed.p[1] = matWorld * tri.p[1];
            triTransformed.p[2] = matWorld * tri.p[2];

            vec3d normal = triTransformed.getNormal();
            vec3d vCameraRay = triTransformed.p[0] - vCamera;
            vCameraRay.Normalize();
            // If visible from camera's point of view
            if (normal.DotProduct(vCameraRay) < 0.0f)
            {
                // Lighting
                float brightness = normal.DotProduct(light_direction);
                if (brightness < 0.05f)
                    triViewed.color = olc::Pixel(20, 20, 40);
                else
                    triViewed.color = olc::Pixel(255 * brightness, 200 * brightness, 200 * brightness);

                triViewed.p[0] = matView * triTransformed.p[0];
                triViewed.p[1] = matView * triTransformed.p[1];
                triViewed.p[2] = matView * triTransformed.p[2];

                // Project from 3D to 2D space
                triProjected.p[0] = matProj * triViewed.p[0];
                triProjected.p[1] = matProj * triViewed.p[1];
                triProjected.p[2] = matProj * triViewed.p[2];
                if (triProjected.p[0].w != 0.0f)
                    triProjected.p[0] /= triProjected.p[0].w;
                if (triProjected.p[1].w != 0.0f)
                    triProjected.p[1] /= triProjected.p[1].w;
                if (triProjected.p[2].w != 0.0f)
                    triProjected.p[2] /= triProjected.p[2].w;
                triProjected.color = triViewed.color;

                // Translate to the center of the screen
                vec3d vOffsetView = {1.0f, 1.0f, 0.0f};
                triProjected.p[0] += vOffsetView;
                triProjected.p[1] += vOffsetView;
                triProjected.p[2] += vOffsetView;

                // Scale to screen width and height
                triProjected.p[0].x *= 0.5f * (float)ScreenWidth();
                triProjected.p[0].y *= 0.5f * (float)ScreenHeight();
                triProjected.p[1].x *= 0.5f * (float)ScreenWidth();
                triProjected.p[1].y *= 0.5f * (float)ScreenHeight();
                triProjected.p[2].x *= 0.5f * (float)ScreenWidth();
                triProjected.p[2].y *= 0.5f * (float)ScreenHeight();

                vTrianglesToRaster.push_back(triProjected);
            }

            // Sorting triangles back to frontby distance from camera

            std::sort(vTrianglesToRaster.begin(), vTrianglesToRaster.end(), [](triangle &A, triangle &B)
                      {
                float z1 = (A.p[0].z + A.p[1].z + A.p[2].z) * 0.333333;
                float z2 = (B.p[0].z + B.p[1].z + B.p[2].z) * 0.333333;
                return z1 < z2; });

            // Drawing sorted triangles
            for (auto &triProjected : vTrianglesToRaster)
            {
                FillTriangle(
                    triProjected.p[0].x, triProjected.p[0].y,
                    triProjected.p[1].x, triProjected.p[1].y,
                    triProjected.p[2].x, triProjected.p[2].y,
                    triProjected.color);

                DrawTriangle(
                    triProjected.p[0].x, triProjected.p[0].y,
                    triProjected.p[1].x, triProjected.p[1].y,
                    triProjected.p[2].x, triProjected.p[2].y,
                    olc::BLACK);
            }
        }

        return true;
    }
};

int main()
{
    FPS game;
    if (game.Construct(360, 360, 2, 2))
        game.Start();

    return 0;
}
