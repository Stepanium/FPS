#define OLC_PGE_APPLICATION

#include "olcPixelGameEngine.h"

#include <algorithm>
#include "MyMath.h"
// g++ -o main main.cpp -lX11 -lGL -lpthread -lpng -lstdc++fs -std=c++17

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
                int sign = 1;
                float color = brightness * brightness * 0.5f;
                brightness = sqrtf((brightness + 1.0f) * 0.5f);
                (sign) ? triViewed.color = olc::Pixel(200 * brightness + 55 * color, 200 * brightness, 200 * brightness) :
                triViewed.color = olc::Pixel(200 * brightness, 200 * brightness, 200 * brightness + 55 * color);

                triViewed.p[0] = matView * triTransformed.p[0];
                triViewed.p[1] = matView * triTransformed.p[1];
                triViewed.p[2] = matView * triTransformed.p[2];

                int nClippedTriangles = 0;
                triangle clipped[2];
                nClippedTriangles = Triangle_ClipAgainstPlane({0.0f, 0.0f, 1.0f}, {0.0f, 0.0f, 1.0f}, triViewed, clipped[0], clipped[1]);

                for (int n = 0; n < nClippedTriangles; n++)
                {
                    // Project from 3D to 2D space
                    triProjected.p[0] = matProj * clipped[n].p[0];
                    triProjected.p[1] = matProj * clipped[n].p[1];
                    triProjected.p[2] = matProj * clipped[n].p[2];
                    if (triProjected.p[0].w != 0.0f)
                        triProjected.p[0] /= triProjected.p[0].w;
                    if (triProjected.p[1].w != 0.0f)
                        triProjected.p[1] /= triProjected.p[1].w;
                    if (triProjected.p[2].w != 0.0f)
                        triProjected.p[2] /= triProjected.p[2].w;
                    triProjected.color = clipped[n].color;

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
            }

            // Sorting triangles back to frontby distance from camera

            std::sort(vTrianglesToRaster.begin(), vTrianglesToRaster.end(), [](triangle &A, triangle &B)
                      {
                float z1 = (A.p[0].z + A.p[1].z + A.p[2].z) * 0.333333;
                float z2 = (B.p[0].z + B.p[1].z + B.p[2].z) * 0.333333;
                return z1 < z2; });

            // Drawing sorted triangles
            for (auto &triToRaster : vTrianglesToRaster)
            {
                triangle clipped[2];
                std::list<triangle> listTriangles;
                listTriangles.push_back(triToRaster);
                int nNewTriangles = 1;

                for (int p = 0; p < 4; p++)
                {
                    int nTrisToAdd = 0;
                    while (nNewTriangles > 0)
                    {
                        triangle test = listTriangles.front();
                        listTriangles.pop_front();
                        nNewTriangles--;

                        switch (p)
                        {
                        case 0:
                            nTrisToAdd = Triangle_ClipAgainstPlane({0.0f, 0.0f, 0.0f}, {0.0f, 1.0f, 0.0f}, test, clipped[0], clipped[1]);
                            break;
                        case 1:
                            nTrisToAdd = Triangle_ClipAgainstPlane({0.0f, (float)ScreenHeight() - 1.0f, 0.0f}, {0.0f, -1.0f, 0.0f}, test, clipped[0], clipped[1]);
                            break;
                        case 2:
                            nTrisToAdd = Triangle_ClipAgainstPlane({0.0f, 0.0f, 0.0f}, {1.0f, 0.0f, 0.0f}, test, clipped[0], clipped[1]);
                            break;
                        case 3:
                            nTrisToAdd = Triangle_ClipAgainstPlane({(float)ScreenWidth() - 1.0f, 0.0f, 0.0f}, {-1.0f, 0.0f, 0.0f}, test, clipped[0], clipped[1]);
                            break;
                        }

                        for (int w = 0; w < nTrisToAdd; w++)
                            listTriangles.push_back(clipped[w]);
                    }
                    nNewTriangles = listTriangles.size();
                }

                // Draw the transformed, viewed, clipped, projected, sorted, clipped triangles
                for (auto &t : listTriangles)
                {
                    FillTriangle(t.p[0].x, t.p[0].y, t.p[1].x, t.p[1].y, t.p[2].x, t.p[2].y, t.color);
                    DrawTriangle(t.p[0].x, t.p[0].y, t.p[1].x, t.p[1].y, t.p[2].x, t.p[2].y, olc::BLACK);
                }
            }
        }

        return true;
    }
};

int main()
{
    FPS game;
    if (game.Construct(640, 480, 1, 1))
        game.Start();

    return 0;
}
