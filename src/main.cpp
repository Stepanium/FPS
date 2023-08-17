#define OLC_PGE_APPLICATION

#include "olcPixelGameEngine.h"
#include <vector>

//g++ -o main main.cpp -lX11 -lGL -lpthread -lpng -lstdc++fs -std=c++17

struct vec3d {
    float x = 0.0f, y = 0.0f, z = 0.0f;

    /*const vec3d operator = (const vec3d rhs) {
        return { x, y, z };
    }*/

    const vec3d operator + (const vec3d rhs) {
        vec3d out;
        out.x = x + rhs.x;
        out.y = y + rhs.y;
        out.z = z + rhs.z;
        return out;
    }

    const vec3d operator - (const vec3d rhs) {
        vec3d out;
        out.x = x - rhs.x;
        out.y = y - rhs.y;
        out.z = z - rhs.z;
        return out;
    }

    void operator += (const vec3d rhs) {
        this->x += rhs.x;
        this->y += rhs.y;
        this->z += rhs.z;
    }

    void operator -= (const vec3d rhs) {
        this->x -= rhs.x;
        this->y -= rhs.y;
        this->z -= rhs.z;
    }

    void operator /= (const float rhs) {
        this->x /= rhs;
        this->y /= rhs;
        this->z /= rhs;
    }

    void operator *= (const float rhs) {
        this->x *= rhs;
        this->y *= rhs;
        this->z *= rhs;
    }

    float getLength() {
        return sqrtf(x * x + y * y + z * z);
    }


    void Normalize() {
        float reverseL = 1.0f / getLength();
        this->x *= reverseL;
        this->y *= reverseL;
        this->z *= reverseL;
    }


};

struct triangle {
    vec3d points[3];

    vec3d normal() {
        vec3d line1 = points[1] - points[0];
        vec3d line2 = points[2] - points[0];
        vec3d normal;
        normal.x = line1.y * line2.z - line1.z * line2.y;
        normal.y = line1.z * line2.x - line1.x * line2.z;
        normal.z = line1.x * line2.y - line1.y * line2.x;
        float reverseL = 1.0f / normal.getLength();
        normal *= reverseL;
        return normal;
    }
};

struct mesh {
    std::vector<triangle> tris;
};

struct mat3x3 {
    float elem[3][3] = { 0.0f };

    vec3d operator * (const vec3d rhs) {
        vec3d output;
        output.x = elem[0][0] * rhs.x + elem[0][1] * rhs.y + elem[0][2] * rhs.z;
        output.y = elem[1][0] * rhs.x + elem[1][1] * rhs.y + elem[1][2] * rhs.z;
        output.z = elem[2][0] * rhs.x + elem[2][1] * rhs.y + elem[2][2] * rhs.z;
        return output;
    }
};

void MakeRotationMatrixAroundX(mat3x3& input, float angle) {
    input.elem[0][0] = 1.0f;
    input.elem[1][1] = cosf(angle);
    input.elem[2][1] = sinf(angle);
    input.elem[2][2] = cosf(angle);
    input.elem[1][2] = -sinf(angle);
}

void MakeRotationMatrixAroundY(mat3x3& input, float angle) {
    input.elem[1][1] = 1.0f;
    input.elem[0][0] = cosf(angle);
    input.elem[0][2] = sinf(angle);
    input.elem[2][2] = cosf(angle);
    input.elem[2][0] = -sinf(angle);
}

void MakeRotationMatrixAroundZ(mat3x3& input, float angle) {
    input.elem[2][2] = 1.0f;
    input.elem[0][0] = cosf(angle);
    input.elem[1][0] = sinf(angle);
    input.elem[1][1] = cosf(angle);
    input.elem[0][1] = -sinf(angle);
}

vec3d ScrProj(vec3d input, float a, float f, float zNear, float q) {
    float x = input.x;
    float y = input.y;
    float z = input.z;
    vec3d output = { a * f * x, f * y, q * (z - zNear) };

    if (z != 0.0f) {
        float Zinverse = 1.0f / z;
        output *= Zinverse;
    }
    return output;
}


class FPS : public olc::PixelGameEngine {
private:
    mesh meshCube;
    float aspectRatio;
    float FoV;
    float f;
    float zNear;
    float zFar;
    float q;

    vec3d camera;

	float fTheta;

public:    
    FPS() {
        sAppName = "Simple and fun FPS";
    }

private:    
    bool OnUserCreate() override {
        meshCube.tris = {

		// SOUTH
		{ 0.0f, 0.0f, 0.0f,    0.0f, 1.0f, 0.0f,    1.0f, 1.0f, 0.0f },
		{ 0.0f, 0.0f, 0.0f,    1.0f, 1.0f, 0.0f,    1.0f, 0.0f, 0.0f },

		// EAST                                                      
		{ 1.0f, 0.0f, 0.0f,    1.0f, 1.0f, 0.0f,    1.0f, 1.0f, 1.0f },
		{ 1.0f, 0.0f, 0.0f,    1.0f, 1.0f, 1.0f,    1.0f, 0.0f, 1.0f },

		// NORTH                                                     
		{ 1.0f, 0.0f, 1.0f,    1.0f, 1.0f, 1.0f,    0.0f, 1.0f, 1.0f },
		{ 1.0f, 0.0f, 1.0f,    0.0f, 1.0f, 1.0f,    0.0f, 0.0f, 1.0f },

		// WEST                                                      
		{ 0.0f, 0.0f, 1.0f,    0.0f, 1.0f, 1.0f,    0.0f, 1.0f, 0.0f },
		{ 0.0f, 0.0f, 1.0f,    0.0f, 1.0f, 0.0f,    0.0f, 0.0f, 0.0f },

		// TOP                                                       
		{ 0.0f, 1.0f, 0.0f,    0.0f, 1.0f, 1.0f,    1.0f, 1.0f, 1.0f },
		{ 0.0f, 1.0f, 0.0f,    1.0f, 1.0f, 1.0f,    1.0f, 1.0f, 0.0f },

		// BOTTOM                                                    
		{ 1.0f, 0.0f, 1.0f,    0.0f, 0.0f, 1.0f,    0.0f, 0.0f, 0.0f },
		{ 1.0f, 0.0f, 1.0f,    0.0f, 0.0f, 0.0f,    1.0f, 0.0f, 0.0f },

		};

        aspectRatio = (float)ScreenHeight() / (float)ScreenWidth();
        
        FoV = 90.0f;
        f = 1.0f / tanf(FoV * 0.5f / 180.f * M_PI);

        zNear = 0.1f;
        zFar = 1000.0f;
        q = zFar / (zFar - zNear);

        fTheta = 0.0f;

        return true;
    }

    bool OnUserUpdate(float fElapsedTime) override {
        fTheta += fElapsedTime;
        // if(fTheta >= M_PI) fTheta -= M_PI;
        
        Clear(olc::BLACK);

        //Drawing triangles

        for(auto tri : meshCube.tris) {
            triangle triRotatedX, triRotatedXZ, triTranslated, triProjected, triScaled;

            mat3x3 matRotX, matRotZ;
            MakeRotationMatrixAroundX(matRotX, fTheta);
            MakeRotationMatrixAroundZ(matRotZ, fTheta * 0.5f);

            triRotatedX = tri;
            triRotatedX.points[0] = matRotX * triRotatedX.points[0];
            triRotatedX.points[1] = matRotX * triRotatedX.points[1];
            triRotatedX.points[2] = matRotX * triRotatedX.points[2];

            triRotatedXZ.points[0] = matRotZ * triRotatedX.points[0];
            triRotatedXZ.points[1] = matRotZ * triRotatedX.points[1];
            triRotatedXZ.points[2] = matRotZ * triRotatedX.points[2];

            triTranslated = triRotatedXZ;

            triTranslated.points[0].z += 3.0f;
            triTranslated.points[1].z += 3.0f;
            triTranslated.points[2].z += 3.0f;

            vec3d normal = triTranslated.normal();
            vec3d ray = triTranslated.points[0] - camera;
            ray.Normalize();
            if((normal.x * ray.x + normal.y * ray.y + normal.z * ray.z) < 0.0f) {
                
                triProjected.points[0] = ScrProj(triTranslated.points[0], aspectRatio, f, zNear, q);
                triProjected.points[1] = ScrProj(triTranslated.points[1], aspectRatio, f, zNear, q);
                triProjected.points[2] = ScrProj(triTranslated.points[2], aspectRatio, f, zNear, q);

                triProjected.points[0].x += 1.0f; triProjected.points[0].y += 1.0f;
                triProjected.points[1].x += 1.0f; triProjected.points[1].y += 1.0f;
                triProjected.points[2].x += 1.0f; triProjected.points[2].y += 1.0f;

                triProjected.points[0].x *= 0.5f * (float)ScreenWidth();
                triProjected.points[0].y *= 0.5f * (float)ScreenHeight();

                triProjected.points[1].x *= 0.5f * (float)ScreenWidth();
                triProjected.points[1].y *= 0.5f * (float)ScreenHeight();

                triProjected.points[2].x *= 0.5f * (float)ScreenWidth();
                triProjected.points[2].y *= 0.5f * (float)ScreenHeight();

                DrawTriangle(
                    triProjected.points[0].x, triProjected.points[0].y, 
                    triProjected.points[1].x, triProjected.points[1].y, 
                    triProjected.points[2].x, triProjected.points[2].y, 
                    olc::WHITE
                );
            }
        }

        return true;
    }
};

int main() {
    FPS game;
    if(game.Construct(256, 256, 4, 4))
        game.Start();

    return 0;
}
