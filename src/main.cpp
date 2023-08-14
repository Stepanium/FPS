#define OLC_PGE_APPLICATION

#include "olcPixelGameEngine.h"
#include <vector>

struct vec3d {
    float x, y, z;

    const vec3d& operator + (const vec3d& rhs) {
        vec3d out;
        out.x = x + rhs.x;
        out.y = y + rhs.y;
        out.z = z + rhs.z;
        return out;
    }

    const vec3d& operator - (const vec3d& rhs) {
        vec3d out;
        out.x = x - rhs.x;
        out.y = y - rhs.y;
        out.z = z - rhs.z;
        return out;
    }

    void operator += (const vec3d& rhs) {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
    }

    void operator -= (const vec3d& rhs) {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
    }

};

struct triangle {
    vec3d points[3];
};

struct mesh {
    std::vector<triangle> tris;
};

class FPS : public olc::PixelGameEngine {
    
public:    
    FPS() {
        sAppName = "Simple and fun FPS";
    }

private:    
    bool OnUserCreate() override {

        return true;
    }

    bool OnUserUpdate(float fElapsedTime) override {

        return true;
    }
};

int main() {
    FPS game;
    if(game.Construct(640, 480, 1, 1))
        game.Start();

    return 0;
}
