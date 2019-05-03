
#include "../src/VistaAtmosphere/VistaAtmosphere.hpp"
#include <gtest/gtest.h>
#include <GL/gl.h>

// ==========================================================================
// tests - dummy

TEST(TestDummy, DummyFunction)
{
	glutInit();
	VistaAtmosphere oAtmosphere(VistaAtmosphere::Preset::EARTH);
	
	oAtmosphere.SetSunDirection(VistaVector3D(1,2,3));
	VistaVector3D ret = oAtmosphere.GetSunDirection();
   
    EXPECT_EQ(ret[0], 1);
	EXPECT_EQ(ret[1], 2);
	EXPECT_EQ(ret[2], 3);
}
