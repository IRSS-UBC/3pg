#include <gtest/gtest.h>
#include <../apps/DataInput.hpp>

// Demonstrate some basic assertions.
TEST(DataInputTests, tryInputParamGoodScalarInputs) {
	DataInput *dataInput = new DataInput();
	EXPECT_TRUE(dataInput->tryAddInputParam("Ratio NPP/GPP", { "0.47" }));
}