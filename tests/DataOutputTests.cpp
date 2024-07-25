#include <gtest/gtest.h>

#include <filesystem>

#include <../apps/GDALRasterImage.hpp>
#include <../apps/DataOutput.hpp>

std::string refgridDir = "test_files/DataOutputTests/refgrid.tif";
std::string outDir = "test_files/DataOutputTests/outputs/";

TEST(DataOutputTests, filenames) {
	//get refgrid
	GDALRasterImage* refgrid = new GDALRasterImage(refgridDir);

	//get output path
	std::filesystem::path outPath = outDir;

	//delete everything in the output folder
	//https://stackoverflow.com/questions/59077670/c-delete-all-files-and-subfolders-but-keep-the-directory-itself
	for (const auto& entry : std::filesystem::directory_iterator(outPath)) {
		std::filesystem::remove_all(entry.path());
	}

	int month1 = 1;
	int year1 = 2020;
	std::string name1 = "w";
	std::filesystem::path testOutput1 = outDir;
	testOutput1 /= name1 + std::to_string(year1) + std::to_string(month1) + ".tif";

	int month2 = 12;
	int year2 = 1;
	std::string name2 = "filename_for_output_parameter";
	std::filesystem::path testOutput2 = outDir;
	testOutput2 /= name2 + std::to_string(year2) + std::to_string(month2) + ".tif";

	int month3 = 6;
	int year3 = 10000;
	std::string name3 = "!@#$";
	std::filesystem::path testOutput3 = outDir;
	testOutput3 /= name3 + std::to_string(year3) + std::to_string(month3) + ".tif";

	DataOutput* dataOutput = new DataOutput(refgrid, outDir);
	dataOutput->setVal(year1, month1, name1, 0, 0);
	dataOutput->setVal(year2, month2, name2, 0, 0);
	dataOutput->setVal(year3, month3, name3, 0, 0);
	dataOutput->writeRow(0);
	dataOutput->writeRow(1);
	delete dataOutput;

	EXPECT_TRUE(std::filesystem::exists(testOutput1));
	EXPECT_TRUE(std::filesystem::exists(testOutput2));
	EXPECT_TRUE(std::filesystem::exists(testOutput3));
}