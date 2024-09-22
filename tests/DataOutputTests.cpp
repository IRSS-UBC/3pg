#include <gtest/gtest.h>

#include <filesystem>

#include <../apps/GDALRasterImage.hpp>
#include <../apps/DataOutput.hpp>

std::string refGridDir = "test_files/DataOutputTests/refgrid.tif";
std::string outDir = "test_files/DataOutputTests/outputs/";

GDALRasterImage img(refGridDir);
RefGridProperties refGrid = img.getRefGrid();
std::unordered_map<std::string, PPPG_OP_VAR> emptyVars;
std::function<void(std::string)> fake_log = [](std::string message) {};


TEST(DataOutputTests, filenames) {
	//get output path
	std::filesystem::path outPath = outDir;

	//delete everything in the output folder
	//https://stackoverflow.com/questions/59077670/c-delete-all-files-and-subfolders-but-keep-the-directory-itself
	for (const auto& entry : std::filesystem::directory_iterator(outPath)) {
		if (entry.path().string().find(".gitignore") == std::string::npos) {
			std::filesystem::remove_all(entry.path());
		}
	}

	//get month, year, name, and filepath
	int month1 = 1;
	int year1 = 2020;
	std::string name1 = "w";
	std::filesystem::path testOutput1 = outDir;
	testOutput1 /= name1 + std::to_string(year1) + std::to_string(month1) + ".tif";

	//get month, year, name, and filepath
	int month2 = 12;
	int year2 = 1;
	std::string name2 = "filename_for_output_parameter";
	std::filesystem::path testOutput2 = outDir;
	testOutput2 /= name2 + std::to_string(year2) + std::to_string(month2) + ".tif";

	//get month, year, name, and filepath
	int month3 = 6;
	int year3 = 10000;
	std::string name3 = "!@#$";
	std::filesystem::path testOutput3 = outDir;
	testOutput3 /= name3 + std::to_string(year3) + std::to_string(month3) + ".tif";

	//create dataOutput, set and write values
	DataOutput dataOutput(refGrid, outDir, emptyVars, fake_log);
	dataOutput.setVal(year1, month1, name1, 0, 0);
	dataOutput.setVal(year2, month2, name2, 0, 0);
	dataOutput.setVal(year3, month3, name3, 0, 0);
	dataOutput.writeRow(0);
	dataOutput.writeRow(1);

	//check the filenames exist that should
	EXPECT_TRUE(std::filesystem::exists(testOutput1));
	EXPECT_TRUE(std::filesystem::exists(testOutput2));
	EXPECT_TRUE(std::filesystem::exists(testOutput3));
}

TEST(DataOutputTests, nans) {
	//get output path
	std::filesystem::path outPath = outDir;

	//delete everything in the output folder
	//https://stackoverflow.com/questions/59077670/c-delete-all-files-and-subfolders-but-keep-the-directory-itself
	for (const auto& entry : std::filesystem::directory_iterator(outPath)) {
		if (entry.path().string().find(".gitignore") == std::string::npos) {
			std::filesystem::remove_all(entry.path());
		}
	}

	//get month, year, name, and filepath
	int nan1Month = 1;
	int nan1Year = 1;
	std::string nan1Name = "nanOne";
	std::filesystem::path nan1path = outDir;
	nan1path /= nan1Name + std::to_string(nan1Year) + std::to_string(nan1Month) + ".tif";

	//get month, year, name, and filepath
	int nan2Month = 1;
	int nan2Year = 1;
	std::string nan2Name = "nanTwo";
	std::filesystem::path nan2path = outDir;
	nan2path /= nan2Name + std::to_string(nan2Year) + std::to_string(nan2Month) + ".tif";

	std::unique_ptr<DataOutput> dataOutput = std::make_unique<DataOutput>(refGrid, outDir, emptyVars, fake_log);
	//only set nan1 values in vertical pattern:
	dataOutput->setVal(nan1Year, nan1Month, nan1Name, 0, 1); //row 0 column 0 = 1
	dataOutput->setVal(nan1Year, nan1Month, nan1Name, 2, 1); //row 1 column 0 = 1
														     //unset column 1 should be set to nan automatically

	//only set nan2 values in horizontal pattern:
	dataOutput->setVal(nan2Year, nan2Month, nan2Name, 2, 1); //row 1 column 0 = 1
	dataOutput->setVal(nan2Year, nan2Month, nan2Name, 3, 1); //row 1 column 1 = 0
															 //unset row 0 should be set to nan automatically
	//write
	dataOutput->writeRow(0);
	dataOutput->writeRow(1);

	//delete
	dataOutput.reset();

	//ensure files exist
	ASSERT_TRUE(std::filesystem::exists(nan1path));
	ASSERT_TRUE(std::filesystem::exists(nan2path));

	GDALRasterImage nan1Tif(nan1path.string());
	GDALRasterImage nan2Tif(nan2path.string());

	//note in the GetVal function, x is the first parameter (which is the column selector)
	//check nan tif 1
	EXPECT_EQ(nan1Tif.GetVal(0, 0), 1);				//[0,0] = 1
	EXPECT_TRUE(std::isnan(nan1Tif.GetVal(1, 0)));		//[0,1] = nan
	EXPECT_EQ(nan1Tif.GetVal(0, 1), 1);				//[1,0] = 1
	EXPECT_TRUE(std::isnan(nan1Tif.GetVal(1, 1)));		//[1,1] = nan

	//check nan tif 2
	EXPECT_TRUE(std::isnan(nan2Tif.GetVal(0, 0)));		//[0,0] = nan
	EXPECT_TRUE(std::isnan(nan2Tif.GetVal(1, 0)));		//[0,1] = nan
	EXPECT_EQ(nan2Tif.GetVal(0, 1), 1);				//[1,0] = 1
	EXPECT_EQ(nan2Tif.GetVal(1, 1), 1);				//[1,1] = 1
}

TEST(DataOutputTests, extremes) {
	//get output path
	std::filesystem::path outPath = outDir;

	//delete everything in the output folder
	//https://stackoverflow.com/questions/59077670/c-delete-all-files-and-subfolders-but-keep-the-directory-itself
	for (const auto& entry : std::filesystem::directory_iterator(outPath)) {
		if (entry.path().string().find(".gitignore") == std::string::npos) {
			std::filesystem::remove_all(entry.path());
		}
	}

	//get month, year, name, and filepath
	int month = 1;
	int year = 1;
	std::string name = "extremes";
	std::filesystem::path path = outDir;
	path /= name + std::to_string(year) + std::to_string(month) + ".tif";

	float smallPos = static_cast<float>(.000001);
	float smallNeg = static_cast<float>(-.000001);

	//create dataOutput, set and write values
	std::unique_ptr<DataOutput> dataOutput = std::make_unique<DataOutput>(refGrid, outDir, emptyVars, fake_log);
	dataOutput->setVal(year, month, name, 0, std::numeric_limits<float>::max());
	dataOutput->setVal(year, month, name, 1, std::numeric_limits<float>::min());
	dataOutput->setVal(year, month, name, 2, smallPos);
	dataOutput->setVal(year, month, name, 3, smallNeg);
	dataOutput->writeRow(0);
	dataOutput->writeRow(1);

	//delete
	dataOutput.reset();

	//ensure file exists
	ASSERT_TRUE(std::filesystem::exists(path));

	GDALRasterImage tif(path.string());

	//check values
	EXPECT_EQ(tif.GetVal(0, 0), std::numeric_limits<float>::max());
	EXPECT_EQ(tif.GetVal(1, 0), std::numeric_limits<float>::min());
	EXPECT_EQ(tif.GetVal(0, 1), smallPos);
	EXPECT_EQ(tif.GetVal(1, 1), smallNeg);
}
