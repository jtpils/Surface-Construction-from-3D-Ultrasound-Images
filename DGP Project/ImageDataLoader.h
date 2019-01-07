#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <set>
#include <string>

#include <cstdlib>
#include <vector>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

class ImageDataLoader{

private:
	size_t numCols;
	size_t numRows;
	
	

public:
	size_t numImages;
	Eigen::MatrixXd readImage(std::string path);
	std::vector<Eigen::MatrixXd> readAllImages(std::string path, size_t numImages);

	std::vector<std::string> imagePath;
};
