#pragma once
#include "glm/glm.hpp"
struct Ray {
	glm::vec3 start, dir;
	Ray(glm::vec3 _start, glm::vec3 _dir) { start = _start; dir = normalize(_dir); }
	glm::vec3 operator-(const glm::vec3& rhs) const { return glm::vec3(this->start - rhs); }
};
