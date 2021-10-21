#pragma once
#include "glm/glm.hpp"
#include "Ray.h"
struct IIntersectable {
	virtual float intersect(const glm::vec3& ray)const = 0;
};

class Sphere : public IIntersectable {
private:
	glm::vec3 center;
	float radius;
public:
	Sphere(glm::vec3 _center, float _radius) : center(_center), radius(_radius) {}
	virtual float intersect(const glm::vec3& ray)const override
	{
		return length(ray - center) - radius;
	}
	glm::vec3 getCenter()const {
		return center;
	}
	void setCenter(const glm::vec3& _v) {
		this->center = _v;
	}

	float getRadius()const {
		return radius;
	}
};

class ControlPoint : public IIntersectable {
private:
	std::shared_ptr<glm::vec3> center;
	float radius;
public:
	ControlPoint(glm::vec3 _center, float _radius) :radius(_radius) {
		center = std::make_shared<glm::vec3>(_center);
	}

	virtual float intersect(const glm::vec3& ray)const override
	{
		return length(ray - *center.get()) - radius;
	}
	glm::vec3 getCenter()const {
		return *center.get();
	}
	void setCenter(const glm::vec3& _v) {
		*center.get() = _v;
	}

	float getRadius()const {
		return radius;
	}
};


class Plane : public IIntersectable {
private:
	glm::vec3 n;

public:
	Plane(glm::vec3 _n) : n(_n) {}
	virtual float intersect(const glm::vec3& ray)const override
	{
		glm::vec3 n2 = glm::normalize(n);
		return dot(n2, ray - glm::vec3(0, -2.0, 0));
	}
};

class Torus : public IIntersectable {
private:
	glm::vec3 c;
	glm::vec2 t;
public:
	Torus(glm::vec3 _c, glm::vec2 _t) : c(_c), t(_t) {}
	virtual float intersect(const glm::vec3& ray)const override
	{
		glm::vec2 q = glm::vec2(length(glm::vec2(c.x, c.z) - glm::vec2(ray.x, ray.z)) - t.x, c.y - ray.y);
		return length(q) - t.y;
	}
};

class Octahedron : public IIntersectable {
private:
	glm::vec3 c;
	float s;
public:
	Octahedron(glm::vec3 _c, float _s) : c(_c), s(_s) {}
	virtual float intersect(const glm::vec3& ray)const override
	{
		glm::vec3 p = abs(ray - c);
		return (p.x + p.y + p.z - s) * 0.57735027f;
	}
};

class BoundingBoxForBezierTransforms : public IIntersectable{
	std::shared_ptr<Sphere> obj;
	glm::vec3 min, max;

public:
	BoundingBoxForBezierTransforms(std::shared_ptr<Sphere> _obj): obj(_obj) {
		max = glm::vec3(1.0f, 2.0f, 1.0f);// obj->getCenter() + glm::vec3(obj->getRadius());
		min = obj->getCenter() - glm::vec3(obj->getRadius());

	}

	bool intersectBox(const glm::vec3& ray) const {
		if (ray.x > min.x && ray.x < max.x && ray.y > min.y && ray.y < max.y && ray.z > min.z && ray.z < max.z)
			return true;
		else
			return false;
	}

	virtual float intersect(const glm::vec3& ray)const override
	{
		//todo lin átalakítás min max alapján
		return obj->intersect(glm::vec3(ray.x, ray.y/2.0f, ray.z));
	}
};

