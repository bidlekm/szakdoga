#include "src/imgui/imgui.h"
#include "src/imgui/imgui_impl_glfw.h"
#include "src/imgui/imgui_impl_opengl3.h"
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstddef>
#include <string>
#include <vector>
#include <memory>
#include <omp.h>
#include "objects.h"
#include "Ray.h"
#include <algorithm>


#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif
#pragma warning(disable : 4996)

const int size = 500;
int iter = 0;
struct ShaderProgramSource
{
	std::string VertexSource;
	std::string FragmentSource;
};

static ShaderProgramSource ParseShader(const std::string& filepath)
{
	std::ifstream stream(filepath);
	std::string line;
	std::stringstream ss[2];
	int type = -1;
	while (getline(stream, line))
	{
		if (line.find("#shader") != std::string::npos)
		{
			if (line.find("vertex") != std::string::npos)
				type = 0;
			else if (line.find("fragment") != std::string::npos)
				type = 1;
		}
		else
		{
			ss[(int)type] << line << '\n';
		}

	}
	return { ss[0].str(), ss[1].str() };
}

static unsigned int CompileShader(unsigned int type, const std::string& source)
{
	unsigned int id = glCreateShader(type);
	const char* src = source.c_str();
	glShaderSource(id, 1, &src, nullptr);
	glCompileShader(id);

	int result;
	glGetShaderiv(id, GL_COMPILE_STATUS, &result);
	if (result == GL_FALSE) {
		int length;
		glGetShaderiv(id, GL_INFO_LOG_LENGTH, &length);
		char* message = (char*)alloca(length * sizeof(char));
		glGetShaderInfoLog(id, length, &length, message);
		std::cout << message << std::endl;
	}
	return id;
}

static int CreateProgram(const std::string& vertexShader, const std::string& fragmentShader)
{
	unsigned int program = glCreateProgram();
	unsigned int vs = CompileShader(GL_VERTEX_SHADER, vertexShader);
	unsigned int fs = CompileShader(GL_FRAGMENT_SHADER, fragmentShader);
	glAttachShader(program, vs);
	glAttachShader(program, fs);
	glLinkProgram(program);
	glValidateProgram(program);
	return program;
}

struct Camera {
	glm::vec3 eye, lookat, right, up;
	float fov;
	void set(glm::vec3 _eye, glm::vec3 _lookat, glm::vec3 vup, float _fov) {
		eye = _eye;
		lookat = _lookat;
		glm::vec3 w = eye - lookat;
		fov = _fov;
		float focus = length(w);
		right = normalize(cross(w, vup)) * focus * tanf(fov / 2);
		up = normalize(cross(right, w)) * focus * tanf(fov / 2);
	}
	Ray getRay(int X, int Y){
		glm::vec3 dir = lookat + right * (2.0f * (X + 0.5f) / size - 1) + up * (2.0f * (Y + 0.5f) / size - 1) - eye;
		return Ray(eye, dir);
	}
};

struct Scene {

	std::vector<std::shared_ptr<IIntersectable>> objects;
	std::shared_ptr<BoundingBoxForBezierTransforms> boundingBox;
	std::vector<std::shared_ptr<ControlPoint>> bezierControlPoints;
	unsigned int textureId = 0;
	std::vector<glm::vec4> image;
	Camera camera;
	size_t controlPointNumberX, controlPointNumberY, controlPointNumberZ;
	std::vector<std::vector<std::vector<glm::vec3>>> points;
	bool deformed, showControlPoints;
	Scene()
	{
		textureId = 0;
	}
	void setControlPoints(int x, int y, int z){
		controlPointNumberX = x;
		controlPointNumberY = y;
		controlPointNumberZ = z;
		points = std::vector<std::vector<std::vector<glm::vec3>>>(controlPointNumberX, std::vector<std::vector<glm::vec3>>(controlPointNumberY, std::vector<glm::vec3>(controlPointNumberZ, glm::vec3())));
		for (size_t i = 0; i < controlPointNumberX; i++)
		{
			for (size_t j = 0; j < controlPointNumberY; j++)
			{
				for (size_t k = 0; k < controlPointNumberZ; k++)
				{
					points[i][j][k] = glm::vec3(static_cast<float>(i) / static_cast<float>(controlPointNumberX - 1),
						static_cast<float>(j) / static_cast<float>(controlPointNumberY - 1),
						static_cast<float>(k) / static_cast<float>(controlPointNumberZ - 1));
					if (j != 0) {
						points[i][j][k].y *= 2.0f;
					}
				}
			}
		}

		bezierControlPoints.clear();
		
		for (size_t i = 0; i < controlPointNumberX; i++)
		{
			for (size_t j = 0; j < controlPointNumberY; j++)
			{
				for (size_t k = 0; k < controlPointNumberZ; k++)
				{
					bezierControlPoints.push_back(std::shared_ptr<ControlPoint>(new ControlPoint(points[i][j][k], 0.05f)));
				}
			}
		}

		points[0][1][0].x = 0.25;
		points[0][1][1].x = 0.25;
		points[0][1][2].x = 0.25;
		points[2][1][0].x = 0.75;
		points[2][1][1].x = 0.75;
		points[2][1][2].x = 0.75;
	}

	Scene(int width, int height, const std::vector<glm::vec4>& image, int sampling = GL_LINEAR) {
		deformed = false;
		showControlPoints = false;
		glm::vec3 eye = glm::vec3(0.0f, 1.0f, -8.0f);
		glm::vec3 up = glm::vec3(0, 1, 0), lookat = glm::vec3(0, 0, 0);
		float fov = 45.0f * M_PI / 180.0f;
		camera.set(eye, lookat, up, fov);
		textureId = 0;

		setControlPoints(3, 2, 3);
		
		objects.push_back(std::shared_ptr<IIntersectable>(new Sphere(glm::vec3(0.5f, 0.5f, 0.5f), 0.45f)));
		objects.push_back(std::shared_ptr<IIntersectable>(new Plane(glm::vec3(0.0f, 1.0f, 0.0f))));
		objects.push_back(std::shared_ptr<IIntersectable>(new Torus(glm::vec3(0.0f, 0.0f, 0.0f), glm::vec2(1.5f, 0.1f))));
		objects.push_back(std::shared_ptr<IIntersectable>(new Octahedron(glm::vec3(2.5f, 0.5f, 0.0f), 1.0f)));
		
		boundingBox = std::shared_ptr<BoundingBoxForBezierTransforms>(new BoundingBoxForBezierTransforms(std::dynamic_pointer_cast<Sphere>(objects[0])));

		create(width, height, image, sampling);
		build();
	}
	void build() {
#pragma omp parallel for
		for (int Y = 0; Y < size; Y++) {
#pragma omp parallel for num_threads(10)
			for (int X = 0; X < size; X++) {
				glm::vec3 color = rayTrace(camera.getRay(X, Y));
				image[Y * size + X] = glm::vec4(color.x, color.y, color.z, 1.0f);
			}
		}
	}
	void render() {
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, size, size, 0, GL_RGBA, GL_FLOAT, &image[0]); // To GPU
	}

	void create(int width, int height, const std::vector<glm::vec4>& image, int sampling = GL_LINEAR) {
		this->image = image;
		if (textureId == 0) glGenTextures(1, &textureId);  				// id generation
		glBindTexture(GL_TEXTURE_2D, textureId);    // binding
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_FLOAT, &image[0]); // To GPU
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, sampling); // sampling
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, sampling);
	}

	float getMinDistanceFromObjects(const glm::vec3& ray, const glm::vec3& dir = glm::vec3()) {
		++iter;
		std::vector<float> distances;
		/*float scale = 2.0f;
		glm::vec3 displacement = glm::vec3(0.0f, 0.0f, 0.0f);
		glm::vec3 fromObjectSpace = ray - displacement;*/
		//mat3 rotationMatrix = rotation(rotationAngles);
		//rotationMatrix = transpose(rotationMatrix); // Orthonormal matrix ---> A^T*A = I
		//fromObjectSpace = fromObjectSpace * rotationMatrix; // apply inverse rotation R^-1 
		//fromObjectSpace *= 1 / scale; //apply scaling
		glm::vec3 csavart = ray * glm::inverse(csavaras(ray, ray.z * 0.3f));
		glm::vec3 sphereRay;


		//TODO: sugár metszi-e?
		//bounding volume hierarchy
		if (!deformed)
			sphereRay = ray;
		else if (boundingBox->intersectBox(ray))
			sphereRay = ray * glm::inverse(freeForm(glm::vec3(ray.x, ray.y / 2.0f, ray.z)));

			
		distances.push_back( objects[0]->intersect(sphereRay));
		distances.push_back(objects[1]->intersect(ray));
		distances.push_back(objects[2]->intersect(ray));
		distances.push_back(objects[3]->intersect(ray));
		if(showControlPoints)
			for (auto iter = bezierControlPoints.begin(); iter != bezierControlPoints.end(); ++iter )
				distances.push_back(iter->get()->intersect(ray));


		float minDist = *std::min_element(distances.begin(), distances.end(),[](const auto& a, const auto& b) { return a < b; });

		//return objects[0];
		//float distortion = sin(3.0 * p.x) * sin(3.0 * p.y) * sin(3.0 * p.z) * 0.15;
		return minDist;// +distortion;
	}

	glm::vec3 calculateNormalFromGradient(const glm::vec3& p)
	{
		float delta = 0.001f;
		float f = getMinDistanceFromObjects(p);
		float gradient_x, gradient_y, gradient_z;
		#pragma omp parallel sections num_threads(3)
		{
#pragma omp section
			gradient_x = getMinDistanceFromObjects(p + glm::vec3(delta, 0, 0)) - f;
#pragma omp section
			gradient_y = getMinDistanceFromObjects(p + glm::vec3(0, delta, 0)) - f;
#pragma omp section
			gradient_z = getMinDistanceFromObjects(p + glm::vec3(0, 0, delta)) - f;
			
		}	
		return glm::normalize(glm::vec3(gradient_x, gradient_y, gradient_z));
	}

	bool sphereTraceShadow(glm::vec3 rayOrigin, glm::vec3 rayDirection, float maxDistance)
	{
		float threshold = 0.01f;
		float t = 0;
		while (t < maxDistance) {
			float minDistance = 10.0f;
			glm::vec3 from = rayOrigin + t * rayDirection;
			if (getMinDistanceFromObjects(from) < minDistance)
			{
				minDistance = getMinDistanceFromObjects(from);
				minDistance += threshold;
				if (minDistance < threshold * t)
					return true;
			}
			// no intersection, move along the ray by minDistance
			t += minDistance;
		}
		return false;
	}

	glm::vec3 shade(glm::vec3 rayOrigin, glm::vec3 rayDirection, float t)
	{
		glm::vec3 p = rayOrigin + t * rayDirection;
		glm::vec3 n = calculateNormalFromGradient(p);
		glm::vec3 R = glm::vec3(0.0f, 0.0f, 0.0f);

		glm::vec3 lightDir = glm::vec3(2.0, 5.0, 0.0) - p;
		if (dot(n, lightDir) > 0.0f) {
			float dist = length(glm::vec3(2.0, 5.0, -3.0) - p);
			lightDir = normalize(lightDir);
			bool shadow = sphereTraceShadow(p, lightDir, dist);
			if (shadow)
				R = R + (glm::vec3(1.0f, 0.0f, 0.0f) * dot(n, lightDir) * 0.5f);
		}

		return R;
	}

	glm::mat3 csavaras(glm::vec3 p, float angle) {
		float fzconst = 0.3f;
		return glm::mat3(cos(angle), -sin(angle), -p.x * sin(angle) * fzconst - p.y * cos(angle) * fzconst,
			sin(angle), cos(angle), p.x * cos(angle) * fzconst - p.y * sin(angle) * fzconst,
			0, 0, 1);
	}

	int fact(int n) {
		int sum = 1;
		for (int i = 1; i <= n; ++i) {
			sum = sum * i;
		}
		return sum;
	}

	float get_b(int d, int i, float t) {
		float combination = static_cast<float>(fact(d) / (fact(i) * fact(d - i)));
		return static_cast<float>(combination * pow(1.0f - t, d - i) * pow(t, i));
	}

	glm::vec3 bezier(int n, int m, int l, glm::vec3 p) {
		glm::vec3 sumVector = glm::vec3(0.0f, 0.0f, 0.0f);
		float bx = 0.0f;
		float by = 0.0f;
		float bz = 0.0f;

		for (int i = 0; i <= n; ++i)
		{
			for (int j = 0; j <= m; ++j)
			{
				for (int k = 0; k <= l; ++k)
				{
					bx = get_b(n, i, p.x);
					by = get_b(m, j, p.y);
					bz = get_b(l, k, p.z);
					if (n == controlPointNumberX-2) {
						sumVector.x += static_cast<float>(controlPointNumberX - 1) * (points[i + 1][j][k].x - points[i][j][k].x) * bx * by * bz;
						sumVector.y += static_cast<float>(controlPointNumberX - 1) * (points[i + 1][j][k].y - points[i][j][k].y) * bx * by * bz;
						sumVector.z += static_cast<float>(controlPointNumberX - 1) * (points[i + 1][j][k].z - points[i][j][k].z) * bx * by * bz;
					}
					else if (m == controlPointNumberY - 2) {
						sumVector.x += static_cast<float>(controlPointNumberY - 1) * (points[i][j + 1][k].x - points[i][j][k].x) * bx * by * bz;
						sumVector.y += static_cast<float>(controlPointNumberY - 1) * (points[i][j + 1][k].y - points[i][j][k].y) * bx * by * bz;
						sumVector.z += static_cast<float>(controlPointNumberY - 1) * (points[i][j + 1][k].z - points[i][j][k].z) * bx * by * bz;
					}
					else if (l == controlPointNumberZ - 2) {
						sumVector.x += static_cast<float>(controlPointNumberZ - 1) * (points[i][j][k + 1].x - points[i][j][k].x) * bx * by * bz;
						sumVector.y += static_cast<float>(controlPointNumberZ - 1) * (points[i][j][k + 1].y - points[i][j][k].y) * bx * by * bz;
						sumVector.z += static_cast<float>(controlPointNumberZ - 1) * (points[i][j][k + 1].z - points[i][j][k].z) * bx * by * bz;
					}
				}
			}
		}
		return sumVector;
	}

	glm::mat3 freeForm(glm::vec3 p) {
		glm::mat3 matrix;
		matrix[0] = bezier(controlPointNumberX -2 , controlPointNumberY -1 , controlPointNumberZ -1, p);
		matrix[1] = bezier(controlPointNumberX - 1, controlPointNumberY - 2, controlPointNumberZ - 1, p);
		matrix[2] = bezier(controlPointNumberX - 1, controlPointNumberY - 1, controlPointNumberZ - 2, p);
		return matrix;

	}

	glm::vec3 rayTrace(Ray ray) {
		float total_distance_traveled = 0.0;
		for (int i = 0; i < 256; ++i)
		{

			//glm::vec3 csavartsugar = ray.start * glm::inverse(csavaras(ray.start, ray.start.z * 0.3f));
			glm::vec3 current_position = ray.start + (total_distance_traveled * ray.dir);
			//current_position += total_distance_traveled *ray.dir * glm::inverse(csavaras(current_position, current_position.z * 0.3f));

			float distance_to_closest = getMinDistanceFromObjects(current_position);

			if (distance_to_closest < 0.001f)
			{
				glm::vec3 normal = calculateNormalFromGradient(current_position);
				glm::vec3 light_position = glm::vec3(10.0, 10.0, -10.0);
				glm::vec3 direction_to_light = normalize(light_position - current_position);
				float diffuse_intensity = (0.2f > dot(normal, direction_to_light)) ? 0.2f : dot(normal, direction_to_light);
				//glm::vec3 shadeing = shade(ray.start, ray.dir, total_distance_traveled);
				//return normal * 0.5 + 0.5;
				return glm::vec3(1.0, 0.0, 0.0) * diffuse_intensity;// -shadeing;
				//return glm::vec3(1.0, 0.0, 0.0);// *diffuse_intensity;
			}
			if (total_distance_traveled > 50.0)
			{
				break;
			}
			total_distance_traveled += distance_to_closest;
		}
		//std::cout << i << currentposition <<std::endl;
		return glm::vec3(0.0f,0.0f,0.0f);
	}

	~Scene() {
		if (textureId > 0) glDeleteTextures(1, &textureId);
	}
};

class FullScreenTexturedQuad {
	unsigned int vao;	// vertex array object id and texture id
	Scene scene;
public:
	FullScreenTexturedQuad(int windowWidth, int windowHeight, std::vector<glm::vec4>& image)
		: scene(windowWidth, windowHeight, image)
	{
		glGenVertexArrays(1, &vao);	// create 1 vertex array object
		glBindVertexArray(vao);		// make it active

		unsigned int vbo;		// vertex buffer objects
		glGenBuffers(1, &vbo);	// Generate 1 vertex buffer objects

		// vertex coordinates: vbo0 -> Attrib Array 0 -> vertexPosition of the vertex shader
		glBindBuffer(GL_ARRAY_BUFFER, vbo); // make it active, it is an array
		float vertexCoords[] = { -1, -1,  1, -1,  1, 1,  -1, 1 };	// two triangles forming a quad
		glBufferData(GL_ARRAY_BUFFER, sizeof(vertexCoords), vertexCoords, GL_STATIC_DRAW);	   // copy to that part of the memory which is not modified 
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, NULL);     // stride and offset: it is tightly packed
	}

	void Draw(unsigned int program) {
		scene.render();
		glBindVertexArray(vao);	// make the vao and its vbos active playing the role of the data source
		int location = glGetUniformLocation(program, "textureUnit");
		if (location >= 0) {
			glUniform1i(location, 0);
			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, scene.textureId);
		}
		
		glDrawArrays(GL_TRIANGLE_FAN, 0, 4);	// draw two triangles forming a quad
	}
	Scene& getScene() {
		return scene;
	}
};

int main()
{
	GLFWwindow* window;
	FullScreenTexturedQuad* fullScreenTexturedQuad;
	/* Initialize the library */
	if (!glfwInit())
		return -1;




	/* Create a windowed mode window and its OpenGL context */
	window = glfwCreateWindow(size, size, "Sphere Tracing", NULL, NULL);
	if (!window)
	{
		glfwTerminate();
		return -1;
	}

	/* Make the window's context current */
	glfwMakeContextCurrent(window);

	IMGUI_CHECKVERSION();
	ImGui::CreateContext();

	// Setup Dear ImGui style
	ImGui::StyleColorsDark();
	//ImGui::StyleColorsClassic();
	const char* glsl_version = "#version 330";
	// Setup Platform/Renderer backends
	ImGui_ImplGlfw_InitForOpenGL(window, true);
	ImGui_ImplOpenGL3_Init(glsl_version);

	//glfwSwapInterval(1);
	if (glewInit() != GLEW_OK)
		std::cout << "Baj van";
	int majorVersion, minorVersion;
	printf("GL Vendor    : %s\n", glGetString(GL_VENDOR));
	printf("GL Renderer  : %s\n", glGetString(GL_RENDERER));
	printf("GL Version (string)  : %s\n", glGetString(GL_VERSION));
	glGetIntegerv(GL_MAJOR_VERSION, &majorVersion);
	glGetIntegerv(GL_MINOR_VERSION, &minorVersion);
	printf("GL Version (integer) : %d.%d\n", majorVersion, minorVersion);
	printf("GLSL Version : %s\n", glGetString(GL_SHADING_LANGUAGE_VERSION));


	ShaderProgramSource source = ParseShader("Shader.shader");
	unsigned int program = CreateProgram(source.VertexSource, source.FragmentSource);
	glUseProgram(program);
	std::vector<glm::vec4> image(size * size);
	fullScreenTexturedQuad = new FullScreenTexturedQuad(size, size, image);

	bool show_demo_window = true;
	bool show_another_window = false;
	ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);

	while (!glfwWindowShouldClose(window)) {
		glfwPollEvents();
		if (GLFW_PRESS == glfwGetKey(window, GLFW_KEY_ESCAPE)) {
			glfwSetWindowShouldClose(window, 1);
		}
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();


		fullScreenTexturedQuad->Draw(program);

		// 2. Show a simple window that we create ourselves. We use a Begin/End pair to created a named window.
		{
			static float f = 0.0f;
			static int counter = 0;
			Camera& camera = fullScreenTexturedQuad->getScene().camera;
			ImGui::Begin("Menu");
			ImGui::Checkbox("Deform", &fullScreenTexturedQuad->getScene().deformed);     
			ImGui::Checkbox("Show control points", &fullScreenTexturedQuad->getScene().showControlPoints);
			float cameraPosition[3] = { camera.eye.x, camera.eye.y, camera.eye.z };
			ImGui::SliderFloat3("Camera Position",cameraPosition, -8.0f, 8.0f);            // Edit 1 float using a slider from 0.0f to 1.0f
			camera.eye.x = cameraPosition[0];
			camera.eye.y = cameraPosition[1];
			camera.eye.z = cameraPosition[2];
			int arr[3] = { fullScreenTexturedQuad->getScene().controlPointNumberX,fullScreenTexturedQuad->getScene().controlPointNumberY ,fullScreenTexturedQuad->getScene().controlPointNumberZ };
			ImGui::SliderInt3("Number of control points on each axis", arr, 2, 4);
			fullScreenTexturedQuad->getScene().setControlPoints(arr[0], arr[1], arr[2]);

			//for (unsigned int i = 0; i < fullScreenTexturedQuad->getScene().bezierControlPoints.size(); ++i)
			//{
			//	std::string label = "Control Point " + std::to_string(i + 1);
			//	ImGui::PushID(i);
			//	glm::vec3 tempPoint = fullScreenTexturedQuad->getScene().bezierControlPoints[i]->getCenter();
			//	float temp[3] = { tempPoint.x, tempPoint.y, tempPoint.z };
			//	ImGui::SliderFloat3(label.c_str(), temp, -3.0f, 3.0f);
			//	fullScreenTexturedQuad->getScene().bezierControlPoints[i]->setCenter(glm::vec3(temp[0], temp[1], temp[2]));
			//	//todo: write back to point
			//	ImGui::PopID();
			//}
			size_t x = fullScreenTexturedQuad->getScene().controlPointNumberX;
			size_t y = fullScreenTexturedQuad->getScene().controlPointNumberY;
			size_t z = fullScreenTexturedQuad->getScene().controlPointNumberZ;

			for (size_t i = 0; i < x; i++)
			{
				for (size_t j = 0; j < y; j++)
				{
					for (size_t k = 0; k < z; k++)
					{
						std::string label = "Control Point " + std::to_string(i * x * y + j * z + k  + 1);
						ImGui::PushID(i);
						glm::vec3 tempPoint = fullScreenTexturedQuad->getScene().points[i][j][k];
						float temp[3] = { tempPoint.x, tempPoint.y, tempPoint.z };
						ImGui::SliderFloat3(label.c_str(), temp, -3.0f, 3.0f);
						fullScreenTexturedQuad->getScene().points[i][j][k] = glm::vec3(temp[0], temp[1], temp[2]);
						fullScreenTexturedQuad->getScene().bezierControlPoints[i*x*y+j*z+k]->setCenter(glm::vec3(temp[0], temp[1], temp[2]));
						//todo: write back to point
						ImGui::PopID();
					}
				}
			}


			camera.set(camera.eye, camera.lookat, camera.up, camera.fov);

			if (ImGui::Button("Render"))                            // Buttons return true when clicked (most widgets return true when edited/activated)
				fullScreenTexturedQuad->getScene().build();

			ImGui::Text("Rendering may take some time\nApplication average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
			ImGui::End();
		}


		ImGui::Render();
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
		glfwSwapBuffers(window);
		
	}
	/*std::cout << iter << std::endl;*/
	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();
	glfwDestroyWindow(window);
	glfwTerminate();
	return 0;
}
