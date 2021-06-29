#include "gltf_loader.h"

#include <common/ml_document/mesh_model.h>
#include <common/mlexception.h>

//declarations
enum GLTF_ATTR_TYPE {POSITION, NORMAL, TEXCOORD_0, INDICES};
const std::array<std::string, 3> GLTF_ATTR_STR {"POSITION", "NORMAL", "TEXCOORD_0"};

void loadMeshPrimitive(
		MeshModel& m,
		const tinygltf::Model& model,
		const tinygltf::Primitive& p);

void loadAttribute(
		MeshModel& m,
		std::vector<CMeshO::VertexPointer>& ivp,
		const tinygltf::Model& model,
		const tinygltf::Primitive& p,
		GLTF_ATTR_TYPE attr,
		int textID = -1);

template <typename Scalar>
void populateAttr(
		GLTF_ATTR_TYPE attr,
		MeshModel&m,
		std::vector<CMeshO::VertexPointer>& ivp,
		const Scalar* array,
		unsigned int number,
		int textID = -1);

template <typename Scalar>
void populateVertices(
		MeshModel& m,
		std::vector<CMeshO::VertexPointer>& ivp,
		const Scalar* posArray,
		unsigned int vertNumber);

template <typename Scalar>
void populateVNormals(
		const std::vector<CMeshO::VertexPointer>& ivp,
		const Scalar* normArray,
		unsigned int vertNumber);

template <typename Scalar>
void populateVTextCoords(
		const std::vector<CMeshO::VertexPointer>& ivp,
		const Scalar* textCoordArray,
		unsigned int vertNumber,
		int textID);

template <typename Scalar>
void populateTriangles(
		MeshModel&m,
		const std::vector<CMeshO::VertexPointer>& ivp,
		const Scalar* triArray,
		unsigned int triNumber);

 /**************
 * Definitions *
 **************/

void gltf::loadMesh(
		MeshModel& m,
		const tinygltf::Mesh& tm,
		const tinygltf::Model& model)
{
	if (!tm.name.empty())
		m.setLabel(QString::fromStdString(tm.name));
	for (const tinygltf::Primitive& p : tm.primitives){
		loadMeshPrimitive(m, model, p);
	}
}

void loadMeshPrimitive(
		MeshModel& m,
		const tinygltf::Model& model,
		const tinygltf::Primitive& p)
{
	int textureImg = -1;
	if (p.material >= 0) {
		const tinygltf::Material& mat = model.materials[p.material];
		auto it = mat.values.find("baseColorTexture");
		if (it != mat.values.end()){ //the material is a texture
			auto it2 = it->second.json_double_value.find("index");
			if (it2 != it->second.json_double_value.end()){
				textureImg = it->second.number_value;
			}
		}
	}
	if (textureImg != -1) {
		m.cm.textures.push_back(model.images[textureImg].uri);
		textureImg = m.cm.textures.size()-1;
	}
	std::vector<CMeshO::VertexPointer> ivp;
	loadAttribute(m, ivp, model, p, POSITION);
	loadAttribute(m, ivp, model, p, NORMAL);
	loadAttribute(m, ivp, model, p, TEXCOORD_0, textureImg);
	loadAttribute(m, ivp, model, p, INDICES);
}

void loadAttribute(
		MeshModel& m,
		std::vector<CMeshO::VertexPointer>& ivp,
		const tinygltf::Model& model,
		const tinygltf::Primitive& p,
		GLTF_ATTR_TYPE attr,
		int textID)
{
	const tinygltf::Accessor* accessor = nullptr;

	if (attr != INDICES) {
		auto it = p.attributes.find(GLTF_ATTR_STR[attr]);

		if (it != p.attributes.end()) {
			accessor = &model.accessors[it->second];
		}
		else if (attr == POSITION) {
			throw MLException("File has not 'Position' attribute");
		}
	}
	else {
		if (p.mode == 4 && p.indices >= 0 &&
				(unsigned int)p.indices < model.accessors.size()) {
			accessor = &model.accessors[p.indices];
		}
	}

	if (accessor) {
		const tinygltf::BufferView& posbw = model.bufferViews[accessor->bufferView];
		const std::vector<unsigned char>& posdata = model.buffers[posbw.buffer].data;
		unsigned int posOffset = posbw.byteOffset + accessor->byteOffset;


		if (accessor->componentType == TINYGLTF_COMPONENT_TYPE_FLOAT) {
			const float* posArray = (const float*) (posdata.data() + posOffset);
			populateAttr(attr, m, ivp, posArray, accessor->count, textID);
		}
		else if (accessor->componentType == TINYGLTF_COMPONENT_TYPE_DOUBLE) {
			const double* posArray = (const double*) (posdata.data() + posOffset);
			populateAttr(attr, m, ivp, posArray, accessor->count, textID);
		}
		else if (accessor->componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_BYTE) {
			const unsigned char* triArray = (const unsigned char*) (posdata.data() + posOffset);
			populateAttr(attr, m, ivp, triArray, accessor->count, textID);
		}
		else if (accessor->componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT) {
			const unsigned short* triArray = (const unsigned short*) (posdata.data() + posOffset);
			populateAttr(attr, m, ivp, triArray, accessor->count, textID);
		}
		else if (accessor->componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT) {
			const unsigned int* triArray = (const unsigned int*) (posdata.data() + posOffset);
			populateAttr(attr, m, ivp, triArray, accessor->count, textID);
		}
	}
	else if (attr == INDICES) {
		populateAttr<unsigned char>(attr, m, ivp, nullptr, 0);
	}

}

template <typename Scalar>
void populateAttr(
		GLTF_ATTR_TYPE attr,
		MeshModel&m,
		std::vector<CMeshO::VertexPointer>& ivp,
		const Scalar* array,
		unsigned int number,
		int textID)
{
	switch (attr) {
	case POSITION:
		populateVertices(m, ivp, array, number); break;
	case NORMAL:
		populateVNormals(ivp, array, number); break;
	case TEXCOORD_0:
		m.enable(vcg::tri::io::Mask::IOM_VERTTEXCOORD);
		populateVTextCoords(ivp, array, number, textID); break;
		break;
	case INDICES:
		populateTriangles(m, ivp, array, number/3); break;
	}
}

template <typename Scalar>
void populateVertices(
		MeshModel&m,
		std::vector<CMeshO::VertexPointer>& ivp,
		const Scalar* posArray,
		unsigned int vertNumber)
{
	ivp.clear();
	ivp.resize(vertNumber);
	CMeshO::VertexIterator vi =
			vcg::tri::Allocator<CMeshO>::AddVertices(m.cm, vertNumber);
	for (unsigned int i = 0; i < vertNumber*3; i+= 3, ++vi){
		ivp[i/3] = &*vi;
		vi->P() = CMeshO::CoordType(posArray[i], posArray[i+1], posArray[i+2]);
	}
}

template <typename Scalar>
void populateVNormals(
		const std::vector<CMeshO::VertexPointer>& ivp,
		const Scalar* normArray,
		unsigned int vertNumber)
{
	for (unsigned int i = 0; i < vertNumber*3; i+= 3){
		ivp[i/3]->N() = CMeshO::CoordType(normArray[i], normArray[i+1], normArray[i+2]);
	}
}

template <typename Scalar>
void populateVTextCoords(
		const std::vector<CMeshO::VertexPointer>& ivp,
		const Scalar* textCoordArray,
		unsigned int vertNumber,
		int textID)
{
	for (unsigned int i = 0; i < vertNumber*2; i+= 2) {
		ivp[i/2]->T() = CMeshO::VertexType::TexCoordType(textCoordArray[i], 1-textCoordArray[i+1]);
		ivp[i/2]->T().N() = textID;
	}
}

template <typename Scalar>
void populateTriangles(
		MeshModel&m,
		const std::vector<CMeshO::VertexPointer>& ivp,
		const Scalar* triArray,
		unsigned int triNumber)
{
	if (triArray != nullptr) {
		CMeshO::FaceIterator fi =
				vcg::tri::Allocator<CMeshO>::AddFaces(m.cm, triNumber);
		for (unsigned int i = 0; i < triNumber*3; i+=3, ++fi) {
			fi->V(0) = ivp[triArray[i]];
			fi->V(1) = ivp[triArray[i+1]];
			fi->V(2) = ivp[triArray[i+2]];
		}
	}
	else {
		CMeshO::FaceIterator fi =
				vcg::tri::Allocator<CMeshO>::AddFaces(m.cm, ivp.size()/3);
		for (unsigned int i = 0; i < ivp.size(); i+=3, ++fi) {
			fi->V(0) = ivp[i];
			fi->V(1) = ivp[i+1];
			fi->V(2) = ivp[i+2];
		}
	}
}


